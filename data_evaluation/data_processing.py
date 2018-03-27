import numpy as np
import os



# Constants:

# Calibration of reference photodiode signal into pulse energy in J/V (due to factor of 1e-6)
# Factor of 1/100 because the calibration factor was measured with an OD2 neutral density filter in
# front of the reference photodiode.
CONSTANTS_calib_photodiode_pulse_energy = np.array([2.8749, 0.0030]) * 1e-6 / 100





def average_ratio_element_wise(data_1, data_2, calib_factor = None):
    """ Computes the statistical average of the element wise ratio of the two arrays data_1 and
        data_2. For this calculation, each element wise ratio value is assumed to be of the same
        significance such that the standard deviation of the average is considered to be the
        standard deviation of the set of all element wise ratio values.
        If a calibration factor is provided, the calibrated average is the uncalibrated average 
        multiplied with the calibration value.

        Input Parameters:
        ------------
        data_1: 1-dim numpy array of single measurements
        data_2: 1-dim numpy array of single measurements, same length as data_1
        calib_factor: Optional, if provided it must be a 1-dim array with two entries. The first
                      being the calibration factor, the second its error.

        Output:
        ------------
        1-dim array of length two. The first entry being the average of the element wise ratio of
        data_1 and data_2, the second being the corresponding standard deviation.
    """

    assert data_1.shape == data_2.shape
    assert len(data_1.shape) == 1
    assert len(data_2.shape) == 1

    ratios = data_1 / data_2
    
    mean = np.mean(ratios)
    std = np.std(ratios, ddof=1)

    if calib_factor is not None:
        assert calib_factor.shape == (2,)
        calibrated_mean = mean / calib_factor[0]
        calibrated_std = np.sqrt(
            (std / calib_factor[0])**2 + (mean / calib_factor[0]**2 * calib_factor[1])**2)
        return np.array([calibrated_mean, calibrated_std])

    return np.array([mean, std])



def average_ratio_array_wise(data_1, data_2, calib_factor = None):
    """ Computes the ratio of the mean of data_1 and data_2 and considers the standard deviation.
        This function basically is the same as average_ratio_element_wise in its behaviour. See
        the corresponding docstring. This function makes more sense if data_1 and data_2 are
        temporarily uncorrelated.
    """

    assert data_1.shape == data_2.shape
    assert len(data_1.shape) == 1
    assert len(data_2.shape) == 1

    data1_mean, data1_std = data_1.mean(), data_1.std()
    data2_mean, data2_std = data_2.mean(), data_2.std()

    ratio = data1_mean / data2_mean
    std = np.sqrt( (data1_std/data2_mean)**2 + (data1_mean*data2_std/data2_mean**2)**2 )

    if calib_factor is not None:
        assert calib_factor.shape == (2,)
        calibrated_mean = ratio / calib_factor[0]
        calibrated_std = np.sqrt(
            (std / calib_factor[0])**2 + (ratio / calib_factor[0]**2 * calib_factor[1])**2)
        return np.array([calibrated_mean, calibrated_std])

    return np.array([ratio, std])





class zScanDataProcessor:

    def __init__(self):
        """ Class to extract calibration factors and transmissions from the photodiode signals.
            Prior to invoking the function "extract_oa_ca_transmissions", these functions need be
            called:
                1. extract_calibration_factors
                2. extract_aperture_transmission (in this order)
            Only then is it sensible to call the function "extract_oa_ca_transmissions".
        """
        self.c_OA = None
        self.c_CA = None
        self.S = None                         # Transmission of aperture
        self.combined_c_CA = None             # S*c_CA

        self.T_OA = None                   # transmission data of open aperture path.
        self.T_CA = None                   # transmission data of closed aperture path.
        self.current_position_step = 0     # Integer indicating next empty transmission array entry.



    def extract_calibration_factors(self, ref_signal, oa_signal, ca_signal):
        """ If the sample is placed into the beam far displaced from the focal spot, i.e. nonlinear
            effects are highly unlikely, we consider this situation as 100% transmission. As such,
            we want to calibrate the photodiode signals to be "equal".
            Hence, for this function to produce reasonable calibration factors, place the probe
            into the beam but far displaced from the from the focal spot. Additionally, open the
            aperture to 100% transmission (S=1).
            
            Input
            ------
            oa_signal: 1-dim numpy array. Nidaq measurements of the open aperture photodiode.
            ca_signal: 1-dim numpy array. Nidaq measurements of the closed aperture photodiode.
            ref_signal: 1-dim numpy array. Nidaq measurements of the reference photodiode.

            Output
            -------
            1-dim array with two entries. First entry is the open, second entry the closed aperture
            calibration factor. Each entry itself is a 1-dim numpy array with the first entry being
            the average calibration factor and the second entry being its corresponding error.
        """
        self.c_OA = average_ratio_element_wise(oa_signal, ref_signal)
        self.c_CA = average_ratio_element_wise(ca_signal, ref_signal)
        self.S = 1
        self.combined_c_CA = self.c_CA * self.S
        return np.array([self.c_OA, self.c_CA])


    def extract_aperture_transmission(self, ref_signal, ca_signal):
        """ If the probe is placed into the beam but is far displaced from the focal spot, i.e.
            nonlinear effects are highly unlikely, and the calibration factor self.c_CA between
            closed aperture photodiode and reference photodiode has already been computed,
            the aperture can be closed to a transmission smaller than unity (self.S<1). This
            function takes the closed aperture photodiode signal and the the reference photodiode
            signal and evaluates the aperture transmission self.S. It stores this value in the
            instance property self.S and returns it in addition.
        """
        assert self.c_CA is not None
        assert self.c_OA is not None
        self.S = average_ratio_element_wise(ca_signal, ref_signal, calib_factor=self.c_CA)
        self.combined_c_CA[0] = self.c_CA[0] * self.S[0]
        self.combined_c_CA[1] = np.sqrt((self.c_CA[0]*self.S[1])**2 + (self.c_CA[1]*self.S[0])**2)
        return self.S


    def extract_oa_ca_transmissions(self, position, ref_signal, oa_signal, ca_signal,
        tot_num_of_pos):
        """ At any one probe position, this function evaluates the transmission in both the open and
            closed aperture path. The result is written into the arrays self.T_CA and self.T_OA in
            a 1-dim numpy array [position, transmission, error on transmission].

            Input:
            tot_num_of_pos  Integer denoting the total number of positions at which the photodiode
                            signals will be measured.
        """
        assert self.c_OA is not None
        assert self.c_CA is not None
        assert self.S is not None
        assert self.combined_c_CA is not None
        assert self.current_position_step < tot_num_of_pos
        
        if self.current_position_step == 0:
            assert self.T_OA is None
            assert self.T_CA is None
            self.T_OA = np.zeros(shape=(tot_num_of_pos, 3))
            self.T_CA = np.zeros(shape=(tot_num_of_pos, 3))

        transmission_oa = average_ratio_element_wise(oa_signal, ref_signal, self.c_OA)
        transmission_ca = average_ratio_element_wise(ca_signal, ref_signal, self.combined_c_CA)
        self.T_OA[self.current_position_step] = np.array([position, *transmission_oa])
        self.T_CA[self.current_position_step] = np.array([position, *transmission_ca])

        self.current_position_step += 1


    def store_transmission_data(self, storage_dir, folder_num, data_file_header):
        """ Stores the transmission data self.T_CA and self.T_OA into a file.

            Input:
            storage_dir         String specifying the directory to store the transmission data file.
            folder_num          Integer specifiying the number suffix of the folder.
            data_file_header    String to be written into the header of the transmission data file.
        """

        # assert that position entries in T_OA and T_CA are identical
        assert (self.T_OA[:,0] == self.T_CA[:,0]).all()
        transmission_array = np.concatenate((self.T_OA, self.T_CA[:,1:]), axis=1)

        try:
            file = os.path.join(storage_dir, "transmission_data_{0:02}.dat".format(folder_num))
            np.savetxt(file, transmission_array, header=data_file_header, fmt="%10.4f")
        except Exception as ex:
            print("Storage of transmission data failed!!!!")
            traceback.print_exc()
            print("\n")


    def reinitialise(self):
        self.T_OA = None
        self.T_CA = None
        self.current_position_step = 0


    def extract_pulse_energy(self, ref_signal_max):
        """ This function extracts from the reference photodiode signal the pulse energy.
            For this to work, I have performed a fit of pulse energy (measured with an energymeter)
            versus reference photodiode signal (measured with PCI Nidaq) at fixed 30Hz laser
            repetition rate.
            Note that there is also a dependence between laser repitition rate, pulse energy and
            photodiode signal which is yet to be measured and included into this function.

            Input:
            ref_signal_max: 1-dim numpy array of arbitrary length > 1. The maximum measured values
            of the reference photodiode in Volt.

            Output:
            1-dim numpy array of length 2, the first entry denoting the pulse_energy, the second
            its error. Both in µJ.

            Note:
            For the calibration measurement, see folder "Masterarbeit/Z-Scan/Setup_calibration/
            Calibration_PulseEnergy_Photodiode/" stored on my computer.
        """
        def pulse_energy(pd_ref_signal):
            """ Fit function. Input be a 1-dim numpy array with two entries, the first entry being
                the signal, the second its error. Both in Volt.
                
                Returns:
                1-dim numpy array of length 2, the first entry denoting the pulse_energy, the second
                its error. Both in J.
            """
            fit_gradient = CONSTANTS_calib_photodiode_pulse_energy

            pulse_e = fit_gradient[0] * pd_ref_signal[0]
            delta_pulse_e = np.sqrt( (fit_gradient[1] * pd_ref_signal[0])**2 + \
                                     (fit_gradient[0] * pd_ref_signal[1])**2 )

            return np.array([pulse_e, delta_pulse_e])


        mean = ref_signal_max.mean()
        std = ref_signal_max.std(ddof=1)

        return pulse_energy(np.array([mean, std]))


    def compute_effective_pulse_energy(self, pulse_energy, refr_index_sample, refr_index_ambient):
        """ Computes the pulse energy that the sample is exposed to, i.e. it computes the fraction
            of incident pulse energy that is not reflected by glass cuvette or by the sample.

            Input:
            pulse_energy        1 dim numpy array of length 2, first entry being the pulse energy,
                                second entry its error.
            refr_index_sample   float
            refr_index_ambient  float

            Output:
            effective pulse energy  1 dim numpy array of length 2, first entry being the eff. pulse
                                    energy, second entry its error.
        """

        transm_ambientMaterial_air = 1 - (
            (refr_index_ambient-1) / (refr_index_ambient+1) )**2
        transm_sample_ambientMaterial = 1 - (
            (refr_index_sample-refr_index_ambient) / (refr_index_sample+refr_index_ambient))**2

        transmission = transm_ambientMaterial_air * transm_sample_ambientMaterial
        return pulse_energy * transmission


    def compute_effective_peak_intensity(self, eff_pulse_energy, pulse_length_fwhm, w0):
        """ Computes the effective peak intensity.
            Input:
            eff_pulse_energy    1 dim numpy array of length 2, first entry being the eff. pulse
                                energy, second entry its error, both in J.
            pulse_length_fwhm   float, the pulse length in seconds.
            w0                  float, beam waist at the focus in metre

            Output:
            effective peak intensity    1 dim numpy array of length 2, first entry being the peak
                                        intensity, second its error, both in units of W/m^2
            """

        pulse_length = pulse_length_fwhm     # seconds        
        eff_pulse_energy = eff_pulse_energy  # in J

        P0 = np.array([eff_pulse_energy[0], eff_pulse_energy[1]]) / pulse_length  # in W
        
        I0 = np.empty(shape=(2))
        I0[0] = 2*P0[0] / (np.pi*w0**2)     # in W/m^2
        I0[1] = 2*P0[1] / (np.pi*w0**2)     # in W/m^2

        return I0


    def compute_alpha(self, transmission, refr_index_sample, refr_index_ambient, geom_length):
        """ Computes the linear absorption coefficient alpha by taking Fresnel reflection losses
            at the boundaries of sample and ambient material and of ambient material and air into
            consideration.

            Input:
            transmission        float, transmission through the sample at low irradiance
            refr_index_sample   float, refractive index of sample
            refr_index_ambient  float, refractive index of the ambient material of the sample
            geom_length         float, geometrical length of the sample in m

            Output:
            alpha   float, in 1/m
        """

        # transmission at each boundary:
        transmission_glass_air = 1 - ((refr_index_ambient-1)/(refr_index_ambient+1))**2
        transmission_sample_ambient = 1 - (
            (refr_index_ambient-refr_index_sample)/(refr_index_ambient+refr_index_sample))**2
        
        # each transmission is squared because the light has to transmit each boundary twice:
        expected_transmission = transmission_sample_ambient**2 * transmission_glass_air**2

        # the transmission through the medium (after being corrected by Fresnel transmission)
        transmission_medium = transmission / expected_transmission

        if transmission_medium >= 1:
            print("alpha would be negative. It is being set to zero manually.")
            return 0

        alpha = -np.log(transmission_medium) / geom_length  # in 1/m
        return alpha


    def compute_effective_length(self, geom_length, alpha):
        """ Computes the effective sample length.
            Input:
            geom_length     float, geometrical sample length in m
            alpha           float, linear absoprtion coefficient in 1/m

            Output:
            effective sample length     float, in m
        """
        if alpha == 0.:
            return geom_length
        else:
            return (1-np.exp(-alpha * geom_length)) / alpha

