import numpy as np


def average_ratio(data_1, data_2, calib_factor = None):
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

    assert data_1.shape == data_2.shape and len(data_1.shape) == 1 and len(data_2.shape) == 1

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




class zScanDataAnalyser:
    def __init__(self, oa_signal, ca_signal, ref_signal, total_steps):
        """ Prior to calling this initialisator, the probe has to be inserted into the beam but
            far displaced from the focal spot, i.e. nonlinear effects are highly unlikely.
            Additionally, the aperture has to be opened.

            Input:
            -------
            oa, ca, ref_signal: 1-dim numpy arrays of Nidaq measurements, all of identical length.
            total_steps: integer, corresponding to the total number of stage positions at which the
                         photodiode signals will be measured.
        """
        self.c_OA = self.extract_calibration_factor(oa_signal, ref_signal)
        self.c_CA = self.extract_calibration_factor(ca_signal, ref_signal)
        self.S = 1  # Transmission of aperture
        self.combined_c_CA = self.c_CA * self.S

        self.T_OA = np.zeros(shape=(total_steps, 2))  # transmission in open aperture path.
        self.T_CA = np.zeros(shape=(total_steps, 2))  # transmission in closed aperture path.

        self.total_steps = total_steps     # The total number of measurement stage positions.
        self.current_position_step = 0     # Integer indicating next empty transmission array entry.



    def extract_calibration_factor(self, pd_signal, ref_signal):
        """ If the probe is placed into the beam far displaced from the focal spot, i.e. nonlinear
            effects are highly unlikely, we consider this situation as 100% transmission. As such,
            we want to calibrate the photodiode signals to be "equal".
            Hence, for this function to produce a reasonable calibration factor, place the probe
            into the beam but far displaced from the from the focal spot. Additionally, open the
            aperture to 100% transmission (S=1).
            
            Input
            ------
            pd_signal: 1-dim numpy array. Nidaq measurements of either the closed or open aperture
            photodiode.
            ref_signal: 1-dim numpy array. Nidaq measurements of the reference photodiode.

            Output
            -------
            calibration_factor: 1-dim numpy array. First entry being the average calibration factor,
            second entry being its corresponding error.
        """
        return average_ratio(pd_signal, ref_signal)


    def extract_aperture_transmission(self, ca_signal, ref_signal):
        """ If the probe is placed into the beam but is far displaced from the focal spot, i.e.
            nonlinear effects are highly unlikely, and the calibration factor self.c_CA between
            closed aperture photodiode and reference photodiode has already been computed,
            the aperture can be closed to a transmission smaller than unity (self.S<1). This
            function takes the closed aperture photodiode signal and the the reference photodiode
            signal and evaluates the aperture transmission self.S. It stores this value in the
            instance property self.S and returns it in addition.
        """
        self.S = average_ratio(ca_signal, ref_signal, calib_factor=self.c_CA)
        self.combined_c_CA[0] = self.c_CA[0] * self.S[0]
        self.combined_c_CA[1] = np.sqrt((self.c_CA[0]*self.S[1])**2 + (self.c_CA[1]*self.S[0])**2)
        return self.S


    def extract_oa_ca_transmissions(self, position, oa_signal, ca_signal, ref_signal):
        """ At any one probe position, this function evaluates the transmission in both the open and
            closed aperture path. The result is written into the arrays self.T_CA and self.T_OA in
            a 1-dim numpy array [position, transmission, error on transmission].
        """
        transmission_oa = average_ratio(oa_signal, ref_signal, self.c_OA)
        transmission_ca = average_ratio(ca_signal, ref_signal, self.combined_c_CA)
        self.T_OA[current_position_step] = np.array([position, *transmission_oa])
        self.T_CA[current_position_step] = np.array([position, *transmission_ca])

        self.current_position_step += 1
        assert self.current_position_step <= total_steps




if __name__ == '__main__':
    a = np.array([.5,.5,.5,.50005,.5])
    b = np.array([1,1,1,1,1])

    calib = average_ratio(a,b)
    print(calib)

    print(average_ratio(a,b,calib))



