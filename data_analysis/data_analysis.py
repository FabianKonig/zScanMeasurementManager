import numpy as np
import datetime
import os
import traceback
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt



CONSTANTS_beam_waist = 19.0537e-6  # m waist of incident beam in vacuum
CONSTANTS_wavelength = 532e-9      # m in vacuum

# Calibration of reference photodiode signal into pulse energy in J/V
CONSTANTS_calib_photodiode_pulse_energy = np.array([16.9274, 0.0212]) * 1e-6 

# Initial guess for fit parameters. The first entry denotes the beam waist position in mm
CONSTANTS_guess_OA = [22,1]  # second entry: dΨ
CONSTANTS_guess_CA = [22,1]  # second entry: dΦ

CONSTANTS_tolerance_z0 = 1.5  # absolute tolerance of beam waist position from z=22mm in mm

CONSTANTS_cuvette_reflectivity = 0.0377  # reflectivity of the glass/air boundary of the cuvette
CONSTANTS_pulse_length_FWHM = 15e-12  # laser pulse length in seconds

# I don't know why, but the fits are much better if I introduce a correction factor for the
# Rayleigh length in vacuum:
CONSTANTS_rayleigh_length_correction_factor = 1.7




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



def T_OA_func(z, z0, zR, dΨ):
    x = (z-z0)/zR
    return 1 - 2*(x**2+3)*dΨ / ((x**2+9) * (x**2+1))

def T_CA_func(z, z0, zR, dΦ, dΨ):
    x = (z-z0)/zR
    return T_OA_func(z, z0, zR, dΨ) + 4*x*dΦ / ((x**2+9) * (x**2+1))




class zScanDataAnalyser:

    def __init__(self, tot_num_of_pos, sample_material, solvent, concentration, laser_rep_rate,
        furtherNotes):
        """ Class to extract calibration factors and transmissions from the photodiode signals.
            The first two functions to be called must be:
                1. extract_calibration_factors
                2. extract_aperture_transmission (in this order)
            Only then is it sensible to call the function "extract_oa_ca_transmissions".

            Input:
            tot_num_of_pos  Integer denoting the total number of positions at which the photodiode
                            signals will be measured.
            sample_material String
            solvent         String
            concentration   Float
            laser_rep_rate  Integer
            furtherNotes    String
        """
        self.c_OA = None
        self.c_CA = None
        self.S = None                         # Transmission of aperture
        self.combined_c_CA = None             # S*c_CA

        self._tot_num_of_pos = tot_num_of_pos  # The total number of measurement stage positions.
        self.current_position_step = 0     # Integer indicating next empty transmission array entry.

        self.T_OA = np.zeros(shape=(self.tot_num_of_pos, 3))  # transmission in open aperture path.
        self.T_CA = np.zeros(shape=(self.tot_num_of_pos, 3))  # transmission in closed aperture path.

        self._sample_material = sample_material
        self._solvent = solvent
        self._concentration = concentration
        self._laser_rep_rate = laser_rep_rate
        self._furtherNotes = furtherNotes
        self._folder = None
        self._folder_number = None  # number suffix of the folder's name
        self._define_folder()

        self.pulse_energy = None  # in J

        self._w0 = CONSTANTS_beam_waist  # m waist of incident beam in vacuum
        self._λ = CONSTANTS_wavelength   # m in vacuum
        # mm Rayleigh length in vacuum
        self.zR = np.pi * self.w0**2 / self.λ * 1e3 * CONSTANTS_rayleigh_length_correction_factor

        self.fit_z0 = None    # fitted results
        self.fit_dΨ = None
        self.fit_dΦ = None

        self._n2 = None   # in units of cm^2/W.


    @property
    def tot_num_of_pos(self):
        return self._tot_num_of_pos

    @tot_num_of_pos.setter
    def tot_num_of_pos(self, value):
        self._tot_num_of_pos = value
        self.T_OA = np.zeros(shape=(self._tot_num_of_pos, 3))
        self.T_CA = np.zeros(shape=(self._tot_num_of_pos, 3))

    @property
    def sample_material(self):
        return self._sample_material

    @sample_material.setter
    def sample_material(self, value):
        self._sample_material = value
        self._define_folder()

    @property
    def solvent(self):
        return self._solvent

    @solvent.setter
    def solvent(self, value):
        self._solvent = value
        self._define_folder()

    @property
    def concentration(self):
        return self._concentration

    @concentration.setter
    def concentration(self, value):
        self._concentration = value

    @property
    def laser_rep_rate(self):
        return self._laser_rep_rate

    @laser_rep_rate.setter
    def laser_rep_rate(self, value):
        self._laser_rep_rate = value

    @property
    def furtherNotes(self):
        return self._furtherNotes

    @furtherNotes.setter
    def furtherNotes(self, value):
        self._furtherNotes = value

    @property
    def w0(self):
        return self._w0

    @w0.setter
    def w0(self, value):
        self._w0 = value
        self.zR = np.pi * self.w0**2 / self.λ * 1e3 * CONSTANTS_rayleigh_length_correction_factor

    @property
    def λ(self):
        return self._λ

    @λ.setter
    def λ(self, value):
        self._λ = value
        self.zR = np.pi * self.w0**2 / self.λ * 1e3 * CONSTANTS_rayleigh_length_correction_factor



    def extract_pulse_energy(self, ref_signal_max):
        """ This function extracts from the signal of the reference photodiode the pulse energy.
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

        self.pulse_energy = pulse_energy(np.array([mean, std]))

        return self.pulse_energy



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
        self.c_OA = average_ratio(oa_signal, ref_signal)
        self.c_CA = average_ratio(ca_signal, ref_signal)
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
        assert self.c_CA is not None and self.c_OA is not None
        self.S = average_ratio(ca_signal, ref_signal, calib_factor=self.c_CA)
        self.combined_c_CA[0] = self.c_CA[0] * self.S[0]
        self.combined_c_CA[1] = np.sqrt((self.c_CA[0]*self.S[1])**2 + (self.c_CA[1]*self.S[0])**2)
        return self.S


    def extract_oa_ca_transmissions(self, position, ref_signal, oa_signal, ca_signal):
        """ At any one probe position, this function evaluates the transmission in both the open and
            closed aperture path. The result is written into the arrays self.T_CA and self.T_OA in
            a 1-dim numpy array [position, transmission, error on transmission].
        """
        assert self.c_OA is not None and self.c_CA is not None and self.S is not None and \
            self.combined_c_CA is not None
        
        assert self.current_position_step < self.tot_num_of_pos

        transmission_oa = average_ratio(oa_signal, ref_signal, self.c_OA)
        transmission_ca = average_ratio(ca_signal, ref_signal, self.combined_c_CA)
        self.T_OA[self.current_position_step] = np.array([position, *transmission_oa])
        self.T_CA[self.current_position_step] = np.array([position, *transmission_ca])

        self.current_position_step += 1


    def _define_folder(self):
        now = datetime.date.today()
        today = "{0:4d}_{1:02d}_{2:02d}".format(now.year, now.month, now.day)

        # Attention, we should take care about the strings we pass to path.join!
        
        if self.solvent != "--":
            name = self.sample_material + "_" + self.solvent
        else:
            name = self.sample_material
            
        folder0 = os.path.join('..', 'Measurements', today, name)

        for i in range(1, 100):
            folder = folder0 + "_{0:02}".format(i)

            if not os.path.exists(folder):
                self.folder = folder
                self._folder_num = i
                break


    def check_and_create_folder(self):
        if not os.path.exists(self.folder):
            os.makedirs(self.folder)


    def store_transmission_data(self):
        """ Stores the transmission data in self.T_CA and self.T_OA into a file.
        """
        assert self.S is not None and self._folder_num is not None

        now = datetime.datetime.today()
        time = "{0:4d}.{1:02d}.{2:02d}  {3:02d}:{4:02d}".format(
            now.year, now.month, now.day, now.hour, now.minute)
        
        header = "Date of measurement: " + time + "\n" + \
                 "Sample material: " + self.sample_material + "\n" + \
                 "Solvent: " + self.solvent + "\n" + \
                 "Concentration: " + str(self.concentration) + "mmol/l\n" + \
                 "Laser rep. rate: " + str(self.laser_rep_rate) + "Hz\n" + \
                 "Pulse energy = ({0:.3f} +- {1:.3f})µJ\n".format(self.pulse_energy[0]*1e6, self.pulse_energy[1]*1e6) + \
                 "Aperture transm. S = {0:.3f} +- {1:.3f}".format(self.S[0], self.S[1]) + "\n" + \
                 "Further notes: " + self.furtherNotes + "\n" + \
                 "\n" + \
                 "Position / mm    T_OA    deltaT_OA    T_CA    deltaT_CA"
        
        # assert that position entries in T_OA and T_CA are identical
        assert (self.T_OA[:,0] == self.T_CA[:,0]).all()
        transmission_array = np.concatenate((self.T_OA, self.T_CA[:,1:]), axis=1)

        try:
            self.check_and_create_folder()
            file = os.path.join(self.folder, "transmission_data_{0:02}.dat".format(self._folder_num))
            np.savetxt(file, transmission_array, header=header, fmt="%10.4f")
        except Exception as ex:
            print("Storage of transmission data failed!!!!")
            traceback.print_exc()
            print("\n")


    def fit_transmission_data(self):

        def fit_OA_transmission(T_OA, T_CA, zR):
            T_OA_fitfunc = lambda z, z0, dΨ: T_OA_func(z, z0, zR, dΨ)
            guess_OA = CONSTANTS_guess_OA

            fit_OA, cov_OA = curve_fit(
                T_OA_fitfunc, T_OA[:,0], T_OA[:,1], sigma=T_OA[:,2], p0=guess_OA)
            std_err_OA = np.sqrt(np.diag(cov_OA))

            fit_z0 = np.array([fit_OA[0], std_err_OA[0]])
            fit_dΨ = np.array([fit_OA[1], std_err_OA[1]])

            return fit_z0, fit_dΨ

            
        def fit_CA_transmission(T_OA, T_CA, z0, zR, dΨ):
            T_CA_fitfunc = lambda z, dΦ: T_CA_func(z, z0, zR, dΦ, dΨ)

            fit_CA, cov_CA = curve_fit(T_CA_fitfunc, T_CA[:,0], T_CA[:,1], sigma=T_CA[:,2])
            std_err_CA = np.sqrt(np.diag(cov_CA))

            fit_dΦ = np.array([fit_CA[0], std_err_CA[0]])
            
            return fit_dΦ

        def fit_CA_transmission_only(T_OA, T_CA, zR):
            T_CA_fitfunc = lambda z, z0, dΦ: T_CA_func(z, z0, zR, dΦ, 0)-T_OA_func(z, z0, zR, 0)+1
            guess_CA = CONSTANTS_guess_CA

            fit_CA, cov_CA = curve_fit(
                T_CA_fitfunc, T_CA[:,0], T_CA[:,1], sigma=T_CA[:,2], p0=guess_CA)
            std_err_CA = np.sqrt(np.diag(cov_CA))

            fit_z0 = np.array([fit_CA[0], std_err_CA[0]])
            fit_dΦ = np.array([fit_CA[1], std_err_CA[1]])

            return fit_z0, fit_dΦ


        assert self.T_OA is not None and self.T_CA is not None


        # firstly, fit T_OA and use self.zR (Rayleigh length in vacuum)
        try:
            fit_z0, fit_dΨ = fit_OA_transmission(self.T_OA, self.T_CA, self.zR)
            # The waist should be close to z=22mm
            assert np.abs(fit_z0[0]-22) < CONSTANTS_tolerance_z0
            assert np.abs(fit_z0[1]/fit_z0[0]) < 1
            assert np.abs(fit_dΨ[1]/fit_dΨ[0]) < 1

            self.fit_z0 = fit_z0
            self.fit_dΨ = fit_dΨ

        except Exception as ex:
            print("Fit of open aperture data failed. Will try to fit closed aperture data only.")
            traceback.print_exc()
            print("\n")
            # I assume there is no nonlinear absorption. Hence I will set dΨ manually to zero:
            self.fit_dΨ = np.array([0,0])

            try:
                fit_z0, fit_dΦ = fit_CA_transmission_only(self.T_OA, self.T_CA, self.zR)
                # The waist should be close to z=22mm
                assert np.abs(fit_z0[0]-22) < CONSTANTS_tolerance_z0
                assert np.abs(fit_z0[1]/fit_z0[0]) < 1
                assert np.abs(fit_dΦ[1]/fit_dΦ[0]) < 1

                self.fit_z0 = fit_z0
                self.fit_dΦ = fit_dΦ
                self.compute_n2(self.fit_dΦ)

                return None

            except Exception as ex:
                print("The fit of closed aperture data only failed, as well!")
                traceback.print_exc()
                print("\n")
                self.fit_z0 = None
                self.fit_dΨ = None
                self.fit_dΦ = None

                return None


        # now, fit T_CA using the data from the fit of T_OA
        try:
            fit_dΦ = fit_CA_transmission(self.T_OA, self.T_CA, fit_z0[0], self.zR, fit_dΨ[0])
            assert np.abs(fit_dΦ[1]/fit_dΦ[0]) < 1

            self.fit_dΦ = fit_dΦ
            self.compute_n2(self.fit_dΦ)

        except Exception as ex:
            print("Fit of closed aperture data failed.")
            traceback.print_exc()
            print("\n")
            self.fit_dΦ = None

        return None


    def compute_n2(self, dΦ):
        """
        Input:
        dΦ: 1-dim numpy array of length 2, first entry being the fitted dΦ, second entry its error.

        Output (OnlyS in place replacement):
        n2: 1-dim numpy array of length 2, first entry being n2, second its error, both 
        in units of 1e-16 cm^2/W.
        """

        # Constants
        R = CONSTANTS_cuvette_reflectivity  #Measured reflectivity of the front boundary glass cuvette/air
        pulse_length = CONSTANTS_pulse_length_FWHM  #seconds
        
        assert self.pulse_energy is not None

        L_eff = 1e-3 # m   TODO: Measure alpha and compute L_eff correctly

        P0 = np.array([self.pulse_energy[0], self.pulse_energy[1]]) / pulse_length*(1-R)
        
        I0 = np.empty(shape=(2))
        I0[0] = 2*P0[0] / (np.pi*self.w0**2)
        I0[1] = 2*P0[1] / (np.pi*self.w0**2)
        k = 2*np.pi / self.λ
        
        n2 = dΦ[0] / (k*I0[0]*L_eff)
        dn2 = np.sqrt( (dΦ[1] / (k*I0[0]*L_eff))**2 + (dΦ[0]*I0[1] / (k*I0[0]*L_eff)**2)**2)
        self._n2 = np.array([n2, dn2]) * 1e4 # in units of cm^2/W


    def store_fit_results(self):

        if self.fit_z0 is not None:
            line0 = "z0: ({0:.3f} +- {1:.3f})mm".format(self.fit_z0[0], self.fit_z0[1])
        else:
            line0 = "z0: Could not be fitted"
        if self.fit_dΨ is not None:
            line1 = "dPsi: ({0:.3f} +- {1:.3f})".format(self.fit_dΨ[0], self.fit_dΨ[1])
        else:
            line1 = "dPsi: Could not be fitted"
        if self.fit_dΦ is not None:
            line2 = "dPhi: ({0:.3f} +- {1:.3f})\n".format(self.fit_dΦ[0], self.fit_dΦ[1])
            
            try:
                n2_exp = self.get_power_of_ten(self._n2[0])

                line2 += "\nn2 = ({0:.2f} +- {1:.2f})e{2} cm^2/W".format(
                    self._n2[0]/10**n2_exp, self._n2[1]/10**n2_exp, n2_exp)
            except Exception as ex:
                print("Problem finding power of ten of n2.")
                traceback.print_exc()
                line2 += "\nFinding the power of ten of n2 raised an exception."

        else:
            line2 = "dPhi: Could not be fitted, hence n2 could not be computed."

        now = datetime.datetime.today()
        time = "{0:4d}.{1:02d}.{2:02d}  {3:02d}:{4:02d}".format(
            now.year, now.month, now.day, now.hour, now.minute)

        header = "# Fit results\n# ---------------------\n# " + time + "\n#\n#\n"

        try:
            self.check_and_create_folder()
            file = os.path.join(self.folder, "fit_results_{0:02}.dat".format(self._folder_num))

            fhandle = open(file, 'w')
            fhandle.write(header)
            fhandle.write(line0 + "\n" + line1 + "\n" + line2 + "\n")

        except Exception as ex:
            print("Storage of fit results file failed!!!!")
            traceback.print_exc()
            print("\n")
        finally:
            fhandle.close()


    def plot_transmission_data(self):

        assert self.folder is not None and self.S is not None


        T_OA = self.T_OA
        T_CA = self.T_CA

        plt.figure(figsize=(8.5, 5.5))
        plt.errorbar(T_OA[:,0], T_OA[:,1], yerr=T_OA[:,2], linestyle="", marker="x", color="red", alpha=0.8, label="OA")
        plt.errorbar(T_CA[:,0], T_CA[:,1], yerr=T_CA[:,2], linestyle="", marker="x", color="black", alpha=0.8, label="CA")

        # Plot the fit functions if fit parameters exist
        if self.fit_z0 is not None and self.fit_dΨ is not None and self.fit_dΦ is not None:
            pos_vals = np.linspace(T_OA[0,0]-.5, T_OA[-1,0]+.5, 200)
            T_OA_vals = T_OA_func(pos_vals, self.fit_z0[0], self.zR, self.fit_dΨ[0])
            T_CA_vals = T_CA_func(pos_vals, self.fit_z0[0], self.zR, self.fit_dΦ[0], self.fit_dΨ[0])
            
            plt.plot(pos_vals, T_OA_vals, color="red")
            plt.plot(pos_vals, T_CA_vals, color="black")

        properties = "Sample: " + self.sample_material + \
            ",     Solvent: " + self.solvent + \
            ",     Concentration = {0}mmol/l".format(self.concentration) + "\n" + \
            r"$E_{Pulse}$" + " = ({0:.3f} $\pm$ {1:.3f})µJ".format(self.pulse_energy[0]*1e6, self.pulse_energy[1]*1e6) + \
            r",     $f_{Laser}$" + " = {0}Hz".format(self.laser_rep_rate) + \
            ",     S = ({0:.2f} $\pm$ {1:.2f})%".format(self.S[0]*100, self.S[1]*100)

        if self._furtherNotes != "---": #default value
            properties += "\nFurther notes: " + self._furtherNotes

        if self._n2 is not None:
            try:
                n2_exp = self.get_power_of_ten(self._n2[0])
                n2string = "$n_2$ = ({0:.2f} $\pm$ {1:.2f})e{2} cm$^2$/W".format(
                                self._n2[0]/10**n2_exp, self._n2[1]/10**n2_exp, n2_exp)

                properties += "\n" + n2string

            except Exception as ex:
                print("Problem finding the power of ten of n2")
                traceback.print_exc()

        plt.title(properties, fontsize=9)
        plt.xlabel("z / mm")
        plt.ylabel("Transmission")
        plt.legend()
        plt.grid()

        try:
            self.check_and_create_folder()
            file = os.path.join(self.folder, "plot_{0:02}.pdf".format(self._folder_num))
            plt.savefig(file, dpi=600)
        except Exception as ex:
            print("Storage of plot failed!!!!")
            traceback.print_exc()
            print("\n")

        plt.show()


    def evaluate_measurement_and_reinitialise(self):
        self.store_transmission_data()
        self.fit_transmission_data()
        self.store_fit_results()
        self.plot_transmission_data()
        self.reinitialise()


    def reinitialise(self):
        self.T_OA = np.zeros(shape=(self.tot_num_of_pos, 3))
        self.T_CA = np.zeros(shape=(self.tot_num_of_pos, 3))

        self.current_position_step = 0
        self.fit_z0 = None
        self.fit_dΨ = None
        self.fit_dΦ = None
        self._n2 = None
        self._define_folder()


    def get_power_of_ten(self, value):
        """ Finds the power of ten of value and returns it as an integer. The power of ten must be
            between -8 and -21!
        """
        assert value is not None

        for i in range(-8,-22,-1):
            if value / 10**i >= 1:
                break

        assert i > -21
        return i





if __name__ == '__main__':
    a = np.array([.5,.5,.5,.50005,.5])
    b = np.array([1,1,1,1,1])

    calib = average_ratio(a,b)
    print(calib)

    print(average_ratio(a,b,calib))



