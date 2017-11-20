import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math
import os
import traceback
import datetime

if __name__ == '__main__':
    import documentation
else:
    from . import documentation



# Constants:

# Correction factor in order to obtain a better fit:
CONSTANTS_rayleighLength_correction_factor = 1.7

# Initial guess for fit parameters. The first entry denotes the beam waist position in mm
CONSTANTS_guess_OA = [22,1]  # second entry: dΨ
CONSTANTS_guess_CA = [22,1]  # second entry: dΦ

CONSTANTS_pulse_length_FWHM = 15e-12  # laser pulse length in seconds



class zScanDataAnalyser:
    """ Class for postprocessing recorded zScan data. """

    def __init__(self, T_OA, T_CA, documentation):
        """ Initialise by passing in open and closed aperture data and a Documentation instance.
        Input:
        T_OA           Open aperture transmission data, 1-dum numpy array of length 3, first entry
                       denoting the sample position, second entry the corresponding transmission
                       and third the error on the transmission.
        T_CA           Closed aperture transmission data, 1-dum numpy array of length 3, first entry
                       denoting the sample position, second entry the corresponding transmission
                       and third the error on the transmission.
        documentation  A Documentation instance storing all relevant information on this measurement
        """
        self.T_OA = T_OA
        self.T_CA = T_CA
        self.T_CA_normalised = None
        self.documentation = documentation


        # Fit results and corresponding computed nonlinear refractive indices, all 2-dim numpy
        # arrays, the second dimension is of length 2, where the first entry denotes the value,
        # second its error.
        self.combined_fit = None          # n2, z0, dΦ, dΨ
        self.independent_CA_fit = None    # n2, z0, dΦ
        self.normalised_CA_fit = None     # n2, z0, dΦ


    @staticmethod
    def init_from_file(self, file):
        """ Initialise by passing in a file that has been stored by the zScanDataProcessor. The file
            will be parsed for all relevant data.
        """
        pass  # Make parsing and then "return zScanDataAnalyser(with the parsed data)".


    @property
    def T_OA(self):
        return self._T_OA


    @T_OA.setter
    def T_OA(self, value):
        self._T_OA = value

        if self.T_CA is not None:
            assert self.T_CA.shape() == value.shape()
            self.T_CA_normalised = self.compute_T_CA_normalised_array()


    @property
    def T_CA(self):
        return self._T_CA


    @T_CA.setter
    def T_CA(self, value):
        self._T_CA = value

        if self.T_OA is not None:
            assert self.T_OA.shape() == value.shape()
            self.T_CA_normalised = self.compute_T_CA_normalised_array()


    def compute_T_CA_normalised_array(self):
        """ Computes from self.T_CA and self.T_OA the normalised closed aperture data (i.e. closed
            aperture data devided by open aperture data).
            Returns that array.
        """
        T_OA = self.T_OA
        T_CA = self.T_CA
        assert T_OA is not None
        assert T_CA is not None

        T_CA_normalised = np.empty(shape=T_OA.shape())
        T_CA_normalised[:,0] = T_OA[:,0]
        T_CA_normalised[:,1] = T_CA[:,0] / T_OA[:,1]
        T_CA_normalised[:,2] = np.sqrt( 
            (T_CA[:,2]/T_OA[:,1])**2 + (T_CA[:,1]*T_OA[:,2] / T_OA[:,1]**2)**2 )

        return T_CA_normalised


    def combined_fit(self):
        """ First, the open aperture data are fitted and z0 and dΨ are retrieved from this fit.
            Then, the closed aperture data are fitted using the results from the first fit. dΦ is
            retrieved from that fit.
        """

        zR_for_fitting = self.documentation.zR * CONSTANTS_rayleighLength_correction_factor

        # Fit open aperture data.    
        try:
            fitfunc = lambda z, z0, dΨ: T_OA_func(z, z0, zR_for_fitting, dΨ)
            fit_z0, fit_dΨ = perform_fit(fitfunc, 
                                         self.T_OA[:,0],
                                         self.T_OA[:,1],
                                         sigma=self.T_OA[:,2],
                                         guess=CONSTANTS_guess_OA)
        except Exception as ex:
            print("Combined fit: Fit of open aperture data failed.")
            traceback.print_exc()
            print("\n")
            return None


        # Now, fit closed aperture data using the results from the open aperture data fit.
        try:
            fitfunc = lambda z, dΦ: T_CA_func(z, fit_z0[0], zR_for_fitting, dΦ, fit_dΨ[0])
            fit_dΦ = perform_fit(fitfunc, 
                                 self.T_CA[:,0],
                                 self.T_CA[:,1],
                                 sigma=self.T_CA[:,2],
                                 guess=CONSTANTS_guess_CA)
        except Exception as ex:
            print("Combined fit: Fit of closed aperture data failed.")
            traceback.print_exc()
            print("\n")
            return None

        self.combined_fit = np.empty(shape=(4,2))
        self.combined_fit[0] = self.compute_n2(fit_dΦ)
        self.combined_fit[1] = fit_z0
        self.combined_fit[2] = fit_dΦ
        self.combined_fit[3] = fit_dΨ


    def independent_ca_fit(self):
        """ The closed aperture data are fitted independently from the open aperture data as in
            self.combined_fit(). This fit retrieves z0 and dΦ. As for the fitfunction, dΨ = 0 is
            assumed.
        """

        zR_for_fitting = self.documentation.zR * CONSTANTS_rayleighLength_correction_factor

        try:
            fitfunc = lambda z, z0, dΦ: T_CA_normalised_func(z, z0, zR_for_fitting, dΦ)
            fit_z0, fit_dΦ = perform_fit(fitfunc, 
                                         self.T_CA[:,0],
                                         self.T_CA[:,1],
                                         sigma=self.T_CA[:,2],
                                         guess=CONSTANTS_guess_CA)
        except Exception as ex:
            print("Independent closed aperture fit failed.")
            traceback.print_exc()
            print("\n")
            return None

        self.independent_CA_fit = np.empty(shape=(3,2))
        self.independent_CA_fit[0] = self.compute_n2(fit_dΦ)
        self.independent_CA_fit[1] = fit_z0
        self.independent_CA_fit[2] = fit_dΦ


    def ca_normalised_wrt_oa_fit(self):
        """ This function fits the closed aperture data which are normalised with respect to the
            open aperture data and retrieves from this fit zo and dΦ.
        """
        zR_for_fitting = self.documentation.zR * CONSTANTS_rayleighLength_correction_factor

        try:
            fitfunc = lambda z, z0, dΦ: T_CA_normalised_func(z, z0, zR_for_fitting, dΦ)
            fit_z0, fit_dΦ = perform_fit(fitfunc, 
                                         self.T_CA_normalised[:,0],
                                         self.T_CA_normalised[:,1],
                                         sigma=self.T_CA_normalised[:,2],
                                         guess=CONSTANTS_guess_CA)
        except Exception as ex:
            print("Normalised closed aperture data fit failed.")
            traceback.print_exc()
            print("\n")
            return None

        self.normalised_CA_fit = np.empty(shape=(3,2))
        self.normalised_CA_fit[0] = self.compute_n2(fit_dΦ)
        self.normalised_CA_fit[1] = fit_z0
        self.normalised_CA_fit[2] = fit_dΦ



    def compute_n2(self, dΦ):
        """
        Input:
        dΦ: 1-dim numpy array of length 2, first entry being the fitted dΦ, second entry its error.
 
        Output:
        n2: 1-dim numpy array of length 2, first entry being n2, second its error, both in units
            of m^2/W.
        """
        pulse_length = CONSTANTS_pulse_length_FWHM              # seconds        
        eff_pulse_energy = self.documentation.eff_pulse_energy  # in J
        eff_length = self.documentation.eff_sample_length       # in m
        w0 = self.documentation.w0                              # m
        λ_vac = self.documentation.λ_vac                        # m


        P0 = np.array([eff_pulse_energy[0], eff_pulse_energy[1]]) / pulse_length  # in W
        
        I0 = np.empty(shape=(2))
        I0[0] = 2*P0[0] / (np.pi*w0**2)     # in W/m^2
        I0[1] = 2*P0[1] / (np.pi*w0**2)     # in W/m^2
        k = 2*np.pi / λ_vac                 # in 1/m
        
        n2 = dΦ[0] / (k*I0[0]*eff_length)
        dn2 = np.sqrt( (dΦ[1]/(k*I0[0]*eff_length))**2 + (dΦ[0]*I0[1]/(k*I0[0]*eff_length)**2)**2)
        
        return np.array([n2, dn2])  # in m^2/W


    def store_fit_results(self, directory, folder_num):
        """ Stores all fit results into a file in the provided directory. """

        now = datetime.datetime.today()
        time = "{0:02d}.{1:02d}.{3:4d}  {3:02d}:{4:02d}".format(
            now.day, now.month, now.year, now.hour, now.minute)
        header = "# Fit results\n# ---------------------\n# " + time + "\n#\n#\n"

        fit_results = np.array([self.combined_fit, self.independent_CA_fit, self.normalised_CA_fit])
        caption = ["combined fit:\n",
            "independent closed aperture fit:\n",
            "closed aperture normalised to open aperture fit:\n"]
        text = ""
        i=0

        for entry in fit_results:
            text += "Results of the " + caption[i]
            i += 1

            if entry is not None:
                n2 = entry[0] * 1e4  # in cm^2/W
                n2_exp = self.get_power_of_ten(n2[0])
                z0 = entry[1]
                dPhi = entry[2]
                if len(entry) > 3:
                    dPsi = entry[3]

                text += "n2: ({0:.3f} +- {1:.3f})e{2} cm^2/W\n".format(
                            n2[0]/10**n2_exp, n2[1]/10**n2_exp, n2_exp) + \
                        "z0: ({0:.3f} +- {1:.3f})mm\n".format(z0[0], z0[1]) + \
                        "dPhi: ({0:.3f} +- {1:.3f})".format(dPhi[0], dPhi[1])

                if len(entry) > 3:
                    text += "\ndPsi:({0:.3f} +- {1:.3f})".format(dPsi[0], dPsi[1])

            else:
                text += "No results."

            text += "#\n#\n"

        text = header + text


        try:
            file = os.path.join(directory, "fit_results_{0:02}.dat".format(folder_num))
            fhandle = open(file, 'w')
            fhandle.write(text)

        except Exception as ex:
            print("Storage of fit results failed!")
            traceback.print_exc()
            print("\n")
        finally:
            fhandle.close()



    def reinitialise(self):
        self.T_OA = None
        self.T_CA = None
        self.T_CA_normalised = None
        self.documentation = None
        self.combined_fit_z0 = None
        self.combined_fit_dΨ = None
        self.combined_fit_dΦ = None
        self.combined_fit_n2 = None
        self.independent_CA_fit_z0 = None
        self.independent_CA_fit_dΦ = None
        self.independent_CA_fit_n2 = None
        self.normalised_CA_fit_z0 = None
        self.normalised_CA_fit_dΦ = None
        self.normalised_CA_fit_n2 = None



# Fit functions:
def T_OA_func(z, z0, zR, dΨ):
    x = (z-z0)/zR
    return 1 - 2*(x**2+3)*dΨ / ((x**2+9) * (x**2+1))
 
def T_CA_func(z, z0, zR, dΦ, dΨ):
    x = (z-z0)/zR
    return T_OA_func(z, z0, zR, dΨ) + 4*x*dΦ / ((x**2+9) * (x**2+1))

def T_CA_normalised_func(z, z0, zR, dΦ):
    x = (z-z0)/zR
    return 1 + 4*x*dΦ / ((x**2+9) * (x**2+1))



def perform_fit(fitfunc, xdata, ydata, sigma=None, guess=None):

    assert len(xdata) == len(ydata)
    if sigma is not None:
        assert len(sigma) == len(xdata)
    if guess is not None:
        assert len(guess) == len(xdata)

    fit_params, cov = curve_fit(fitfunc, xdata, ydata, sigma=sigma, p0=guess)
    std = np.sqrt(np.diag(cov))

    number_of_parameters = len(fit_params)
    result = np.empty(shape=(number_of_parameters, 2))
    for i in range(number_of_parameters):
        assert np.abs(std[i]/fit_params[i]) < 1
        result[i] = np.array([fit_params[i], std[i]])

    return result


def get_power_of_ten(value):
    """ Finds the power of ten of value and returns it as an integer. The power of ten must be
        between 20 and -21!
    """
    assert value < 1e21
    assert value > 1e-21

    for i in range(20,-22,-1):
        if np.abs(value) / 10**i >= 1:
            return i



if __name__ == '__main__':
    pass
