import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math
import os

if __name__ == '__main__':
    import documentation
else:
    from . import documentation



# Constants:
CONSTANTS_rayleighLength_correction_factor = 1.7

CONSTANTS_pulse_length_FWHM = 15e-12  # laser pulse length in seconds



class zScanDataAnalyser:
    """ Class for postprocessing recorded zScan data. """

    def __init__(self, T_OA, T_CA, documentation, storage_dir, folder_num):
        """ Initialise by passing in open and closed aperture data and a Documentation instance.
        Input:
        T_OA           Open aperture transmission data, 1-dum numpy array of length 3, first entry
                       denoting the sample position, second entry the corresponding transmission
                       and third the error on the transmission.
        T_CA           Closed aperture transmission data, 1-dum numpy array of length 3, first entry
                       denoting the sample position, second entry the corresponding transmission
                       and third the error on the transmission.
        documentation  A Documentation instance storing all relevant information on this measurement
        storage_dir    String denoting the directory in which to store fit results and plots.
        folder_num     Integer denoting the number suffix of the folder's name.
        """
        self.T_OA = T_OA
        self.T_CA = T_CA
        self.documentation = documentation
        self.storage_dir = storage_dir
        self.folder_num = folder_num


    @staticmethod
    def init_from_file(self, file):
        """ Initialise by passing in a file that has been stored by the zScanDataProcessor. The file
            will be parsed for all relevant data.
        """
        pass  # To be done





    def combined_fit(self):

        zR = CONSTANTS_rayleighLength * CONSTANTS_rayleighLength_correction_factor
        
        fitfunc = lambda z, z0, dΨ: T_OA_func(z, z0, zR, dΨ)





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
        k = 2*np.pi / λ_vac                     # in 1/m
        
        n2 = dΦ[0] / (k*I0[0]*eff_length)
        dn2 = np.sqrt( (dΦ[1]/(k*I0[0]*eff_length))**2 + (dΦ[0]*I0[1]/(k*I0[0]*eff_length)**2)**2)
        
        return np.array([n2, dn2])  # in m^2/W





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
