import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math
import os
import traceback
import datetime
import re

if __name__ == '__main__':
    import sys
    sys.dont_write_bytecode = True
    import documentation_manually as documentation
else:
    from . import documentation



# Constants:

# Correction factor in order to obtain a better fit:
CONSTANTS_rayleighLength_correction_factor = 1.0

# Initial guess for fit parameters. The first entry denotes the beam waist position in mm
CONSTANTS_guess_OA = [23,1]  # second entry: q0
CONSTANTS_guess_CA = [23,1]  # second entry: dΦ



class zScanDataAnalyser:
    """ Class for postprocessing recorded zScan data. """

    def __init__(self, T_OA, T_CA, doc):
        """ Initialise by passing in open and closed aperture data and a Documentation instance.
        Input:
        T_OA           Open aperture transmission data, 1-dum numpy array of length 3, first entry
                       denoting the sample position, second entry the corresponding transmission
                       and third the error on the transmission.
        T_CA           Closed aperture transmission data, 1-dum numpy array of length 3, first entry
                       denoting the sample position, second entry the corresponding transmission
                       and third the error on the transmission.
        doc            A Documentation instance storing all relevant information on this measurement
        """
        self._T_OA = T_OA
        self._T_CA = T_CA
        self.T_CA_normalised = self.compute_T_CA_normalised_array()
        self.T_CA_reconstructed = self.compute_T_CA_reconstructed_array()
        self.doc = doc


        # Fit results and corresponding computed nonlinear refractive indices, all 2-dim numpy
        # arrays, the second dimension is of length 2, where the first entry denotes the value,
        # second its error.
        self.combined_fit = None          # n2, z0, dΦ, alpha2, q0
        self.independent_CA_fit = None    # n2, z0, dΦ
        self.normalised_CA_fit = None     # n2, z0, dΦ
        self.reconstructed_CA_fit = None  # n2, z0, dΦ


    @staticmethod
    def init_from_file(file):
        """ Initialise by passing in a file that has been stored by the zScanDataProcessor. The file
            will be parsed for all relevant data.
        """
        data = np.genfromtxt(file)
        T_OA = np.empty(shape=(len(data), 3))
        T_CA = np.empty(shape=(len(data), 3))
        T_OA[:,0] = data[:,0]
        T_OA[:,1] = data[:,1]
        T_OA[:,2] = data[:,2]
        T_CA[:,0] = data[:,0]    
        T_CA[:,1] = data[:,3]
        T_CA[:,2] = data[:,4]

        doc = documentation.Documentation.init_from_file(file)
        
        return zScanDataAnalyser(T_OA, T_CA, doc)


    @property
    def T_OA(self):
        return self._T_OA


    @T_OA.setter
    def T_OA(self, value):
        self._T_OA = value

        if self.T_CA is not None and value is not None:
            assert self.T_CA.shape == value.shape
            self.T_CA_normalised = self.compute_T_CA_normalised_array()
            self.T_CA_reconstructed = self.compute_T_CA_reconstructed_array()


    @property
    def T_CA(self):
        return self._T_CA


    @T_CA.setter
    def T_CA(self, value):
        self._T_CA = value

        if self.T_OA is not None and value is not None:
            assert self.T_OA.shape == value.shape
            self.T_CA_normalised = self.compute_T_CA_normalised_array()
            self.T_CA_reconstructed = self.compute_T_CA_reconstructed_array()


    def compute_T_CA_normalised_array(self):
        """ Computes from self.T_CA and self.T_OA the normalised closed aperture data (i.e. closed
            aperture data devided by open aperture data).
            Returns that array.
        """
        T_OA = self.T_OA
        T_CA = self.T_CA
        assert T_OA is not None
        assert T_CA is not None

        T_CA_normalised = np.empty(shape=T_OA.shape)
        T_CA_normalised[:,0] = T_OA[:,0]
        T_CA_normalised[:,1] = T_CA[:,1] / T_OA[:,1]
        T_CA_normalised[:,2] = np.sqrt( 
            (T_CA[:,2]/T_OA[:,1])**2 + (T_CA[:,1]*T_OA[:,2] / T_OA[:,1]**2)**2 )

        return T_CA_normalised


    def compute_T_CA_reconstructed_array(self):
        """ Reconstructs from self.T_CA and self.T_OA the purely refractive data (i.e. closed
            aperture data minus open aperture data).
            Returns that array.
        """
        T_OA = self.T_OA
        T_CA = self.T_CA
        assert T_OA is not None
        assert T_CA is not None

        T_CA_reconstructed = np.empty(shape=T_OA.shape)
        T_CA_reconstructed[:,0] = T_OA[:,0]
        T_CA_reconstructed[:,1] = T_CA[:,1] - T_OA[:,1] + 1
        T_CA_reconstructed[:,2] = np.sqrt( T_CA[:,2]**2 + T_OA[:,2]**2 )

        return T_CA_reconstructed


    def perform_combined_fit(self):
        """ First, the open aperture data are fitted and z0 and q0 are retrieved from this fit.
            Then, the closed aperture data are fitted using the results from the first fit. dΦ is
            retrieved from that fit.
        """

        zR_for_fitting = self.doc.zR * CONSTANTS_rayleighLength_correction_factor * 1e3  #mm

        # Fit open aperture data.    
        try:
            fitfunc = lambda z, z0, q0: T_OA_func(z, z0, zR_for_fitting, q0)
            fit_z0, fit_q0 = perform_fit(fitfunc, 
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
            fitfunc = lambda z, dΦ: T_CA_func(z, fit_z0[0], zR_for_fitting, dΦ, fit_q0[0])
            fit_dΦ, = perform_fit(fitfunc, 
                                 self.T_CA[:,0],
                                 self.T_CA[:,1],
                                 sigma=self.T_CA[:,2],
                                 guess=CONSTANTS_guess_CA[1]) # only provide guess for dΦ
        except Exception as ex:
            print("Combined fit: Fit of closed aperture data failed.")
            traceback.print_exc()
            print("\n")
            return None

        self.combined_fit = np.empty(shape=(5,2))
        self.combined_fit[0] = self.compute_n2(fit_dΦ)
        self.combined_fit[1] = fit_z0
        self.combined_fit[2] = fit_dΦ
        self.combined_fit[3] = self.compute_alpha2(fit_q0)
        self.combined_fit[4] = fit_q0


    def perform_ca_fit(self, xdata, ydata, sigma, guess, zR):
        """ Only helper function for the functions self.independent_ca_fit() and
            self.ca_normalised_wrt_oa_fit(). See the docstrings of those functions.
        """
        try:
            fitfunc = lambda z, z0, dΦ: T_CA_normalised_func(z, z0, zR, dΦ)
            fit_z0, fit_dΦ = perform_fit(fitfunc, 
                                         xdata,
                                         ydata,
                                         sigma=sigma,
                                         guess=guess)

        except Exception as ex:
            print("Closed aperture data fit failed.")
            traceback.print_exc()
            print("\n")
            return None

        result_array = np.empty(shape=(3,2))
        result_array[0] = self.compute_n2(fit_dΦ)
        result_array[1] = fit_z0
        result_array[2] = fit_dΦ

        return result_array


    def perform_independent_ca_fit(self):
        """ The closed aperture data are fitted independently from the open aperture data as in
            self.combined_fit(). This fit retrieves z0 and dΦ. As for the fitfunction, q0 = 0 is
            assumed.
        """
        zR_for_fitting = self.doc.zR * CONSTANTS_rayleighLength_correction_factor * 1e3  #mm

        self.independent_CA_fit = self.perform_ca_fit(self.T_CA[:,0],
                                                      self.T_CA[:,1],
                                                      self.T_CA[:,2],
                                                      CONSTANTS_guess_CA,
                                                      zR_for_fitting)


    def perform_ca_normalised_wrt_oa_fit(self):
        """ This function fits the closed aperture data which are normalised with respect to the
            open aperture data and retrieves from this fit z0 and dΦ.
        """
        zR_for_fitting = self.doc.zR * CONSTANTS_rayleighLength_correction_factor * 1e3  #mm

        self.normalised_CA_fit = self.perform_ca_fit(self.T_CA_normalised[:,0],
                                                     self.T_CA_normalised[:,1],
                                                     self.T_CA_normalised[:,2],
                                                     CONSTANTS_guess_CA,
                                                     zR_for_fitting)

    def perform_ca_reconstructed_fit(self):
        """ This function fits the closed aperture data which are reconstructed by subtracting the
            open aperture data from it and retrieves from this fit z0 and dΦ.
        """
        zR_for_fitting = self.doc.zR * CONSTANTS_rayleighLength_correction_factor * 1e3  #mm

        self.reconstructed_CA_fit = self.perform_ca_fit(self.T_CA_reconstructed[:,0],
                                                        self.T_CA_reconstructed[:,1],
                                                        self.T_CA_reconstructed[:,2],
                                                        CONSTANTS_guess_CA,
                                                        zR_for_fitting)


    def compute_n2(self, dΦ):
        """
        Input:
        dΦ: 1-dim numpy array of length 2, first entry being the fitted dΦ, second entry its error.
 
        Output:
        n2: 1-dim numpy array of length 2, first entry being n2, second its error, both in units
            of m^2/W.
        """
        eff_length = self.doc.eff_sample_length       # in m
        λ_vac = self.doc.λ_vac                        # m
        I0 = self.doc.eff_peak_intensity              # in W/m^2
        
        k = 2*np.pi / λ_vac                           # in 1/m
        
        n2 = dΦ[0] / (k*I0[0]*eff_length)
        dn2 = np.sqrt( (dΦ[1]/(k*I0[0]*eff_length))**2 + (dΦ[0]*I0[1]/(k*I0[0]**2*eff_length))**2 )
        
        return np.array([n2, dn2])  # in m^2/W


    def compute_alpha2(self, q0):
        """
        Input:
        q0: 1-dim numpy array of length 2, first entry being the fitted q0, second entry its error.
 
        Output:
        alpha2: 1-dim numpy array of length 2, first entry being alpha2, second its error, both in
                units of m/W.
        """
        eff_length = self.doc.eff_sample_length     # in m
        I0 = self.doc.eff_peak_intensity            # in W/m^2

        alpha2 = q0[0] / (I0[0]*eff_length)
        dalpha2 = np.sqrt( (q0[1]/(I0[0]*eff_length))**2 + (q0[0]*I0[1]/(I0[0]**2*eff_length))**2 )

        return np.array([alpha2, dalpha2])  # in m/W


    def store_fit_results(self, directory, folder_num):
        """ Stores all fit results into a file in the provided directory. """

        now = datetime.datetime.today()
        time = "{0:02d}.{1:02d}.{2:4d}  {3:02d}:{4:02d}".format(
            now.day, now.month, now.year, now.hour, now.minute)
        header = "Fit results\n---------------------\n" + time + "\n\n\n"

        caption = ["combined fit:\n",
            "independent closed aperture fit:\n",
            "closed aperture normalised to open aperture fit:\n",
            "reconstructed pure refraction data fit (CA-OA+1):\n"]
        text = ""
        i=0

        for fit_result in [self.combined_fit, self.independent_CA_fit, self.normalised_CA_fit,
                           self.reconstructed_CA_fit]:
            text += "Results of the " + caption[i]
            i += 1

            if fit_result is not None:
                n2 = fit_result[0] * 1e4  # in cm^2/W
                n2_exp = get_power_of_ten(n2[0])
                z0 = fit_result[1]
                dPhi = fit_result[2]
                if len(fit_result) > 3:  # combined fit
                    alpha2 = fit_result[3] * 1e2  # in cm/W
                    alpha2_exp = get_power_of_ten(alpha2[0])
                    q0 = fit_result[4]

                text += "n2: ({0:.4f} +- {1:.4f})e{2} cm^2/W\n".format(
                            n2[0]/10**n2_exp, n2[1]/10**n2_exp, n2_exp) + \
                        "z0: ({0:.4f} +- {1:.4f})mm\n".format(z0[0], z0[1]) + \
                        "dPhi: ({0:.4f} +- {1:.4f})".format(dPhi[0], dPhi[1])

                if len(fit_result) > 3:
                    text += "\nalpha2: ({0:.4f} +- {1:.4f})e{2} cm/W\n".format(
                                alpha2[0]/10**alpha2_exp, alpha2[1]/10**alpha2_exp, alpha2_exp) + \
                            "q0:({0:.4f} +- {1:.4f})".format(q0[0], q0[1])

            else:
                text += "No results."  # combined fit

            text += "\n\n"

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


    def plot_data(self, directory, folder_num):

        T_OA = self.T_OA
        T_CA = self.T_CA
        T_CA_normalised = self.T_CA_normalised
        T_CA_reconstructed = self.T_CA_reconstructed

        fig0, axes0 = plt.subplots(1, 1, figsize=(8.5, 5.5))
        fig1, axes1 = plt.subplots(1, 1, figsize=(8.5, 5.5))
        fig2, axes2 = plt.subplots(1, 1, figsize=(8.5, 5.5))
        fig3, axes3 = plt.subplots(1, 1, figsize=(8.5, 5.5))
        axes = [axes0, axes1, axes2, axes3]
        figures = [fig0, fig1, fig2, fig3]

        axes[0].errorbar(T_OA[:,0], T_OA[:,1], yerr=T_OA[:,2], linestyle="", marker="x", color="red", alpha=0.8, label="OA")
        axes[0].errorbar(T_CA[:,0], T_CA[:,1], yerr=T_CA[:,2], linestyle="", marker="x", color="black", alpha=0.8, label="CA")
        axes[1].errorbar(T_OA[:,0], T_OA[:,1], yerr=T_OA[:,2], linestyle="", marker="x", color="red", alpha=0.8, label="OA")
        axes[1].errorbar(T_CA[:,0], T_CA[:,1], yerr=T_CA[:,2], linestyle="", marker="x", color="black", alpha=0.8, label="CA")
        axes[2].errorbar(T_OA[:,0], T_OA[:,1], yerr=T_OA[:,2], linestyle="", marker="x", color="red", alpha=0.8, label="OA")
        axes[2].errorbar(T_CA_normalised[:,0], T_CA_normalised[:,1], yerr=T_CA_normalised[:,2], linestyle="", marker="x", color="black", alpha=0.8, label="CA/OA")
        axes[3].errorbar(T_OA[:,0], T_OA[:,1], yerr=T_OA[:,2], linestyle="", marker="x", color="red", alpha=0.8, label="OA")
        axes[3].errorbar(T_CA_reconstructed[:,0], T_CA_reconstructed[:,1], yerr=T_CA_reconstructed[:,2], linestyle="", marker="x", color="black", alpha=0.8, label="CA-OA+1")

        n2_header_string = ["", "", "", ""]
        alpha2_header_string = ["", "", "", ""]
        i = 0

        for fit_result in [self.combined_fit, self.independent_CA_fit, self.normalised_CA_fit,
                           self.reconstructed_CA_fit]:
            if fit_result is not None:
                n2 = fit_result[0] * 1e4  # in cm^2/W
                n2_exp = get_power_of_ten(n2[0])
                n2_header_string[i] = "\n$n_2$ = ({0:.2f} $\pm$ {1:.2f})e{2} cm$^2$/W".format(
                                        n2[0]/10**n2_exp, n2[1]/10**n2_exp, n2_exp)
                
                z0 = fit_result[1]
                dΦ = fit_result[2]
                if len(fit_result) == 5:  # combined fit plot
                    alpha2 = fit_result[3] * 1e2  # in cm/W
                    alpha2_exp = get_power_of_ten(alpha2[0])
                    alpha2_header_string[i] = r"   $\alpha_2$" + \
                                    " = ({0:.2f} $\pm$ {1:.2f})e{2} cm/W".format(
                                    alpha2[0]/10**alpha2_exp, alpha2[1]/10**alpha2_exp, alpha2_exp)
                    q0 = fit_result[4]

                pos_vals = np.linspace(T_OA[0,0]-.5, T_OA[-1,0]+.5, 500)
                zR_for_fitting = self.doc.zR * CONSTANTS_rayleighLength_correction_factor * 1e3 #mm

                if len(fit_result) == 5:  # combined fit plot
                    T_OA_vals = T_OA_func(pos_vals, z0[0], zR_for_fitting, q0[0])
                    T_CA_vals = T_CA_func(pos_vals, z0[0], zR_for_fitting, dΦ[0], q0[0])
                    axes[i].plot(pos_vals, T_OA_vals, color="red")
                    axes[i].plot(pos_vals, T_CA_vals, color="black")
                else:
                    T_CA_vals = T_CA_normalised_func(pos_vals, z0[0], zR_for_fitting, dΦ[0])
                    axes[i].plot(pos_vals, T_CA_vals, color="black")

            i += 1


        header = self.doc.get_plot_header()
        for i in range(len(axes)):
            axes[i].set_title(header + n2_header_string[i] + alpha2_header_string[i], fontsize=9)
            axes[i].set_xlabel("z / mm")
            axes[i].set_ylabel("Transmission")
            axes[i].legend()
            axes[i].grid()
            try:
                file = os.path.join(directory, "plot_{0:02}_{1}.pdf".format(folder_num, i))
                figures[i].savefig(file, dpi=600, bbox_inches='tight')
            except Exception as ex:
                print("Storage of plot {1} failed!!!!".format(i))
                traceback.print_exc()
                print("\n")

        plt.show()


    def reinitialise(self):
        self.T_OA = None
        self.T_CA = None
        self.T_CA_normalised = None
        self.doc = None
        self.combined_fit = None
        self.independent_CA_fit = None
        self.normalised_CA_fit = None



# Fit functions:
def T_OA_func(z, z0, zR, q0):
    x = (z-z0)/zR
    res = 0
    for m in range(0, 2): # only summation over m=0,1 seems to be sufficient for most applications
                          # (2 seems to be a very good compromise: If q<1, the series converges
                          #  and a large m would be better, but 2 yields already very good results.
                          #  The function actually is only approximated for q<<1. For q~1, the
                          #  function looks quite distorted and the series diverges for q>1. In
                          #  these cases m_max=1 seems to give quite reasonable results without the
                          #  series beginning to diverge!)
        res += (-q0)**m / ( (m+1)**1.5 * (1+x**2)**m )
    return res
 
def T_CA_func(z, z0, zR, dΦ, q0):
    x = (z-z0)/zR
    return T_OA_func(z, z0, zR, q0) + 4*x*dΦ / ((x**2+9) * (x**2+1))

def T_CA_normalised_func(z, z0, zR, dΦ):
    x = (z-z0)/zR
    return 1 + 4*x*dΦ / ((x**2+9) * (x**2+1))



def perform_fit(fitfunc, xdata, ydata, sigma=None, guess=None):

    assert len(xdata) == len(ydata)
    if sigma is not None:
        assert len(sigma) == len(xdata)

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
    try:
        assert np.abs(value) < 1e21
        assert np.abs(value) > 1e-21
    except Exception as ex:
        print("value is to large or to small for this function to find the power of ten!")
        traceback.print_exc()
        print("\n")
        return 0

    for i in range(20,-22,-1):
        if np.abs(value) / 10**i >= 1:
            return i



if __name__ == '__main__':

    # find the transmission data file in the current working directory:
    for line in os.listdir():
        res = re.search("transmission_data_(\d*).dat", line)
        if res:
            file, folder_num = res.group(), int(res.group(1))

    directory = os.getcwd()
    data_analyser = zScanDataAnalyser.init_from_file(file)
    data_analyser.perform_combined_fit()
    data_analyser.perform_independent_ca_fit()
    data_analyser.perform_ca_normalised_wrt_oa_fit()
    data_analyser.perform_ca_reconstructed_fit()
    data_analyser.store_fit_results(directory, folder_num)
    data_analyser.plot_data(directory, folder_num)
    data_analyser.reinitialise()
