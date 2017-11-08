import numpy as np
import datetime
import os
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


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
    def __init__(self):
        """ Class to extract calibration factors and transmissions from the photodiode signals.
            The first two functions to be called must be:
                1. extract_calibration_factors
                2. extract_aperture_transmission (in this order)
            Only then is it sensible to call the function "extract_oa_ca_transmissions".
        """
        self.c_OA = None
        self.c_CA = None
        self.S = None  # Transmission of aperture
        self.combined_c_CA = None  # S*c_CA

        self.T_OA = np.zeros(shape=(tot_num_of_pos, 3))  # transmission in open aperture path.
        self.T_CA = np.zeros(shape=(tot_num_of_pos, 3))  # transmission in closed aperture path.

        assert tot_num_of_pos > 2
        self.tot_num_of_pos = tot_num_of_pos   # The total number of measurement stage positions.
        self.current_position_step = 0     # Integer indicating next empty transmission array entry.

        self.sample_material = "default"
        self.solvent = "default"
        self.concentration = "default"
        self.laserfreq = "default"
        self.folder = self.define_folder()   # The folder to store results in

        self.w0 = 20.0299e-6  #m waist of incident beam in vacuum
        self.λ = 532e-9  #m in vacuum
        self.zR = np.pi * self.w0**2 / self.λ * 1e3   #mm Rayleigh length in vacuum

        self.fit_z0 = None   # fitted results
        self.fit_dΨ = None
        self.fit_dΦ = None


    def extract_calibration_factors(self, ref_signal, oa_signal, ca_signal):
        """ If the probe is placed into the beam far displaced from the focal spot, i.e. nonlinear
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


    def store_transmission_data(self, note):
        """ Stores the transmission data in self.T_CA and self.T_OA into a file.
        """
        now = datetime.datetime.today()
        time = "{0:4d}.{1:02d}.{2:02d}  {3:02d}:{4:02d}".format(
            now.year, now.month, now.day, now.hour, now.minute)
        header = note + " \n" + time + "\n\nPosition / mm    T_OA    deltaT_OA    T_CA    deltaT_CA"
        
        # assert that position entries in T_OA and T_CA are identical
        assert (self.T_OA[:,0] == self.T_CA[:,0]).all()
        transmission_array = np.concatenate((self.T_OA, self.T_CA[:,1:]), axis=1)

        try:
            file = os.path.join(self.folder, "data.csv")
            np.savetxt(file, transmission_array, header=header, fmt="%10.4f")
        except Exception as ex:
            print("Storage of transmission data failed!!!!")
            print(ex)
        

    def define_folder(self):
        """ 
        """
        now = datetime.date.today()
        today = "{0:4d}_{1:02d}_{2:02d}".format(now.year, now.month, now.day)
        directory = os.path.join('..', 'Measurements', today, self.material)  # Attention, we should take care about the strings we pass to path.join!

        if not os.path.exists(directory):
            os.makedirs(directory)
        self.folder = directory


    def fit_transmission_data(self):

        # abbreviation :-)
        T_OA = self.T_OA
        T_CA = self.T_CA

        # firstly, fit T_OA and use self.zR (Rayleigh length in vacuum)
        T_OA_fitfunc = lambda z, z0, dΨ: T_OA_func(z, z0, self.zR, dΨ)
        fit_OA, cov_OA = curve_fit(T_OA_fitfunc, T_OA[:,0], T_OA[:,1], sigma=T_OA[:,2])
        std_err_OA = np.sqrt(np.diag(cov_OA))

        fit_z0 = np.array([fit_OA[0], std_err_OA[0]])
        fit_dΨ = np.array([fit_OA[1], std_err_OA[1]])
        del fit_OA, cov_OA, std_err_OA


        # Here, we have to check whether the first fit was successful! There might be samples where
        # no nonlinear absorption occurs. Obviously, the fitting procedure will definitively fail in
        # those samples.


        # now, fit T_CA using the data from the fit to T_OA
        T_CA_fitfunc = lambda z, dΦ : T_CA_func(z, fit_z0[0], self.zR, dΦ, fit_dΨ[0])
        fit_CA, cov_CA = curve_fit(T_CA_fitfunc, T_CA[:,0], T_CA[:,1], sigma=T_CA[:,2])
        std_err_CA = np.sqrt(np.diag(cov_CA))

        fit_dΦ = np.array([fit_CA[0], std_err_CA[0]])
        del fit_CA, cov_CA, std_err_CA

        # store results in instance variables:
        self.fit_z0 = fit_z0
        self.fit_dΨ = fit_dΨ
        self.fit_dΦ = fit_dΦ



    def plot_transmission_data(self):

        assert self.folder is not None
        file = os.path.join(self.folder, "plot.pdf")

        T_OA = self.T_OA
        T_CA = self.T_CA
        plt.errorbar(T_OA[:,0], T_OA[:,1], yerr=T_OA[:,2], linestyle="", marker="x", color="black", label="OA")
        plt.errorbar(T_CA[:,0], T_CA[:,1], yerr=T_CA[:,2], linestyle="", marker="x", color="red", label="CA")

        # Plot the fit functions if fit parameters exist
        if self.fit_z0 is not None and self.fit_dΨ is not None and self.fit_dΦ is not None:
            pos_vals = np.linspace(T_OA[0,0]-.5, T_OA[-1,0]+.5, 200)
            T_OA_vals = T_OA_func(pos_vals, self.fit_z0[0], self.zR, self.fit_dΨ[0])
            T_CA_vals = T_CA_func(pos_vals, self.fit_z0[0], self.zR, self.fit_dΦ[0], self.fit_dΨ[0])
            
            plt.plot(pos_vals, T_OA_vals, color="black")
            plt.plot(pos_vals, T_CA_vals, color="red")

        plt.grid()
        plt.legend()
        plt.savefig(file, dpi=600)
        plt.show()


    def evaluate_measurement_and_reinitialise(self, note):
        self.define_folder(note)
        self.store_transmission_data(note)
        self.fit_transmission_data()
        self.plot_transmission_data()
        self.reinitialise()


    def reinitialise(self):
        self.T_OA = np.zeros(shape=(self.tot_num_of_pos, 3))
        self.T_CA = np.zeros(shape=(self.tot_num_of_pos, 3))

        self.current_position_step = 0
        self.folder = None





if __name__ == '__main__':
    a = np.array([.5,.5,.5,.50005,.5])
    b = np.array([1,1,1,1,1])

    calib = average_ratio(a,b)
    print(calib)

    print(average_ratio(a,b,calib))



