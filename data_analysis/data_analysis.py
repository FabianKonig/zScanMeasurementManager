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

    assert data_1.shape == data_2.shape and len(data_1.shape) == len(data_2.shape) == 1

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
    def __init__(self, oa_signal, ca_signal, ref_signal):
        self.c_OA = self.extract_calibration_factor(oa_signal, ref_signal)
        self.c_CA = self.extract_calibration_factor(ca_signal, ref_signal)
        self.S = 1  # Transmission of aperture



    def extract_calibration_factor(pd_signal, ref_signal):
        """ Bla bla
        """
        return average_ratio(pd_signal, ref_signal)



if __name__ == '__main__':
    a = np.array([.5,.5,.5,.50005,.5])
    b = np.array([1,1,1,1,1])

    calib = average_ratio(a,b)
    print(calib)

    print(average_ratio(a,b,calib))



