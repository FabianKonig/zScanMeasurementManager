import nidaqmx
import numpy as np


# Settings
default_sampling_rate = 15993        # 48000 is maximum, we use 3 channels, thus 16000. We use
                                     # an uneven number to avoid possible "ringing effects".
default_num_samp_per_chan = 50000
default_num_of_max_values = 3        # must be greater than 1

# Reference photodiode
pd_ref_channel = "Dev1/ai0"
pd_oa_channel = "Dev1/ai1"
pd_ca_channel = "Dev1/ai2"

channels = [pd_ref_channel, pd_oa_channel, pd_ca_channel]



def read_nidaq(sampling_rate=default_sampling_rate, num_samples_per_chan=default_num_samp_per_chan):
    """ Retrieves num_samples_per_chan samples per channel. The rate at which these data are
        acquired from the NIDAQ is given in units of Hz and cannot exceed 48kHz/3channels=16000Hz.
        Output: 2-dim numpy array, the first dimension denoting the channel, the second dimension
                contains num_samples_per_chan measurements.
    """
    with nidaqmx.Task() as task:
        for channel in channels:
            task.ai_channels.add_ai_voltage_chan(channel)

        task.timing.cfg_samp_clk_timing(
            sampling_rate,
            sample_mode = nidaqmx.constants.AcquisitionType.FINITE,
            samps_per_chan = num_samples_per_chan)
            
        signals = np.array(task.read(number_of_samples_per_channel=num_samples_per_chan))

    return signals


def get_nidaq_measurement_max_values(sampling_rate=default_sampling_rate,
    num_samples_per_chan=default_num_samp_per_chan, number_of_max_values=default_num_of_max_values):
    """ I figured out that the best way of obtaining reproducable data from the Nidaq measurments is
        by acquiring data from it for approx 1 second or longer with maximum sample rate
        (approx corresponds to a sample rate of 16000 and number of samples of 25000) and taking the
        maximum measured value for each channel. This function thus queries the Nidaq signals,
        stores the maximum value of each channel and repeats this measurement number_of_max_values
        often.
        Output: 2-dim numpy array with the first dimension denoting the channel and the second
        dimension is number_of_max_values entries long, each storing the maximum value of that
        measurement. The order of the channels is 1: ref, 2: oa, 3: ca.
    """
    
    max_signal_values = np.empty(shape=(len(channels), number_of_max_values))

    for value_index in range(number_of_max_values):
        signals = read_nidaq(sampling_rate, num_samples_per_chan)
        for channel_index in range(len(channels)):
            max_signal_values[channel_index, value_index] = signals[channel_index].max()

    return max_signal_values


def filter_nidaq_signal_peaks(signals):
    """ The photodiode signals retrieved with the NIDAQ usually consist of many entries which have
        been measured at a moment of time where no optical laser pulse was incident. As such, lot of
        the measured values have to be discarded. We discard all measured entries where one or more
        of the measured signals drop below a certian threshold with respect to their overall maximum
        value. Example: If, for any given measurement of reference diode, closed aperture and open
        aperture diode, the measured value of, e.g., the open aperture photodiode signal drops below
        threshold*maximum(entire measured set of open aperture signal), this measurement is
        discarded.

        Input:
        2-dim array, first dimension being the measured channels, second dimension being the actual
        measurement.

        Output:
        2-dim array as input array, with discarded entries removed
    """
    threshold = 0.5

    accepted_signals_mask = np.empty(shape=signals.shape[1], dtype=bool)

    for signal in signals:
        accepted_signals_mask *= signal > threshold*signal.max()

    signals = signals[:, np.where(accepted_signals_mask)[0]]
    return signals


def get_filtered_nidaq_signal(sampling_rate=default_sampling_rate,
    num_samples_per_chan=default_num_samp_per_chan):
    """ Retrieves num_samples_per_chan samples for the channels with a sampling_rate.
        It then filters the measurements and returns them. See the docstrings of the corresponding
        functions.
        Returns a 2-dim numpy array with the first dimension denoting the channel and the second
        dimension denoting the measurements.
    """

    # We want at least 10 useful measurements:
    while(True):
        signals = read_nidaq(sampling_rate, num_samples_per_chan)
        signals = filter_nidaq_signal_peaks(signals)
        if signals.shape[1] >= 10:
            break

    return signals



if __name__ == '__main__':
    
    import matplotlib.pyplot as plt
    signals = read_nidaq()
    print(signals[0].max(), signals[0].sum(), signals[0].max() / signals[0].sum())
    print(signals[1].max(), signals[1].sum(), signals[1].max() / signals[1].sum())
    print(signals[2].max(), signals[2].sum(), signals[2].max() / signals[2].sum())

    plt.plot(signals[0], alpha=0.5, linestyle="", marker="x", label="Ref")
    plt.plot(signals[1], alpha=0.5, linestyle="", marker="x", label="OA")    
    plt.plot(signals[2], alpha=0.5, linestyle="", marker="x", label="CA")
    plt.legend()
    plt.show()