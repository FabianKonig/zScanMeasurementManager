import nidaqmx
import numpy as np


# Reference photodiode.
# Dev2 denotes the fast PCI Nidaq device.
pd_ref_channel = "Dev2/ai0"
pd_oa_channel = "Dev2/ai1"
pd_ca_channel = "Dev2/ai2"

channels = [pd_ref_channel, pd_oa_channel, pd_ca_channel]



class NidaqReader:

    def __init__(self, sampling_rate, num_samples_per_chan, iterations):
        self.sampling_rate = sampling_rate
        self.num_samples_per_chan = num_samples_per_chan
        self.iterations = iterations


    def read_nidaq(self):
        """ Retrieves self.num_samples_per_chan samples per channel. The rate at which these data
            are acquired from the NIDAQ is given in units of Hz and cannot exceed
            250kHz/3channels=83333Hz.
            Output: 2-dim numpy array, the first dimension denoting the channel, the second
                    dimension contains self.num_samples_per_chan measurements.
        """
        with nidaqmx.Task() as task:
            for channel in channels:
                task.ai_channels.add_ai_voltage_chan(
                    channel, name_to_assign_to_channel="",
                    terminal_config=nidaqmx.constants.TerminalConfiguration.RSE, min_val=-10.0,
                    max_val=1.0, units=nidaqmx.constants.VoltageUnits.VOLTS, custom_scale_name="")

            task.timing.cfg_samp_clk_timing(
                self.sampling_rate,
                sample_mode = nidaqmx.constants.AcquisitionType.FINITE,
                samps_per_chan = self.num_samples_per_chan)
                
            signals = np.array(task.read(number_of_samples_per_channel=self.num_samples_per_chan))

        return signals


    def read_nidaq_one_channel_after_the_other(self):
        """ Same as read_nidaq, but I figured out that it is better to query each channel separately
            as otherwise there will be electronic reflections between the channels inside the Nidaq
            introducing offset  voltages.
        """

        signals = np.empty(shape=(len(channels), self.num_samples_per_chan))
        channel_index = 0

        for channel in channels:
            with nidaqmx.Task() as task:
                task.ai_channels.add_ai_voltage_chan(
                    channel, name_to_assign_to_channel="",
                    terminal_config=nidaqmx.constants.TerminalConfiguration.RSE, min_val=-10.0,
                    max_val=1.0, units=nidaqmx.constants.VoltageUnits.VOLTS, custom_scale_name="")

                task.timing.cfg_samp_clk_timing(
                    self.sampling_rate,
                    sample_mode = nidaqmx.constants.AcquisitionType.FINITE,
                    samps_per_chan = self.num_samples_per_chan)
                
                signals[channel_index] = np.array(
                    task.read(number_of_samples_per_channel=self.num_samples_per_chan))
                channel_index += 1

        return signals


    def get_nidaq_measurement_max_values(self):
        """ I figured out that the best way of obtaining reproducable data from the Nidaq
            measurments is by acquiring approx 20000 data from with its maximum sample rate
            (sample rate of 8333 and number of samples of 20000) and taking the maximum measured
            value for each channel. This function thus queries the Nidaq signals, stores the maximum
            value of each channel and repeats this measurement self.iterations times.

            Output: 
            2-dim numpy array with the first dimension denoting the channel and the second
            dimension is self.iterations entries long, each storing the maximum value of that
            measurement. The order of the channels is 1: ref, 2: oa, 3: ca.
        """
        
        max_signal_values = np.empty(shape=(len(channels), self.iterations))

        for value_index in range(self.iterations):
            signals = self.read_nidaq_one_channel_after_the_other()
            for channel_index in range(len(channels)):
                max_signal_values[channel_index, value_index] = signals[channel_index].max()

        return max_signal_values


    def filter_nidaq_signal_peaks(self, signals):
        """ Currently not in use.
            The photodiode signals retrieved with the NIDAQ usually consist of many entries which
            have been measured at a moment of time where no optical laser pulse was incident. As
            such, lot of the measured values have to be discarded. We discard all measured entries
            where one or more of the measured signals drop below a certian threshold with respect to
            their channel's max value. Example: If, for any given measurement of reference diode,
            closed aperture and open aperture diode, the measured value of, e.g., the open aperture
            photodiode signal drops below threshold*maximum(entire measured set of open aperture
            signal), this measurement is discarded.

            Input:
            2-dim array, first dimension being the measured channels, second dimension being the
            actual measurement.

            Output:
            2-dim array as input array, with discarded entries removed
        """
        threshold = 0.8

        accepted_signals_mask = np.empty(shape=signals.shape[1], dtype=bool)

        for signal in signals:
            accepted_signals_mask *= signal > threshold*signal.max()

        signals = signals[:, np.where(accepted_signals_mask)[0]]
        return signals


    def get_filtered_nidaq_signal(self):
        """ Currently not in use.
            Retrieves self.num_samples_per_chan samples for the channels with a self.sampling_rate.
            It then filters the measurements and returns them. See the docstrings of the
            corresponding functions.
            Returns a 2-dim numpy array with the first dimension denoting the channel and the second
            dimension denoting the measurements.
        """

        # We want at least 10 useful measurements:
        while(True):
            signals = self.read_nidaq()
            signals = self.filter_nidaq_signal_peaks(signals)
            if signals.shape[1] >= 10:
                break

        return signals



if __name__ == '__main__':
    
    import matplotlib.pyplot as plt

    nr = NidaqReader(250000, 120000, 8)
    signals = nr.get_nidaq_measurement_max_values()
    print(signals[0].mean(), signals[0].std())
    print(signals[1].mean(), signals[1].std())
    print(signals[2].mean(), signals[2].std())
    print((signals[1] / signals[0]).mean(), (signals[1] / signals[0]).std())
    print((signals[2] / signals[0]).mean(), (signals[2] / signals[0]).std())

    plt.plot(signals[0], alpha=0.5, linestyle="", marker="x", label="Ref")
    plt.plot(signals[1], alpha=0.5, linestyle="", marker="x", label="OA")    
    plt.plot(signals[2], alpha=0.5, linestyle="", marker="x", label="CA")
    plt.legend()
    plt.show()