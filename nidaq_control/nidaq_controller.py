import nidaqmx
import numpy as np


# Reference photodiode.
# Dev2 denotes the fast PCI Nidaq device.
pd_ref_channel = "Dev2/ai0"
pd_oa_channel = "Dev2/ai1"
pd_ca_channel = "Dev2/ai2"

channels = [pd_ref_channel, pd_oa_channel, pd_ca_channel]



class NidaqReader:

    def __init__(self, sampling_rate, num_samples_per_chan):
        self.sampling_rate = sampling_rate
        self.num_samples_per_chan = num_samples_per_chan


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
                    terminal_config=nidaqmx.constants.TerminalConfiguration.RSE, min_val=-0.5,
                    max_val=10.0, units=nidaqmx.constants.VoltageUnits.VOLTS, custom_scale_name="")

            task.timing.cfg_samp_clk_timing(
                self.sampling_rate,
                sample_mode = nidaqmx.constants.AcquisitionType.FINITE,
                samps_per_chan = self.num_samples_per_chan)
                
            signals = np.array(task.read(number_of_samples_per_channel=self.num_samples_per_chan))

        return signals


    def read_nidaq_one_channel_after_the_other(self):
        """ Same as read_nidaq, but I figured out that it is better to query each channel separately
            as otherwise there will be electronic reflections between the channels inside the Nidaq
            introducing offset voltages.
        """

        signals = np.empty(shape=(len(channels), self.num_samples_per_chan))
        channel_index = 0

        for channel in channels:
            with nidaqmx.Task() as task:
                task.ai_channels.add_ai_voltage_chan(
                    channel, name_to_assign_to_channel="",
                    terminal_config=nidaqmx.constants.TerminalConfiguration.RSE, min_val=-0.5,
                    max_val=10.0, units=nidaqmx.constants.VoltageUnits.VOLTS, custom_scale_name="")

                task.timing.cfg_samp_clk_timing(
                    self.sampling_rate,
                    sample_mode = nidaqmx.constants.AcquisitionType.FINITE,
                    samps_per_chan = self.num_samples_per_chan)
                
                signals[channel_index] = np.array(
                    task.read(number_of_samples_per_channel=self.num_samples_per_chan))
                channel_index += 1


        # correct for possible channel offsets:
        signals = self.correct_offset(signals)

        return signals


    def correct_offset(self, signals):
        """ At large repetition rates, the photdiode signals do not have sufficient time to fall
            down to zero until the next pulse is incident. This causes an offset for the signals,
            artifically increasing the peak values. I epmirically figured out, that it is sufficient
            to modify the signals such that their minimum value is zero. This does neglect the
            first 2000 points of a NIDAQ measurement where there is some weird electronics going on
            which causes the signals to be very low.
            Input:
            signals, a 2-dim numpy array, the first dimension denoting the channel, the second
                     denotes the measurment point
            
            Output:
            signals (as in input) with minimum value being zero (neglecting the first 2000 points).
        """

        for i in range(len(signals)):
            signals[i] = signals[i] - signals[i, 2000::].min()

        return signals


    def peak_finder(self, rtn_peak_positions=False, rtn_raw_nidaq_signal=False):
        """ This function queries the Nidaq once. It then takes the signal of each channel,
            finds each peak in that signal and finds for each peak the maximum measured value.

            Output:
            2-dim numpy array. The first dimension denotes the channel, the second dimension stores
            all peak values found for that channel.
            If rtn_peak_positions==True: Returns a second 2-dim numpy array, the first dimension
            denoting the channel, the second the peak position.
            If rtn_raw_nidaq_signal==True: Returns a third 2-dim numpy array, the first dimension
            denoting the channel, the second the entire measured Nidaq signal.
        """

        # Firstly, read out Nidaq.
        signals = self.read_nidaq_one_channel_after_the_other()
        assert len(signals) == len(channels)

        # Allocate an array for each channel, each being as long as we could possibly find peaks.
        # This length is the maximum possible repetition rate of the pulsed laser multiplied by the
        # time the Nidaq is queried. Additional 100 for safety :)
        max_possible_peaks = 1100 * self.num_samples_per_chan // self.sampling_rate + 100
        peaks = np.empty(shape=(len(signals), max_possible_peaks), dtype=float)
        found_peak_index = np.zeros(shape=len(signals), dtype=int) # actual number of found peaks

        if rtn_peak_positions:
            peak_positions = np.empty(shape=(len(signals), max_possible_peaks), dtype=int)

        for channel_index in range(len(signals)):
            channel_length = len(signals[channel_index])
            assert channel_length == self.num_samples_per_chan

            threshold = signals[channel_index].max() * 0.95

            point_index = 0
            number_of_found_peaks = 0

            # Fo each measured point in the channel signal, test whether the point's value exceeds
            # the threshold. If so, it is a candidate for a peak value.
            while True:
                value = signals[channel_index, point_index]

                if  value > threshold:
                    last_found_peak, position = value, point_index

                    # A possible peak value is found. Test the next points whether they might be
                    # even larger than the one that is already found. If so, update last_found_peak.
                    # Examine all points like this until they fall below the threshold, meaning that
                    # we have found the maximum value of this peak for sure. Then store it and 
                    # go on finding the next peak in the upper while loop.
                    while True:
                        point_index += 1
                        if point_index >= channel_length:
                            break
                            
                        value = signals[channel_index, point_index]
                        if value > last_found_peak:
                             last_found_peak, position = value, point_index

                        # check if we are far from the interesting region: threshold*0.5
                        elif value < threshold*0.5 or point_index >= channel_length:
                            peaks[channel_index, found_peak_index[channel_index]] = last_found_peak
                    
                            if rtn_peak_positions:
                                peak_positions[channel_index, found_peak_index[channel_index]] \
                                    = position

                            found_peak_index[channel_index] += 1
                            break

                point_index += 1
                if point_index >= channel_length:
                    break

        # Finally, store all found peaks in an array of the correct length, i.e. the minimum 
        # number of found peaks in one of the channels.
        min_num_of_found_peaks = found_peak_index.min()
        peaks_result = np.empty(shape=(len(signals), min_num_of_found_peaks))

        if rtn_peak_positions:
            peak_positions_result = np.empty(shape=(len(signals), min_num_of_found_peaks))

        for i in range(min_num_of_found_peaks):
            for channel_index in range(len(signals)):
                peaks_result[channel_index, i] = peaks[channel_index, i]

                if rtn_peak_positions:
                    peak_positions_result[channel_index, i] = peak_positions[channel_index, i]

        if rtn_peak_positions and rtn_raw_nidaq_signal:
            return peaks_result, peak_positions_result, signals
        elif rtn_peak_positions and not rtn_raw_nidaq_signal:
            return peaks_result, peak_positions_result
        elif not rtn_peak_positions and rtn_raw_nidaq_signal:
            return peaks_result, signals
        else:
            return peaks_result


    def get_nidaq_measurement_max_values(self, iterations):
        """ I figured out that the best way of obtaining reproducable data from the Nidaq
            measurements is by acquiring approx 70000 data from it with its maximum sample rate
            (sample rate of 250k and number of samples of 70000) and taking the maximum measured
            value for each channel. This function thus queries the Nidaq signals, stores the maximum
            value of each channel and repeats this measurement iterations times.

            Output: 
            2-dim numpy array with the first dimension denoting the channel and the second
            dimension is iterations entries long, each storing the maximum value of that
            measurement. The order of the channels is 1: ref, 2: oa, 3: ca.
        """
        
        max_signal_values = np.empty(shape=(len(channels), iterations))

        for iteration_index in range(iterations):
            signals = self.read_nidaq_one_channel_after_the_other()
            for channel_index in range(len(channels)):
                max_signal_values[channel_index, iteration_index] = signals[channel_index].max()

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

    nr = NidaqReader(250000, 70000)
    #peaks, peak_positions, signals = nr.peak_finder(rtn_peak_positions=True, 
    #                                                rtn_raw_nidaq_signal=True)

    signals = nr.read_nidaq_one_channel_after_the_other()

    #plt.plot(peak_positions[0], peaks[0], color="brown", linestyle="", marker="+", markersize=8)
    #plt.plot(peak_positions[1], peaks[1], color="brown", linestyle="", marker="+", markersize=8)
    #plt.plot(peak_positions[2], peaks[2], color="brown", linestyle="", marker="+", markersize=8)


    print(signals[0,2000::].min())
    print(signals[1,2000::].min())
    print(signals[2,2000::].min())

    #plt.plot(signals[0]-signals[0,2000::].min()+(-.000177219), alpha=0.5, linestyle="", marker="x", label="Ref")
    #plt.plot(signals[1], alpha=0.5, linestyle="", marker="x", label="OA")    
    #plt.plot(signals[2], alpha=0.5, linestyle="", marker="x", label="CA")
    #plt.legend()
    #plt.show()

