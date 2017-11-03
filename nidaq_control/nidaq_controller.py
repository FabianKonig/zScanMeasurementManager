if __name__ == '__main__':
    import nidaqmx

else:
    from . import nidaqmx

import numpy as np


# Reference photodiode
pd_ref_channel = "Dev1/ai0"
pd_oa_channel = "Dev1/ai1"
pd_ca_channel = "Dev1/ai2"

channels = [pd_ref_channel, pd_oa_channel, pd_ca_channel]



def read_nidaq(self, sampling_rate, num_samples_per_chan):
    """ Retrieves num_samples_per_chan samples per channel. The rate at which these data are
        acquired from the NIDAQ is given in units of Hz and cannot exceed 24000Hz.
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
