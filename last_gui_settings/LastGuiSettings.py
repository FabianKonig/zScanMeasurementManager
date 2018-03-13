import traceback
import pickle
import os



CONSTANTS_gui_settings_persistence_file = os.path.join('.', 
                                                       'last_gui_settings',
                                                       'last_gui_settings.pkl')


class LastGuiSettings:
    def __init__(self, sample, solvent, concentration, laser_rep_rate, geom_sample_length,
        refr_index_sample, refr_index_ambient, furtherNotes, numPositions, samplingRate,
        samplesPerChannel, iterations, i_i0, alpha, eff_sample_length, attenuation_pd_ref,
        attenuation_sample, oa_pmt_exponent, ca_pmt_exponent):

        self.sample = sample
        self.solvent = solvent
        self.concentration = concentration
        self.laser_rep_rate = laser_rep_rate
        self.geom_sample_length = geom_sample_length
        self.refr_index_sample = refr_index_sample
        self.refr_index_ambient = refr_index_ambient
        self.furtherNotes = furtherNotes
        self.numPositions = numPositions
        self.samplingRate = samplingRate
        self.samplesPerChannel = samplesPerChannel
        self.iterations = iterations
        self.i_i0 = i_i0
        self.alpha = alpha
        self.eff_sample_length = eff_sample_length
        self.attenuation_pd_ref = attenuation_pd_ref
        self.attenuation_sample = attenuation_sample
        self.oa_pmt_exponent = oa_pmt_exponent
        self.ca_pmt_exponent = ca_pmt_exponent


def get_last_settings():
    try:
        with open(CONSTANTS_gui_settings_persistence_file, 'rb') as input:
            last_settings = pickle.load(input)
            return last_settings
    except FileNotFoundError as fnfe:
            print("No persisted documentation pickle file found.")
            print(fnfe)
            return None


def persist_last_settings(last_settings):
    try:
        with open(CONSTANTS_gui_settings_persistence_file, 'wb') as output:
            pickle.dump(last_settings, output, pickle.HIGHEST_PROTOCOL)
    except Exception as ex:
        print("Could not persist last GUI settings.")
        traceback.print_exc()

