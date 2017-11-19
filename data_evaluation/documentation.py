import datetime
import os

class Documentation:

    def __init__(self, sample, solvent, concentration, laser_rep_rate, furtherNotes,
        refr_index_sample, refr_index_ambient, alpha, geom_sample_length, eff_sample_length,
        pulse_energy, eff_pulse_energy, S):
        
        self.sample = sample
        self.solvent = solvent
        self.concentration = concentration              # in mol/l
        self.laser_rep_rate = laser_rep_rate            # in Hz
        self.furtherNotes = furtherNotes
        self.refr_index_sample = refr_index_sample
        self.refr_index_ambient = refr_index_ambient
        self.alpha = alpha                              # in 1/m
        self.geom_sample_length = geom_sample_length    # in m
        self.eff_sample_length = eff_sample_length      # in m

        self.pulse_energy = pulse_energy                # in J
        self.eff_pulse_energy = eff_pulse_energy        # in J
        self.S = S

        self.folder = None
        self.folder_num = None                          # number suffix of folder name
        self.define_folder()


    @property
    def sample(self):
        return self._sample
 

    @sample_material.setter
    def sample(self, value):
        self._sample = value
        self.define_folder()
 

    @property
    def solvent(self):
        return self._solvent
 

    @solvent.setter
    def solvent(self, value):
        self._solvent = value
        self.define_folder()


    def define_folder(self):
        now = datetime.date.today()
        today = "{0:4d}_{1:02d}_{2:02d}".format(now.year, now.month, now.day)

        # Attention, we should take care about the strings we pass to path.join!
        if self.solvent != "--":
            name = self.sample + "_" + self.solvent
        else:
            name = self.sample
            
        folder0 = os.path.join('..', 'Measurements', today, name)

        for i in range(1, 100):
            folder = folder0 + "_{0:02}".format(i)

            if not os.path.exists(folder):
                self.folder = folder
                self.folder_num = i
                break


    def get_folder_for_storage(self):
        if not os.path.exists(self.folder):
            os.makedirs(self.folder)

        return self.folder, self.folder_num


    def get_data_file_header(self):
        now = datetime.datetime.today()
        time = "{0:02d}.{1:02d}.{2:4d}  {3:02d}:{4:02d}".format(
            now.day, now.month, now.year, now.hour, now.minute)

        header = "Date of measurement:      " + time + "\n" + \
                 "Sample:                   " + self.sample + "\n" + \
                 "Solvent:                  " + self.solvent + "\n" + \
                 "Concentration:            {0:.2f}mmol/l\n".format(self.concentration*1e-3) + \
                 "Laser rep. rate:          {0}Hz\n".format(self.laser_rep_rate) + \
                 "Further notes:            " + furtherNotes + "\n" + \
                 "Pulse energy:             ({0:.3f} +- {1:.3f})µJ\n".format(self.pulse_energy[0]*1e6, self.pulse_energy[1]*1e6) + \
                 "Eff. pulse energy:        ({0:.3f} +- {1:.3f})µJ\n".format(self.eff_pulse_energy[0]*1e6, self.eff_pulse_energy[1]*1e6) + \
                 "Aperture transm. S:       {0:.3f} +- {1:.3f}\n".format(self.S[0], self.S[1]) + \
                 "Linear refractive index:  {0:.3f}\n".format(self.refr_index_sample) + \
                 "Ambient refractive index: {0:.3f}\n".format(self.refr_index_ambient) + \
                 "alpha:                    {0:.3f} mm^-1\n".format(self.alpha*1e3) + \
                 "Geometric sample length:  {0:.3f}mm\n".format(self.geom_sample_length*1e3) + \
                 "Effective sample length:  {0:.3f}mm\n".format(self.eff_sample_length*1e3) + \
                 "\n" + \
                 "Position / mm    T_OA    deltaT_OA    T_CA    deltaT_CA"

        return header