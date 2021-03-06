import numpy as np
import datetime
import os
import re
import glob



class Documentation:

    def __init__(self, sample, solvent, concentration, laser_rep_rate, furtherNotes,
        refr_index_sample, refr_index_ambient, alpha, geom_sample_length, eff_sample_length,
        pulse_energy, eff_pulse_energy, eff_peak_intensity, S, λ_vac, w0):
        
        self.sample = sample                            # string
        self.solvent = solvent                          # string
        self.concentration = concentration              # string
        self.laser_rep_rate = laser_rep_rate            # in Hz, float
        self.furtherNotes = furtherNotes                # string
        self.refr_index_sample = refr_index_sample      # float
        self.refr_index_ambient = refr_index_ambient    # float
        self.alpha = alpha                              # in 1/m, float
        self.geom_sample_length = geom_sample_length    # in m, float
        self.eff_sample_length = eff_sample_length      # in m, float

        self.pulse_energy = pulse_energy                # in J, 1-dim np.array
        self.eff_pulse_energy = eff_pulse_energy        # in J, 1-dim np.array
        self.eff_peak_intensity = eff_peak_intensity    # in W/m^2, 1-dim np.array
        self.S = S                                      # aperture transmission, 1-dim np.array

        self._λ_vac = λ_vac                             # in m, float, vacuum wavelength of beam
        self._w0 = w0                                   # in m, float, beam waist of incident beam
        self.zR = np.pi*self.w0**2/self.λ_vac           # in m, float, Rayleigh length in vacuum


    @staticmethod
    def empty(λ_vac, w0):
        return Documentation(None, None, None, None, None, None, None, None, None, None, None, None,
            None, None, λ_vac, w0)


    @staticmethod
    def init_from_file(file):

        attributes = {}

        with open(file, 'r', encoding='iso-8859-1') as f:
            for line in f:
                result = re.search("Sample: *(.*)", line)
                if result:
                    attributes["sample"] = result.group(1)
                    continue

                result = re.search("Solvent: *(.*)", line)
                if result:
                    attributes["solvent"] = result.group(1)
                    continue

                result = re.search("Concentration: *(.*)", line)
                if result:
                    attributes["concentration"] = result.group(1)
                    continue

                result = re.search("Laser rep. rate: *(\d*)Hz", line)
                if result:
                    attributes["laser_rep_rate"] = int(result.group(1))
                    continue

                result = re.search("Further notes: *(.*)", line)
                if result:
                    attributes["furtherNotes"] = result.group(1)
                    continue

                result = re.search("Aperture transm. S: *([\d.]*) \+\- ([\d.]*)", line)
                if result:
                    attributes["S"] = np.array([float(result.group(1)), float(result.group(2))])
                    continue

                result = re.search("Linear refractive index: *(.*)", line)
                if result:
                    attributes["refr_index_sample"] = float(result.group(1))
                    continue

                result = re.search("Ambient refractive index: *(.*)", line)
                if result:
                    attributes["refr_index_ambient"] = float(result.group(1))
                    continue

                result = re.search("Geometric sample length: *(.*)mm", line)
                if result:
                    attributes["geom_sample_length"] = float(result.group(1))*1e-3
                    continue

                result = re.search("Effective sample length: *(.*)mm", line)
                if result:
                    attributes["eff_sample_length"] = float(result.group(1))*1e-3
                    continue

                result = re.search("alpha: *(.*) mm\^\-1", line)
                if result:
                    attributes["alpha"] = float(result.group(1))*1e3
                    continue

                result = re.search("Pulse energy: *\(([\d.]*) \+\- ([\d.]*)\)µJ", line)
                if result:
                    attributes["pulse_energy"] = np.array([float(result.group(1)), 
                                                           float(result.group(2))]) * 1e-6
                    continue

                result = re.search("Eff. pulse energy: *\(([\d.]*) \+\- ([\d.]*)\)µJ", line)
                if result:
                    attributes["eff_pulse_energy"] = np.array([float(result.group(1)), 
                                                               float(result.group(2))]) * 1e-6
                    continue

                result = re.search("Eff. peak intensity: *\(([\d.]*) \+\- ([\d.]*)\)MW/cm\^2", line)
                if result:
                    attributes["eff_peak_intensity"] = np.array([float(result.group(1)), 
                                                                 float(result.group(2))]) * 1e10
                    continue

                result = re.search("Wavelength vacuum: *([\d.]*)nm", line)
                if result:
                    attributes["λ_vac"] = float(result.group(1))*1e-9
                    continue

                result = re.search("Beam waist: *([\d.]*)µm", line)
                if result:
                    attributes["w0"] = float(result.group(1))*1e-6
                    continue

        doc = Documentation.empty(attributes["λ_vac"], attributes["w0"])

        for key in attributes:
            setattr(doc, key, attributes[key])

        return doc


    @property
    def w0(self):
        return self._w0
 

    @w0.setter
    def w0(self, value):
        self._w0 = value
        self.zR = np.pi*self.w0**2/self.λ_vac


    @property
    def λ_vac(self):
        return self._λ_vac
 

    @λ_vac.setter
    def λ_vac(self, value):
        self._λ_vac = value
        self.zR = np.pi*self.w0**2/self.λ_vac


    def define_directory(self):
        now = datetime.date.today()
        today = "{0:4d}_{1:02d}_{2:02d}".format(now.year, now.month, now.day)

        # Attention, we should take care about the strings we pass to path.join!
        if self.solvent != "--":
            name = self.sample + "_" + self.solvent
        else:
            name = self.sample
            
        directory0 = os.path.join('..', 'Measurements', today, name)

        for i in range(1, 100):
            directory = directory0 + "_{0:02}".format(i)

            if not Documentation.exists_dir(directory):
                break

        return directory, i


    def get_new_directory_for_storage(self):
        directory, folder_num = self.define_directory()

        if not Documentation.exists_dir(directory):
            os.makedirs(directory)

        return directory, folder_num


    @staticmethod
    def exists_dir(directory):
        """ Checks whether the given directory exists. Returns True even if a directory exists
            with the given name and an additional suffix: full_directory = directory + "AnySuffix".
            Input:
            directory    String

            Output: 
            True or False, depending on whether the directory exists or not.
        """

        # Get the number of directories that match the pattern: directory + "*".
        pattern_match_number = len(glob.glob(directory + "*"))

        # In our case, each directory (which includes a folder number) should exist at maximum once,
        # as we want each folder number to occur only once in order to aid as an identifier.
        assert pattern_match_number == 0 or pattern_match_number == 1

        if pattern_match_number == 0:
            return False
        elif pattern_match_number == 1:
            return True


    def get_data_file_header(self):
        now = datetime.datetime.today()
        time = "{0:02d}.{1:02d}.{2:4d}  {3:02d}:{4:02d}".format(
            now.day, now.month, now.year, now.hour, now.minute)

        header = "Date of measurement:      " + time + "\n" + \
                 "Sample:                   " + self.sample + "\n" + \
                 "Solvent:                  " + self.solvent + "\n" + \
                 "Concentration:            " + self.concentration + "\n" + \
                 "Laser rep. rate:          {0}Hz\n".format(self.laser_rep_rate) + \
                 "Further notes:            " + self.furtherNotes + "\n" + \
                 "\n" + \
                 "Aperture transm. S:       {0:.4f} +- {1:.4f}\n".format(self.S[0], self.S[1]) + \
                 "Linear refractive index:  {0:.3f}\n".format(self.refr_index_sample) + \
                 "Ambient refractive index: {0:.3f}\n".format(self.refr_index_ambient) + \
                 "Geometric sample length:  {0:.3f}mm\n".format(self.geom_sample_length*1e3) + \
                 "Effective sample length:  {0:.3f}mm\n".format(self.eff_sample_length*1e3) + \
                 "alpha:                    {0:.8f} mm^-1\n".format(self.alpha*1e-3) + \
                 "\n" + \
                 "Pulse energy:             ({0:.6f} +- {1:.6f})µJ\n".format(
                    self.pulse_energy[0]*1e6, self.pulse_energy[1]*1e6) + \
                 "Eff. pulse energy:        ({0:.6f} +- {1:.6f})µJ\n".format(
                    self.eff_pulse_energy[0]*1e6, self.eff_pulse_energy[1]*1e6) + \
                 "Eff. peak intensity:      ({0:.3f} +- {1:.3f})MW/cm^2\n".format(
                    self.eff_peak_intensity[0]*1e-10, self.eff_peak_intensity[1]*1e-10) + \
                 "\n" + \
                 "Wavelength vacuum:        {0:.3f}nm\n".format(self.λ_vac*1e9) + \
                 "Beam waist:               {0:.3f}µm\n".format(self.w0*1e6) + \
                 "Rayleigh length vacuum:   {0:.3f}mm\n".format(self.zR*1e3) + \
                 "\n" + \
                 "\n" + \
                 "Position / mm    T_OA    deltaT_OA    T_CA    deltaT_CA"

        return header


    def get_plot_header(self):
        header = "Sample: " + self.sample + \
                 ",     Solvent: " + self.solvent + \
                 ",     Concentration = " + self.concentration + "\n" + \
                 "$I_{0}^{eff}$" + " = ({0:.3f} $\pm$ {1:.3f})MW/cm$^2$".format(
                    self.eff_peak_intensity[0]*1e-10, self.eff_peak_intensity[1]*1e-10) + \
                 ",     $f_{Laser}$" + " = {0}Hz".format(self.laser_rep_rate) + \
                 ",     S = ({0:.2f} $\pm$ {1:.2f})%".format(self.S[0]*100, self.S[1]*100)

        if self.furtherNotes != "---":  # default value
            header += "\nFurther notes: " + self.furtherNotes

        return header


    def assert_completeness(self):
        assert self.sample is not None
        assert self.solvent is not None
        assert self.concentration is not None
        assert self.laser_rep_rate is not None
        assert self.furtherNotes is not None
        assert self.refr_index_sample is not None
        assert self.refr_index_ambient is not None
        assert self.alpha is not None
        assert self.geom_sample_length is not None
        assert self.eff_sample_length is not None
        assert self.pulse_energy is not None
        assert self.eff_pulse_energy is not None
        assert self.eff_peak_intensity is not None
        assert self.S is not None
        assert self.λ_vac is not None
        assert self.w0 is not None
