import datetime
import os


CONSTANTS_beam_waist = 19.0537e-6  # waist of incident beam in vacuum in m
CONSTANTS_wavelength = 532e-9      # wavelength of incident beam in vacuum in m


class Documentation:

    def __init__(self, sample, solvent, concentration, laser_rep_rate, furtherNotes,
        refr_index_sample, refr_index_ambient, alpha, geom_sample_length, eff_sample_length,
        pulse_energy, eff_pulse_energy, S, λ_vac, w0):
        
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
        self.S = S                                      # aperture transmission, 1-dim np.array

        self.λ_vac = λ_vac                              # in m, float, vacuum wavelength of beam
        self.w0 = w0                                    # in m, float, beam waist of incident beam
        self.zR = np.pi*self.w0**2/self.λ_vac           # in m, float, Rayleigh length in vacuum


    @staticmethod
    def empty(self):
        self.__init__(None, None, None, None, None, None, None, None, None, None, None, None, None,
            None, None)


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

            if not os.path.exists(directory):
                break

        return directory, i


    def get_new_directory_for_storage(self):

        directory, folder_num = self.define_directory()

        if not os.path.exists(directory):
            os.makedirs(directory)

        return directory, folder_num


    def get_data_file_header(self):
        now = datetime.datetime.today()
        time = "{0:02d}.{1:02d}.{2:4d}  {3:02d}:{4:02d}".format(
            now.day, now.month, now.year, now.hour, now.minute)

        header = "Date of measurement:      " + time + "\n" + \
                 "Sample:                   " + self.sample + "\n" + \
                 "Solvent:                  " + self.solvent + "\n" + \
                 "Concentration:            " + self.concentration + "\n" + \
                 "Laser rep. rate:          {0}Hz\n".format(self.laser_rep_rate) + \
                 "Further notes:            " + furtherNotes + "\n" + \
                 "Pulse energy:             ({0:.3f} +- {1:.3f})µJ\n".format(
                    self.pulse_energy[0]*1e6, self.pulse_energy[1]*1e6) + \
                 "Eff. pulse energy:        ({0:.3f} +- {1:.3f})µJ\n".format(
                    self.eff_pulse_energy[0]*1e6, self.eff_pulse_energy[1]*1e6) + \
                 "Aperture transm. S:       {0:.3f} +- {1:.3f}\n".format(self.S[0], self.S[1]) + \
                 "Linear refractive index:  {0:.3f}\n".format(self.refr_index_sample) + \
                 "Ambient refractive index: {0:.3f}\n".format(self.refr_index_ambient) + \
                 "alpha:                    {0:.3f} mm^-1\n".format(self.alpha*1e3) + \
                 "Geometric sample length:  {0:.3f}mm\n".format(self.geom_sample_length*1e3) + \
                 "Effective sample length:  {0:.3f}mm\n".format(self.eff_sample_length*1e3) + \
                 "Wavelength vacuum:        {0:.3f}nm\n".format(self.λ_vac*1e9) + \
                 "Beam waist:               {0:.3f}µm\n".format(self.w0*1e6) + \
                 "Rayleigh length vacuum:   {0:.3f}mm\n".format(self.zR*1e3)
                 "\n" + \
                 "\n" + \
                 "Position / mm    T_OA    deltaT_OA    T_CA    deltaT_CA"

        return header


    def get_plot_header(self):
        header = "Sample: " + self.sample + \
                 ",     Solvent: " + self.solvent + \
                 ",     Concentration = " + self.concentration + "\n" + \
                 "$E_{Pulse}^{eff}$" + " = ({0:.3f} $\pm$ {1:.3f})µJ".format(
                    self.eff_pulse_energy[0]*1e6, self.eff_pulse_energy[1]*1e6) + \
                 ",     $f_{Laser}$" + " = {0}Hz".format(self.laser_rep_rate) + \
                 ",     S = ({0:.2f} $\pm$ {1:.2f})%".format(self.S[0]*100, self.S[1]*100)

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
        assert self.S is not None
        assert self.λ_vac is not None
        assert self.w0 is not None
