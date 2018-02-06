from PyQt5 import QtWidgets
from math import isclose
import sys
import numpy as np
import last_gui_settings
import gui_design
import stage_control
import nidaq_control
import data_evaluation



# TODO:
# -----------------------------------
# - A preselect box for different materials would be helpful: Refractive index, alpha,
#   geom sample length, ambient refractive index...



# Constants:
CONSTANTS_beam_waist = 19.0537e-6     # waist of incident beam in vacuum in m
CONSTANTS_wavelength = 532e-9         # wavelength of incident beam in vacuum in m
CONSTANTS_pulse_length_FWHM = 15e-12  # laser pulse length in seconds




class Window(QtWidgets.QMainWindow, gui_design.Ui_MainWindow):
    def __init__(self):
        super().__init__()    # call __init__ of QtWidgets.QMainWindow
        self.setupUi(self)    # call setupUI of gui_design.Ui_MainWindow (generated with QtDesigner)
        self.defineSignalsSlots()

        self.doc = data_evaluation.Documentation.empty(CONSTANTS_wavelength, CONSTANTS_beam_waist)
        self.doc.alpha = 0  # initially
        self.doc.eff_sample_length = self.doubleSpinBox_geomSampleLength.value()*1e-3
        self.compile_documentation()

        self.data_processor = data_evaluation.zScanDataProcessor()

        self.nidaq_reader = nidaq_control.NidaqReader(
            self.spinBox_samplingRate.value(), 
            self.spinBox_samplesPerChannel.value())

        QtWidgets.QMessageBox.information(self, "Stage position initialisation",
            "The stages will now be moved home and subsequently to their initial positions." +
            " Make sure the stages can move unhindered!")
        
        self.stage_controller = stage_control.APT_Controller()

        # New function: Load input values of last GUI access:
        self.load_last_gui_settings()


    def defineSignalsSlots(self):
        self.pushButton_computeAlpha.clicked.connect(self.onClick_computeAlpha)
        self.pushButton_calibratePDs.clicked.connect(self.onClick_calibratePhotodiodes)
        self.pushButton_measureAperture.clicked.connect(self.onClick_measureApertureTransmission)
        self.pushButton_startMeasurement.clicked.connect(self.onClick_startMeasurement)

        self.lineEdit_sample.textChanged.connect(self.compile_documentation)
        self.lineEdit_solvent.textChanged.connect(self.compile_documentation)
        self.lineEdit_concentration.textChanged.connect(self.compile_documentation)
        self.spinBox_laserRepRate.valueChanged.connect(self.compile_documentation)
        self.doubleSpinBox_geomSampleLength.valueChanged.connect(self.compile_documentation)
        self.lineEdit_furtherNotes.textChanged.connect(self.compile_documentation)
        self.doubleSpinBox_refrIndexSample.valueChanged.connect(self.compile_documentation)
        self.doubleSpinBox_refrIndexAmbient.valueChanged.connect(self.compile_documentation)

        self.doubleSpinBox_geomSampleLength.valueChanged.connect(self.reinitAlphaAndEffLength)
        self.doubleSpinBox_refrIndexSample.valueChanged.connect(self.reinitAlphaAndEffLength)
        self.doubleSpinBox_refrIndexAmbient.valueChanged.connect(self.reinitAlphaAndEffLength)

        self.spinBox_samplingRate.valueChanged.connect(self.onNidaqParamsChange)
        self.spinBox_samplesPerChannel.valueChanged.connect(self.onNidaqParamsChange)
        self.spinBox_iterations.valueChanged.connect(self.onNidaqParamsChange)


    def compile_documentation(self):
        if self.lineEdit_solvent.text() == "":
            self.lineEdit_solvent.setText("--")  # default value

        if self.lineEdit_furtherNotes.text() == "":
            self.lineEdit_furtherNotes.setText("---")  # default value

        self.doc.sample = self.lineEdit_sample.text()
        self.doc.solvent = self.lineEdit_solvent.text()
        self.doc.concentration = self.lineEdit_concentration.text()
        self.doc.laser_rep_rate = self.spinBox_laserRepRate.value()
        self.doc.geom_sample_length = self.doubleSpinBox_geomSampleLength.value()*1e-3
        self.doc.furtherNotes = self.lineEdit_furtherNotes.text()
        self.doc.refr_index_sample = self.doubleSpinBox_refrIndexSample.value()
        self.doc.refr_index_ambient = self.doubleSpinBox_refrIndexAmbient.value()
        self.doc.λ_vac = CONSTANTS_wavelength
        self.doc.w0 = CONSTANTS_beam_waist


    def load_last_gui_settings(self):
        last_settings = last_gui_settings.get_last_settings()
        if last_settings is not None:
            self.update_gui_with_last_settings(last_settings)


    def update_gui_with_last_settings(self, last_settings):
        self.lineEdit_sample.setText(last_settings.sample)
        self.lineEdit_solvent.setText(last_settings.solvent)
        self.lineEdit_concentration.setText(last_settings.concentration)
        self.spinBox_laserRepRate.setValue(last_settings.laser_rep_rate)
        self.doubleSpinBox_geomSampleLength.setValue(last_settings.geom_sample_length)
        self.doubleSpinBox_refrIndexSample.setValue(last_settings.refr_index_sample)
        self.doubleSpinBox_refrIndexAmbient.setValue(last_settings.refr_index_ambient)
        self.lineEdit_furtherNotes.setText(last_settings.furtherNotes)
        self.spinBox_numPositions.setValue(last_settings.numPositions)
        self.spinBox_samplingRate.setValue(last_settings.samplingRate)
        self.spinBox_samplesPerChannel.setValue(last_settings.samplesPerChannel)
        self.spinBox_iterations.setValue(last_settings.iterations)
        self.doubleSpinBox_II0.setValue(last_settings.i_i0)
        self.label_alphaValue.setText(last_settings.alpha)
        self.label_effSampleLengthValue.setText(last_settings.eff_sample_length)
        self.doubleSpinBox_attenuationPdRef.setValue(last_settings.attenuation_pd_ref)

        self.doc.alpha = float(self.label_alphaValue.text())
        self.doc.eff_sample_length = float(self.label_effSampleLengthValue.text())*1e-3


    def onNidaqParamsChange(self):
        self.nidaq_reader.sampling_rate = self.spinBox_samplingRate.value()
        self.nidaq_reader.num_samples_per_chan = self.spinBox_samplesPerChannel.value()


    def reinitAlphaAndEffLength(self):
        """ If the geometrical sample length, or one of the refractive indices has been changed by
            the user, alpha and effective length have to be computed again with these new values.
            Before having done so, alpha is set to 0 again and the effective length to the geometric
            length.
        """
        self.label_effSampleLengthValue.setText("{0:.3f}".format(
            self.doubleSpinBox_geomSampleLength.value()))
        self.doc.eff_sample_length = self.doubleSpinBox_geomSampleLength.value() * 1e-3

        self.label_alphaValue.setText("{0:.3f}".format(0.))
        self.doc.alpha = 0


    def onClick_computeAlpha(self):
        transmission = self.doubleSpinBox_II0.value()
        refr_index_sample = self.doc.refr_index_sample
        refr_index_ambient = self.doc.refr_index_ambient
        geom_sample_length = self.doc.geom_sample_length  # in m

        alpha = self.data_processor.compute_alpha(transmission, refr_index_sample,
            refr_index_ambient, geom_sample_length)  # in 1/m
        self.doc.alpha = alpha
        self.label_alphaValue.setText("{0:.3f}".format(alpha))  # in 1/m
        
        eff_sample_length = self.data_processor.compute_effective_length(geom_sample_length,
            alpha)  # in m
        self.doc.eff_sample_length = eff_sample_length
        self.label_effSampleLengthValue.setText("{0:.3f}".format(eff_sample_length * 1e3))  # in mm


    def onClick_calibratePhotodiodes(self):
        signals = self.nidaq_reader.get_nidaq_measurement_max_values(
            self.spinBox_iterations.value())
        
        calib_factors = list(self.data_processor.extract_calibration_factors(*signals))
        self.label_cOAValue.setText("{0:.3f} +- {1:.3f}".format(*calib_factors[0]))
        self.label_cCAValue.setText("{0:.3f} +- {1:.3f}".format(*calib_factors[1]))
        
        pulse_energy = self.data_processor.extract_pulse_energy(
            signals[0] * self.doubleSpinBox_attenuationPdRef.value())
        self.doc.pulse_energy = pulse_energy
        self.label_pulseEnergyValue.setText("{0:.3f} +- {1:.3f}".format(*pulse_energy*1e6))  # in µJ

        eff_pulse_energy = self.data_processor.compute_effective_pulse_energy(
            pulse_energy,
            self.doc.refr_index_sample,
            self.doc.refr_index_ambient)
        self.doc.eff_pulse_energy = eff_pulse_energy
        self.label_effPulseEnergyValue.setText("{0:.3f} +- {1:.3f}".format(
            *eff_pulse_energy*1e6))  # in units of µJ

        eff_peak_intensity = self.data_processor.compute_effective_peak_intensity(
            eff_pulse_energy,
            CONSTANTS_pulse_length_FWHM,
            self.doc.w0)
        self.doc.eff_peak_intensity = eff_peak_intensity
        self.label_peakIntensityEffectiveValue.setText("{0:.0f} +- {1:.0f}".format(
            *eff_peak_intensity*1e-10))  # in units of MW/cm^2

        if self.labelApertureTransmittanceValue.text() == "":
            self.groupBox_Aperture.setEnabled(True)
        else:
            self.groupBox_Measurement.setEnabled(False)  # redo the aperture measurement!
            self.labelApertureTransmittanceValue.clear()


    def onClick_measureApertureTransmission(self):
        signals = self.nidaq_reader.get_nidaq_measurement_max_values(
            self.spinBox_iterations.value())

        S = self.data_processor.extract_aperture_transmission(signals[0], signals[2])
        self.doc.S = S
        self.labelApertureTransmittanceValue.setText("{0:.3f} +- {1:.3f}".format(*list(S)))

        self.groupBox_Measurement.setEnabled(True)


    def onClick_startMeasurement(self):
        """ Assumption: Stages are located at their initial position! """

        if self.assert_measurement_ready_to_start() != "everything good to go!":
            return None
        

        # abbreviation
        tot_num_of_pos = self.spinBox_numPositions.value()

        """
        # Possibilitly one of two
        # Move both stages in tiny steps:
        for pos_index in range(tot_num_of_pos):
            signals = self.nidaq_reader.get_nidaq_measurement_max_values(
            self.spinBox_iterations.value())
            
            # Position with respect to beam:
            # If the physical stage position is zero, it is actually behind the focal spot (this is
            # because the stages are aligned such that their max position is in front and their zero
            # position behind the focal spot). Hence, we invert the positions now:
            position_wrt_beam = self.stage_controller.total_travel_distance - \
                                    self.stage_controller.combined_position

            self.data_processor.extract_oa_ca_transmissions(position_wrt_beam, *signals,
                tot_num_of_pos)
            # Don't move the last time because stages are already at their maximum positions:
            if pos_index < tot_num_of_pos-1:
                self.stage_controller.move_in_steps(tot_num_of_pos, "backward")
        """
        
        # Possibility two of two
        # Firstly, move the first stage, and only if necessary the second stage and so on:
        for position in np.linspace(self.stage_controller.total_travel_distance, 0, tot_num_of_pos):
            self.stage_controller.move_to_position(position)
            signals = self.nidaq_reader.get_nidaq_measurement_max_values(
            self.spinBox_iterations.value())
            position_wrt_beam = self.stage_controller.total_travel_distance - \
                                    self.stage_controller.combined_position
            self.data_processor.extract_oa_ca_transmissions(position_wrt_beam, *signals,
                tot_num_of_pos)


        storage_directory, folder_num = self.doc.get_new_directory_for_storage()

        data_file_header = self.doc.get_data_file_header()
        self.data_processor.store_transmission_data(storage_directory, folder_num, data_file_header)
        
        self.stage_controller.initialise_stages()


        # data postprocessing:
        data_analyser = data_evaluation.zScanDataAnalyser(self.data_processor.T_OA,
                                                          self.data_processor.T_CA,
                                                          self.doc)

        if self.checkBox_wantFit.isChecked():
            data_analyser.perform_combined_fit()
            data_analyser.perform_independent_ca_fit()
            data_analyser.perform_ca_normalised_wrt_oa_fit()
            data_analyser.perform_ca_reconstructed_fit()
            data_analyser.store_fit_results(storage_directory, folder_num)
        data_analyser.plot_data(storage_directory, folder_num)
        data_analyser.reinitialise()

        self.data_processor.reinitialise()


    def assert_measurement_ready_to_start(self):
        # Assert that stages are either at their home position or at their maximum position with a
        # precision of 0.02mm=20µm.
        assert isclose(self.stage_controller.combined_position, 0, abs_tol=0.02) or \
        isclose(self.stage_controller.combined_position,
            self.stage_controller.total_travel_distance,
            abs_tol=0.02)

        # Make sure a sample has been specified, otherwise return to the main loop.
        if self.lineEdit_sample.text() == "--": # default String
            QtWidgets.QMessageBox.information(self, "Specify sample",
            "Please specify the sample under examination.")
            return None


        if self.doc.S[0] > 0.40:
            buttonReply = QtWidgets.QMessageBox.question(self, "Large aperture",
                "The aperture transmission is larger than S = 40%. Do you still wish to start " + \
                "the measurement?")
            if buttonReply == QtWidgets.QMessageBox.No:
                return None


        # Assert that documentation is complete
        self.doc.assert_completeness()

        return "everything good to go!"


    def exit_handler(self):
        # Persist last settings
        last_settings = last_gui_settings.LastGuiSettings(
            self.lineEdit_sample.text(),
            self.lineEdit_solvent.text(),
            self.lineEdit_concentration.text(),
            self.spinBox_laserRepRate.value(),
            self.doubleSpinBox_geomSampleLength.value(),
            self.doubleSpinBox_refrIndexSample.value(),
            self.doubleSpinBox_refrIndexAmbient.value(),
            self.lineEdit_furtherNotes.text(),
            self.spinBox_numPositions.value(),
            self.spinBox_samplingRate.value(),
            self.spinBox_samplesPerChannel.value(),
            self.spinBox_iterations.value(),
            self.doubleSpinBox_II0.value(),
            self.label_alphaValue.text(),
            self.label_effSampleLengthValue.text(),
            self.doubleSpinBox_attenuationPdRef.value())

        last_gui_settings.persist_last_settings(last_settings)
        sys.exit(0)



if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    window = Window()
    app.aboutToQuit.connect(window.exit_handler)
    window.show()
    app.exec_()
    