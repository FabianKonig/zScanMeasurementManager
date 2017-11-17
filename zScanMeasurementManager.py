import sys
import numpy as np
from PyQt5 import QtWidgets
import gui_design
import data_analysis
import stage_control
import nidaq_control
from math import isclose


# TODO:
# -----------------------------------
# - Condensates of Light Anmeldung.
# - Allow more decimal digits in I/I0 field in GUI.
# - When changing the geometrical length, the effective length must be changed as well in GUI

# - Make more measurements with Rhodamine with low power and high power and also with different repetition rates.
# - Try to fit Julians "5.dat" measurement of RH6G in Ethylenglykol with both curves separately.
#   Do I obtain identical (at least similar) results?
# - Make a measurement with ZnSe and Rhodamine-Ethylenglykol and water.
# - With those measurements, make sure the behaviour of Nitrobenzole is not due to alignment, but
#   really only due to Nitrobenzole itself.
# - What if you decrease the time the sample is exposed to radiation? Does the Rayleigh length decrease?
#   Then it might be a thermal effect!
# 
# - Read the paper sent by Martin!
#
# - A fit function to convert the reference photodiode signal to the input pulse energy is necessary.
#   It should also take the pulse repetition rate into account! For high pulse rep rates, the photo
#   diode signals increase, however, the power actually decreases! Make this calibration measurement!
#   A change of the pulse rep rate should erase pulse energy label that might have been measured
#   before changing the pulse rep rate.
# - The absorption (alpha coefficient) measurement should be taken care of.                             DONE. Check it.
# - When the measurement is started, the measurement parameters section should be disabled.
#   At the moment no problem as the GUI freezes anyway.
# - Take care of multithreading.


class Window(QtWidgets.QMainWindow, gui_design.Ui_MainWindow):
    def __init__(self):
        super().__init__()    # call __init__ of QtWidgets.QMainWindow
        self.setupUi(self)    # call setupUI of gui_design.Ui_MainWindow (generated with QtDesigner)
        self.defineSignalsSlots()


        self.data_analyser = data_analysis.zScanDataAnalyser(
            self.spinBox_numPositions.value(),
            self.lineEdit_sampleMaterial.text(),
            self.lineEdit_solvent.text(),
            self.doubleSpinBox_concentration.value(),
            self.spinBox_laserRepRate.value(),
            self.lineEdit_furtherNotes.text(),
            self.doubleSpinBox_refrIndexMaterial.value(),
            self.doubleSpinBox_refrIndexAmbient.value(),
            self.doubleSpinBox_geomSampleLength.value(),
            float(self.label_alphaValue.text()))


        self.nidaq_reader = nidaq_control.NidaqReader(
            self.spinBox_samplingRate.value(), 
            self.spinBox_samplesPerChannel.value(), 
            self.spinBox_iterations.value())

        QtWidgets.QMessageBox.information(self, "Stage position initialisation",
            "The stages will now be moved home and subsequently to their initial positions." +
            " Make sure the stages can move unhindered!")
        
        self.stage_controller = stage_control.APT_Controller()



    def defineSignalsSlots(self):
        self.pushButton_calibratePDs.clicked.connect(self.onClick_calibratePhotodiodes)
        self.pushButton_measureAperture.clicked.connect(self.onClick_measureApertureTransmission)
        self.pushButton_startMeasurement.clicked.connect(self.onClick_startMeasurement)
        self.pushButton_computeAlpha.clicked.connect(self.onClick_computeAlpha)

        self.lineEdit_sampleMaterial.textChanged.connect(self.onNotesChange)
        self.lineEdit_solvent.textChanged.connect(self.onNotesChange)
        self.doubleSpinBox_concentration.valueChanged.connect(self.onNotesChange)
        self.spinBox_laserRepRate.valueChanged.connect(self.onNotesChange)
        self.spinBox_numPositions.valueChanged.connect(self.onNotesChange)
        self.lineEdit_furtherNotes.textChanged.connect(self.onNotesChange)
        self.doubleSpinBox_geomSampleLength.valueChanged.connect(self.onNotesChange)

        self.spinBox_samplingRate.valueChanged.connect(self.onNidaqParamsChange)
        self.spinBox_samplesPerChannel.valueChanged.connect(self.onNidaqParamsChange)
        self.spinBox_iterations.valueChanged.connect(self.onNidaqParamsChange)

        self.doubleSpinBox_refrIndexMaterial.valueChanged.connect(self.onNotesChange)
        self.doubleSpinBox_refrIndexAmbient.valueChanged.connect(self.onNotesChange)


    def onNotesChange(self):
        self.data_analyser.sample_material = self.lineEdit_sampleMaterial.text()
        self.data_analyser.solvent = self.lineEdit_solvent.text()
        self.data_analyser.concentration = self.doubleSpinBox_concentration.value()
        self.data_analyser.laser_rep_rate = self.spinBox_laserRepRate.value()
        self.data_analyser.tot_num_of_pos = self.spinBox_numPositions.value()
        self.data_analyser.furtherNotes = self.lineEdit_furtherNotes.text()
        self.data_analyser.refr_index_material = self.doubleSpinBox_refrIndexMaterial.value()
        self.data_analyser.refr_index_ambient = self.doubleSpinBox_refrIndexAmbient.value()
        self.data_analyser.geom_length = self.doubleSpinBox_geomSampleLength.value()


    def onNidaqParamsChange(self):
        self.nidaq_reader.sampling_rate = self.spinBox_samplingRate.value()
        self.nidaq_reader.num_samples_per_chan = self.spinBox_samplesPerChannel.value()
        self.nidaq_reader.iterations = self.spinBox_iterations.value()


    def onClick_computeAlpha(self):
        transmission = self.doubleSpinBox_II0.value()
        refr_index_material = self.doubleSpinBox_refrIndexMaterial.value()
        refr_index_ambient = self.doubleSpinBox_refrIndexAmbient.value()
        geom_length = self.doubleSpinBox_geomSampleLength.value()  # given in mm

        transmission_glass_air = 1 - ((refr_index_ambient-1)/(refr_index_ambient+1))**2
        transmission_ambient_material = 1 - ((refr_index_ambient-refr_index_material)/(refr_index_ambient+refr_index_material))**2
        # each transmission is squared because the light has to transmit each boundary twice:
        expected_transmission = transmission_ambient_material**2 * transmission_glass_air**2

        # the transmission through the medium (after being corrected by Fresnel transmission)
        transmission_medium = transmission / expected_transmission

        alpha = -np.log(transmission_medium) / geom_length  # in 1/mm

        self.label_alphaValue.setText("{0:.3f}".format(alpha))
        self.data_analyser.alpha = alpha
        eff_length = self.data_analyser.compute_effective_length() * 1e3  # in mm
        self.label_effLengthValue.setText("{0:.3f}".format(eff_length))


    def onClick_calibratePhotodiodes(self):
        signals = self.nidaq_reader.get_nidaq_measurement_max_values()
        pulse_energy = self.data_analyser.extract_pulse_energy(
            signals[0] * self.doubleSpinBox_attenuationPdRef.value())
        self.label_pulseEnergyValue.setText("{0:.3f} +- {1:.3f}".format(*pulse_energy*1e6))

        calib_factors = list(self.data_analyser.extract_calibration_factors(*signals))
        self.label_cOAValue.setText("{0:.3f} +- {1:.3f}".format(*calib_factors[0]))
        self.label_cCAValue.setText("{0:.3f} +- {1:.3f}".format(*calib_factors[1]))


        if self.labelApertureTransmittanceValue.text() == "":
            self.groupBox_Aperture.setEnabled(True)
        else:
            self.groupBox_Measurement.setEnabled(False)  # redo the aperture measurement!
            self.labelApertureTransmittanceValue.clear()


    def onClick_measureApertureTransmission(self):
        signals = self.nidaq_reader.get_nidaq_measurement_max_values()

        s = list(self.data_analyser.extract_aperture_transmission(signals[0], signals[2]))
        self.labelApertureTransmittanceValue.setText("{0:.3f} +- {1:.3f}".format(*s))

        self.groupBox_Measurement.setEnabled(True)


    def onClick_startMeasurement(self):
        """ Assumption: Stages are located at their initial position!
        """
        # Assert that stages are either at their home position or at their maximum position with a
        # precision of 0.02mm=20µm.
        assert isclose(self.stage_controller.combined_position, 0, abs_tol=0.02) or \
        isclose(self.stage_controller.combined_position,
            self.stage_controller.total_travel_distance,
            abs_tol=0.02)


        # Make sure the material of the measurement has been typed in, otherwise return to the main loop.
        if self.data_analyser.sample_material == "--": # default String
            QtWidgets.QMessageBox.information(self, "Specify material",
            "Please specify the sample material.")
            return None

        # abbreviation
        tot_num_of_pos = self.data_analyser.tot_num_of_pos

        """
        # Possibilitly one of two
        # Move both stages in tiny steps:
        for pos_index in range(tot_num_of_pos):
            signals = self.nidaq_reader.get_nidaq_measurement_max_values()
            
            # Position with respect to beam:
            # If the physical stage position is zero, it is actually behind the focal spot (this is
            # because the stages are aligned such that their max position is in front and their zero
            # position behind the focal spot). Hence, we invert the positions now:
            position_wrt_beam = self.stage_controller.total_travel_distance - \
                                    self.stage_controller.combined_position

            self.data_analyser.extract_oa_ca_transmissions(position_wrt_beam, *signals)
            # Don't move the last time because stages are already at their maximum positions:
            if pos_index < tot_num_of_pos-1:
                self.stage_controller.move_in_steps(tot_num_of_pos, "backward")
        """
        
        # Possibilitly two of two
        # Firstly, move the first stage, and only if necessary the second stage and so on:
        for position in np.linspace(self.stage_controller.total_travel_distance, 0, tot_num_of_pos):
            self.stage_controller.move_to_position(position)
            signals = self.nidaq_reader.get_nidaq_measurement_max_values()
            position_wrt_beam = self.stage_controller.total_travel_distance - \
                                    self.stage_controller.combined_position
            self.data_analyser.extract_oa_ca_transmissions(position_wrt_beam, *signals)
        

        self.data_analyser.evaluate_measurement_and_reinitialise(self.checkBox_wantFit.isChecked())
        self.stage_controller.initialise_stages()




if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    window = Window()
    window.show()
    sys.exit(app.exec_())
