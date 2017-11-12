import sys
import numpy as np
from PyQt5 import QtWidgets
import gui_design
import data_analysis
import stage_control
import nidaq_control



# TODO:
# -----------------------------------
# - Check how the fits would perform with another zR. Check how they would perform if I deleted the
#   "middle points."
# - Read the paper sent by Martin.
# - A fit function to convert the reference photodiode signal to the input pulse energy is necessary.
#   It should also take the pulse repetition rate into account! For high pulse rep rates, the photo
#   diode signals increase, however, the power actually decreases! Make this calibration measurement!
#   A change of the pulse rep rate should erase the calibration and aperture values that might have
#   been measured before changing the pulse rep rate.
# - The absorption (alpha coefficient) measurement should be taken care of.
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
            self.lineEdit_furtherNotes.text())


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

        self.lineEdit_sampleMaterial.textChanged.connect(self.onNotesChange)
        self.lineEdit_solvent.textChanged.connect(self.onNotesChange)
        self.doubleSpinBox_concentration.valueChanged.connect(self.onNotesChange)
        self.spinBox_laserRepRate.valueChanged.connect(self.onNotesChange)
        self.spinBox_numPositions.valueChanged.connect(self.onNotesChange)
        self.lineEdit_furtherNotes.textChanged.connect(self.onNotesChange)

        self.spinBox_samplingRate.valueChanged.connect(self.onNidaqParamsChange)
        self.spinBox_samplesPerChannel.valueChanged.connect(self.onNidaqParamsChange)
        self.spinBox_iterations.valueChanged.connect(self.onNidaqParamsChange)


    def onNotesChange(self):
        self.data_analyser.sample_material = self.lineEdit_sampleMaterial.text()
        self.data_analyser.solvent = self.lineEdit_solvent.text()
        self.data_analyser.concentration = self.doubleSpinBox_concentration.value()
        self.data_analyser.spinBox_laserRepRate = self.spinBox_laserRepRate.value()
        self.data_analyser.tot_num_of_pos = self.spinBox_numPositions.value()
        self.data_analyser.furtherNotes = self.lineEdit_furtherNotes.text()


    def onNidaqParamsChange(self):
        self.nidaq_reader.sampling_rate = self.spinBox_samplingRate.value()
        self.nidaq_reader.num_samples_per_chan = self.spinBox_samplesPerChannel.value()
        self.nidaq_reader.iterations = self.spinBox_iterations.value()


    def onClick_calibratePhotodiodes(self):
        signals = self.nidaq_reader.get_nidaq_measurement_max_values()
        pulse_energy = self.data_analyser.extract_pulse_energy(signals[0])
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

        material = self.data_analyser.sample_material

        # Make sure a note on the measurement has been typed in, otherwise return to the main loop.
        if material == "--": # default String
            self.lineEdit_sampleMaterial.setStyleSheet("color: rgb(255,0,0)")
            return None
        else:
            self.lineEdit_sampleMaterial.setStyleSheet("color: rgb(0,0,0)")

        # abbreviation
        tot_num_of_pos = self.data_analyser.tot_num_of_pos

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

        self.data_analyser.evaluate_measurement_and_reinitialise()
        self.stage_controller.initialise_stages()




if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    window = Window()
    window.show()
    sys.exit(app.exec_())
