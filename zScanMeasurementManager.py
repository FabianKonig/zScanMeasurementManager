import sys
import numpy as np
from PyQt5 import QtWidgets
import gui_design
import data_analysis
import stage_control
import nidaq_control
import matplotlib.pyplot as plt


class Window(QtWidgets.QMainWindow, gui_design.Ui_MainWindow):
    def __init__(self):
        super().__init__()    # call __init__ of QtWidgets.QMainWindow
        self.data_analyser = data_analysis.zScanDataAnalyser(tot_num_of_pos=30)
        self.stage_controller = None

        self.setupUi(self)    # call setupUI of gui_design.Ui_MainWindow (generated with QtDesigner)
        self.label_progress.setVisible(False)
        self.defineSignalsSlots()


        QtWidgets.QMessageBox.information(self, "Stage position initialisation",
            "The stages will now be moved home and subsequently to their initial positions." +
            " Make sure the stages can move unhindered!")
        
        print("self.stage_controller is being initialised, I will freeze")
        self.stage_controller = stage_control.APT_Controller()
        # Activate pushButtonCalibrate_PDs only after the stage initialisation has completed.


    def defineSignalsSlots(self):
        self.pushButtonCalibrate_PDs.clicked.connect(self.onClick_calibrate_photodiodes)
        self.pushButton_MeasureAperture.clicked.connect(self.onClick_measure_aperture_transmission)
        self.pushButtonStartStopMeasurement.clicked.connect(self.onClick_start_stop_measurement)


    def onClick_calibrate_photodiodes(self):
        print("Get Nidaq Signal, I will need some time and want the progress symbol to appear meanwhile.")
        signals = nidaq_control.get_nidaq_measurement_max_values()

        calib_factors = list(self.data_analyser.extract_calibration_factors(*signals))
        self.label_cOAValue.setText("{0:.3f} +- {1:.3f}".format(*calib_factors[0]))
        self.label_cCAValue.setText("{0:.3f} +- {1:.3f}".format(*calib_factors[1]))

        self.groupBox_Aperture.setEnabled(True)


    def onClick_measure_aperture_transmission(self):
        print("Get Nidaq Signal, I will need some time and want the progress symbol to appear meanwhile.")
        signals = nidaq_control.get_nidaq_measurement_max_values()

        s = list(self.data_analyser.extract_aperture_transmission(signals[0], signals[2]))
        self.labelApertureTransmittanceValue.setText("{0:.3f} +- {1:.3f}".format(*s))

        self.groupBox_Measurement.setEnabled(True)


    def onClick_start_stop_measurement(self):
        """
        Assumption: Stages are located at their initial position!
        """

        note = self.notesLineEdit.text()
        # Make sure a note on the measurement has been typed in, otherwise return to the main loop.
        if note == "Type in a note here (e.g. sample material and concentration)!":
            self.notesLineEdit.setStyleSheet("color: rgb(255,0,0)")
            return None

        else:
            self.notesLineEdit.setStyleSheet("color: rgb(0,0,0)")


        print("Will start measurement, will take plenty of time!")
        tot_num_of_pos = self.data_analyser.tot_num_of_pos

        for pos_index in range(tot_num_of_pos):
            signals = nidaq_control.get_nidaq_measurement_max_values()

            # Position with respect to beam:
            # If the physical stage position is zero, it is actually behind the focal spot (this is
            # because the stages are aligned such that their max position is in front and their zero
            # position behind the focal spot). Hence, we invert the positions now:
            position_wrt_beam = self.stage_controller.total_travel_distance -
                                    self.stage_controller.combined_position

            self.data_analyser.extract_oa_ca_transmissions(position_wrt_beam, *signals)
            # Don't move the last time because stages are already at their maximum positions:
            if pos_index < tot_num_of_pos-1:
                self.stage_controller.move_in_steps(tot_num_of_pos, "backward")


        self.data_analyser.store_transmission_data(note)
        self.data_analyser.plot_transmission()

        self.stage_controller.reinitialise_stages()
        self.data_analyser.reinitialise()




if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    window = Window()
    window.show()
    sys.exit(app.exec_())
