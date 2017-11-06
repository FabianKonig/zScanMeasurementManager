import sys
import numpy as np
from PyQt5 import QtWidgets
import gui_design
import data_analysis
#import stage_control
#import nidaq_control


class Window(QtWidgets.QMainWindow, gui_design.Ui_MainWindow):
    def __init__(self):
        super().__init__()    # call __init__ of QtWidgets.QMainWindow
        self.data_analyser = data_analysis.zScanDataAnalyser(tot_num_of_pos=2)
        self.stage_controller = None

        self.setupUi(self)    # call setupUI of gui_design.Ui_MainWindow (generated with QtDesigner)
        self.label_progress.setVisible(False)
        self.defineSignalsSlots()


        QtWidgets.QMessageBox.information(self, "Stage position initialisation",
            "The stages will now be moved home and subsequently to their initial positions." +
            " Make sure the stages can move unhindered!")
        
        print("self.stage_controller is being initialised, I will freeze")
        #self.stage_controller = stage_control.APT_Controller()
        # Activate pushButtonCalibrate_PDs only after the stage initialisation has completed.


    def defineSignalsSlots(self):
        self.pushButtonCalibrate_PDs.clicked.connect(self.onClick_calibrate_photodiodes)
        self.pushButton_MeasureAperture.clicked.connect(self.onClick_measure_aperture_transmission)
        self.pushButtonStartStopMeasurement.clicked.connect(self.onClick_start_stop_measurement)


    def onClick_calibrate_photodiodes(self):
        print("Get Nidaq Signal, I will need some time and want the progress symbol to appear meanwhile.")
        #signals = nidaq_control.get_filtered_nidaq_signal()
        signals = np.array([[12,8,2,10], [8,5,2,6.1], [12,8,2,10]])  # temporary

        calib_factors = list(self.data_analyser.extract_calibration_factors(*signals))
        self.label_cOAValue.setText("{0:.3f} +- {1:.3f}".format(*calib_factors[0]))
        self.label_cCAValue.setText("{0:.3f} +- {1:.3f}".format(*calib_factors[1]))

        self.groupBox_Aperture.setEnabled(True)


    def onClick_measure_aperture_transmission(self):
        print("Get Nidaq Signal, I will need some time and want the progress symbol to appear meanwhile.")
        #signals = nidaq_control.get_filtered_nidaq_signal()
        signals = np.array([[12,8,2,10], [8,5,2,6.1], [12*0.8,8*0.8,2*0.8,10*0.8]])  # temporary

        s = list(self.data_analyser.extract_aperture_transmission(signals[0], signals[2]))
        self.labelApertureTransmittanceValue.setText("{0:.3f} +- {1:.3f}".format(*s))

        self.groupBox_Measurement.setEnabled(True)


    def onClick_start_stop_measurement(self):
        """
        Assumption: Stages are located at their initial position!
        """
        print("Will start measurement, will take plenty of time!")
        tot_num_of_pos = self.data_analyser.tot_num_of_pos

        for pos_index in range(tot_num_of_pos):
            #signals = get_filtered_nidaq_signal()
            signals = np.array([[12,8,2,10], [8,5,2,6.1], [12*0.8,8*0.8,2*0.8,10*0.8]])  # temporary
            comb_pos = 3  #temporary
            #self.data_analyser.extract_oa_ca_transmissions(self.stage_controller.combined_position, *signals)
            self.data_analyser.extract_oa_ca_transmissions(comb_pos, *signals)
            # Don't move the last time because stages are already at their maximum positions:
            if pos_index < tot_num_of_pos-1:
                pass
                #self.stage_controller.move_in_steps(tot_num_of_pos, "backward")

        print(self.data_analyser.T_CA)

        #self.stage_controller.reinitialise_stages()
        self.data_analyser.reinitialise()



if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    window = Window()
    window.show()
    sys.exit(app.exec_())
