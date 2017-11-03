import sys
from PyQt5 import QtWidgets
import gui_design
#import stage_control as sc
#import nidaq_control as nc
import data_analysis as da


class Window(QtWidgets.QMainWindow, gui_design.Ui_MainWindow):
    def __init__(self):
        super().__init__()    # call __init__ of QtWidgets.QMainWindow

        self.setupUi(self)    # call setupUI of gui_design.Ui_MainWindow (generated with QtDesigner)
        self.defineSignals()

    def defineSignals(self):
        self.pushButtonCalibrate_PDs.clicked.connect(self.btnClicked)

    def btnClicked(self):
        self.label_OK.setText("JUHUUU")



if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    window = Window()
    window.show()
    sys.exit(app.exec_())
