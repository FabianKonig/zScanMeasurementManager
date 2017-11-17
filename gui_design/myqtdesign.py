# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'design.ui'
#
# Created by: PyQt5 UI code generator 5.9
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(597, 802)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MainWindow.sizePolicy().hasHeightForWidth())
        MainWindow.setSizePolicy(sizePolicy)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.centralwidget)
        self.verticalLayout.setObjectName("verticalLayout")
        self.groupBox4Notes = QtWidgets.QGroupBox(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox4Notes.sizePolicy().hasHeightForWidth())
        self.groupBox4Notes.setSizePolicy(sizePolicy)
        self.groupBox4Notes.setObjectName("groupBox4Notes")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.groupBox4Notes)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.gridLayout_measurementParams = QtWidgets.QGridLayout()
        self.gridLayout_measurementParams.setObjectName("gridLayout_measurementParams")
        self.spinBox_numPositions = QtWidgets.QSpinBox(self.groupBox4Notes)
        self.spinBox_numPositions.setMinimum(2)
        self.spinBox_numPositions.setMaximum(150)
        self.spinBox_numPositions.setProperty("value", 45)
        self.spinBox_numPositions.setObjectName("spinBox_numPositions")
        self.gridLayout_measurementParams.addWidget(self.spinBox_numPositions, 6, 1, 1, 1)
        self.lineEdit_sampleMaterial = QtWidgets.QLineEdit(self.groupBox4Notes)
        self.lineEdit_sampleMaterial.setObjectName("lineEdit_sampleMaterial")
        self.gridLayout_measurementParams.addWidget(self.lineEdit_sampleMaterial, 0, 1, 1, 1)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout_measurementParams.addItem(spacerItem, 1, 2, 1, 1)
        self.label_solvent = QtWidgets.QLabel(self.groupBox4Notes)
        self.label_solvent.setObjectName("label_solvent")
        self.gridLayout_measurementParams.addWidget(self.label_solvent, 0, 3, 1, 1)
        self.label_sampleMaterial = QtWidgets.QLabel(self.groupBox4Notes)
        self.label_sampleMaterial.setObjectName("label_sampleMaterial")
        self.gridLayout_measurementParams.addWidget(self.label_sampleMaterial, 0, 0, 1, 1)
        self.label_laserRepRate = QtWidgets.QLabel(self.groupBox4Notes)
        self.label_laserRepRate.setObjectName("label_laserRepRate")
        self.gridLayout_measurementParams.addWidget(self.label_laserRepRate, 1, 3, 1, 1)
        self.lineEdit_solvent = QtWidgets.QLineEdit(self.groupBox4Notes)
        self.lineEdit_solvent.setObjectName("lineEdit_solvent")
        self.gridLayout_measurementParams.addWidget(self.lineEdit_solvent, 0, 4, 1, 1)
        self.label_concentration = QtWidgets.QLabel(self.groupBox4Notes)
        self.label_concentration.setObjectName("label_concentration")
        self.gridLayout_measurementParams.addWidget(self.label_concentration, 1, 0, 1, 1)
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout_measurementParams.addItem(spacerItem1, 0, 2, 1, 1)
        self.doubleSpinBox_concentration = QtWidgets.QDoubleSpinBox(self.groupBox4Notes)
        self.doubleSpinBox_concentration.setDecimals(2)
        self.doubleSpinBox_concentration.setMaximum(10000.0)
        self.doubleSpinBox_concentration.setObjectName("doubleSpinBox_concentration")
        self.gridLayout_measurementParams.addWidget(self.doubleSpinBox_concentration, 1, 1, 1, 1)
        self.spinBox_laserRepRate = QtWidgets.QSpinBox(self.groupBox4Notes)
        self.spinBox_laserRepRate.setMinimum(10)
        self.spinBox_laserRepRate.setMaximum(1100)
        self.spinBox_laserRepRate.setSingleStep(10)
        self.spinBox_laserRepRate.setObjectName("spinBox_laserRepRate")
        self.gridLayout_measurementParams.addWidget(self.spinBox_laserRepRate, 1, 4, 1, 1)
        spacerItem2 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout_measurementParams.addItem(spacerItem2, 6, 2, 1, 1)
        self.doubleSpinBox_refrIndexMaterial = QtWidgets.QDoubleSpinBox(self.groupBox4Notes)
        self.doubleSpinBox_refrIndexMaterial.setDecimals(2)
        self.doubleSpinBox_refrIndexMaterial.setMinimum(1.0)
        self.doubleSpinBox_refrIndexMaterial.setMaximum(100.0)
        self.doubleSpinBox_refrIndexMaterial.setSingleStep(0.1)
        self.doubleSpinBox_refrIndexMaterial.setProperty("value", 1.46)
        self.doubleSpinBox_refrIndexMaterial.setObjectName("doubleSpinBox_refrIndexMaterial")
        self.gridLayout_measurementParams.addWidget(self.doubleSpinBox_refrIndexMaterial, 4, 1, 1, 1)
        self.label_numPositions = QtWidgets.QLabel(self.groupBox4Notes)
        self.label_numPositions.setObjectName("label_numPositions")
        self.gridLayout_measurementParams.addWidget(self.label_numPositions, 6, 0, 1, 1)
        self.labe_geomSampleLength = QtWidgets.QLabel(self.groupBox4Notes)
        self.labe_geomSampleLength.setObjectName("labe_geomSampleLength")
        self.gridLayout_measurementParams.addWidget(self.labe_geomSampleLength, 2, 0, 1, 1)
        self.doubleSpinBox_geomSampleLength = QtWidgets.QDoubleSpinBox(self.groupBox4Notes)
        self.doubleSpinBox_geomSampleLength.setSingleStep(0.1)
        self.doubleSpinBox_geomSampleLength.setProperty("value", 1.0)
        self.doubleSpinBox_geomSampleLength.setObjectName("doubleSpinBox_geomSampleLength")
        self.gridLayout_measurementParams.addWidget(self.doubleSpinBox_geomSampleLength, 2, 1, 1, 1)
        self.doubleSpinBox_refrIndexAmbient = QtWidgets.QDoubleSpinBox(self.groupBox4Notes)
        self.doubleSpinBox_refrIndexAmbient.setDecimals(2)
        self.doubleSpinBox_refrIndexAmbient.setMinimum(1.0)
        self.doubleSpinBox_refrIndexAmbient.setMaximum(100.0)
        self.doubleSpinBox_refrIndexAmbient.setSingleStep(0.1)
        self.doubleSpinBox_refrIndexAmbient.setProperty("value", 1.46)
        self.doubleSpinBox_refrIndexAmbient.setObjectName("doubleSpinBox_refrIndexAmbient")
        self.gridLayout_measurementParams.addWidget(self.doubleSpinBox_refrIndexAmbient, 4, 4, 1, 1)
        self.label_RefrIndexAmbient = QtWidgets.QLabel(self.groupBox4Notes)
        self.label_RefrIndexAmbient.setObjectName("label_RefrIndexAmbient")
        self.gridLayout_measurementParams.addWidget(self.label_RefrIndexAmbient, 4, 3, 1, 1)
        spacerItem3 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout_measurementParams.addItem(spacerItem3, 4, 2, 1, 1)
        self.label_RefrIndexMaterial = QtWidgets.QLabel(self.groupBox4Notes)
        self.label_RefrIndexMaterial.setObjectName("label_RefrIndexMaterial")
        self.gridLayout_measurementParams.addWidget(self.label_RefrIndexMaterial, 4, 0, 1, 1)
        self.label_furtherNotes = QtWidgets.QLabel(self.groupBox4Notes)
        self.label_furtherNotes.setObjectName("label_furtherNotes")
        self.gridLayout_measurementParams.addWidget(self.label_furtherNotes, 5, 0, 1, 1)
        self.lineEdit_furtherNotes = QtWidgets.QLineEdit(self.groupBox4Notes)
        self.lineEdit_furtherNotes.setObjectName("lineEdit_furtherNotes")
        self.gridLayout_measurementParams.addWidget(self.lineEdit_furtherNotes, 5, 1, 1, 4)
        self.verticalLayout_4.addLayout(self.gridLayout_measurementParams)
        self.verticalLayout.addWidget(self.groupBox4Notes)
        self.line_1 = QtWidgets.QFrame(self.centralwidget)
        self.line_1.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_1.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_1.setObjectName("line_1")
        self.verticalLayout.addWidget(self.line_1)
        self.groupBox_nidaqParams = QtWidgets.QGroupBox(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_nidaqParams.sizePolicy().hasHeightForWidth())
        self.groupBox_nidaqParams.setSizePolicy(sizePolicy)
        self.groupBox_nidaqParams.setObjectName("groupBox_nidaqParams")
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout(self.groupBox_nidaqParams)
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.label_samplingRate = QtWidgets.QLabel(self.groupBox_nidaqParams)
        self.label_samplingRate.setObjectName("label_samplingRate")
        self.horizontalLayout_3.addWidget(self.label_samplingRate)
        self.spinBox_samplingRate = QtWidgets.QSpinBox(self.groupBox_nidaqParams)
        self.spinBox_samplingRate.setMinimum(1000)
        self.spinBox_samplingRate.setMaximum(250000)
        self.spinBox_samplingRate.setSingleStep(10000)
        self.spinBox_samplingRate.setProperty("value", 250000)
        self.spinBox_samplingRate.setObjectName("spinBox_samplingRate")
        self.horizontalLayout_3.addWidget(self.spinBox_samplingRate)
        spacerItem4 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem4)
        self.label_samplesPerChannel = QtWidgets.QLabel(self.groupBox_nidaqParams)
        self.label_samplesPerChannel.setObjectName("label_samplesPerChannel")
        self.horizontalLayout_3.addWidget(self.label_samplesPerChannel)
        self.spinBox_samplesPerChannel = QtWidgets.QSpinBox(self.groupBox_nidaqParams)
        self.spinBox_samplesPerChannel.setMinimum(1000)
        self.spinBox_samplesPerChannel.setMaximum(500000)
        self.spinBox_samplesPerChannel.setSingleStep(1000)
        self.spinBox_samplesPerChannel.setProperty("value", 70000)
        self.spinBox_samplesPerChannel.setObjectName("spinBox_samplesPerChannel")
        self.horizontalLayout_3.addWidget(self.spinBox_samplesPerChannel)
        spacerItem5 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem5)
        self.label_iterations = QtWidgets.QLabel(self.groupBox_nidaqParams)
        self.label_iterations.setObjectName("label_iterations")
        self.horizontalLayout_3.addWidget(self.label_iterations)
        self.spinBox_iterations = QtWidgets.QSpinBox(self.groupBox_nidaqParams)
        self.spinBox_iterations.setMinimum(2)
        self.spinBox_iterations.setProperty("value", 4)
        self.spinBox_iterations.setObjectName("spinBox_iterations")
        self.horizontalLayout_3.addWidget(self.spinBox_iterations)
        self.verticalLayout.addWidget(self.groupBox_nidaqParams)
        self.line = QtWidgets.QFrame(self.centralwidget)
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.verticalLayout.addWidget(self.line)
        self.groupBox = QtWidgets.QGroupBox(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox.sizePolicy().hasHeightForWidth())
        self.groupBox.setSizePolicy(sizePolicy)
        self.groupBox.setObjectName("groupBox")
        self.verticalLayout_5 = QtWidgets.QVBoxLayout(self.groupBox)
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.gridLayout_2 = QtWidgets.QGridLayout()
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.pushButton_computeAlpha = QtWidgets.QPushButton(self.groupBox)
        self.pushButton_computeAlpha.setObjectName("pushButton_computeAlpha")
        self.gridLayout_2.addWidget(self.pushButton_computeAlpha, 1, 0, 1, 2)
        self.label_II0 = QtWidgets.QLabel(self.groupBox)
        self.label_II0.setObjectName("label_II0")
        self.gridLayout_2.addWidget(self.label_II0, 0, 0, 1, 1)
        self.doubleSpinBox_II0 = QtWidgets.QDoubleSpinBox(self.groupBox)
        self.doubleSpinBox_II0.setProperty("value", 1.0)
        self.doubleSpinBox_II0.setObjectName("doubleSpinBox_II0")
        self.gridLayout_2.addWidget(self.doubleSpinBox_II0, 0, 1, 1, 1)
        self.line_4 = QtWidgets.QFrame(self.groupBox)
        self.line_4.setFrameShape(QtWidgets.QFrame.VLine)
        self.line_4.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_4.setObjectName("line_4")
        self.gridLayout_2.addWidget(self.line_4, 0, 2, 2, 1)
        self.label_alphaValue = QtWidgets.QLabel(self.groupBox)
        self.label_alphaValue.setObjectName("label_alphaValue")
        self.gridLayout_2.addWidget(self.label_alphaValue, 1, 4, 1, 1)
        self.label_alpha = QtWidgets.QLabel(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_alpha.sizePolicy().hasHeightForWidth())
        self.label_alpha.setSizePolicy(sizePolicy)
        self.label_alpha.setMaximumSize(QtCore.QSize(65, 16777215))
        self.label_alpha.setObjectName("label_alpha")
        self.gridLayout_2.addWidget(self.label_alpha, 1, 3, 1, 1)
        self.verticalLayout_5.addLayout(self.gridLayout_2)
        self.verticalLayout.addWidget(self.groupBox)
        self.groupBox_Calibration = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox_Calibration.setEnabled(True)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_Calibration.sizePolicy().hasHeightForWidth())
        self.groupBox_Calibration.setSizePolicy(sizePolicy)
        self.groupBox_Calibration.setObjectName("groupBox_Calibration")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.groupBox_Calibration)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.widgetRefractiveIndexCalibrateBtn = QtWidgets.QWidget(self.groupBox_Calibration)
        self.widgetRefractiveIndexCalibrateBtn.setObjectName("widgetRefractiveIndexCalibrateBtn")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.widgetRefractiveIndexCalibrateBtn)
        self.verticalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.gridLayout_refrIndexes = QtWidgets.QGridLayout()
        self.gridLayout_refrIndexes.setObjectName("gridLayout_refrIndexes")
        self.label_attenuationPdRef = QtWidgets.QLabel(self.widgetRefractiveIndexCalibrateBtn)
        self.label_attenuationPdRef.setObjectName("label_attenuationPdRef")
        self.gridLayout_refrIndexes.addWidget(self.label_attenuationPdRef, 2, 0, 1, 1)
        self.doubleSpinBox_attenuationPdRef = QtWidgets.QDoubleSpinBox(self.widgetRefractiveIndexCalibrateBtn)
        self.doubleSpinBox_attenuationPdRef.setMinimum(1.0)
        self.doubleSpinBox_attenuationPdRef.setMaximum(10000.0)
        self.doubleSpinBox_attenuationPdRef.setSingleStep(100.0)
        self.doubleSpinBox_attenuationPdRef.setProperty("value", 1000.0)
        self.doubleSpinBox_attenuationPdRef.setObjectName("doubleSpinBox_attenuationPdRef")
        self.gridLayout_refrIndexes.addWidget(self.doubleSpinBox_attenuationPdRef, 2, 1, 1, 1)
        self.verticalLayout_2.addLayout(self.gridLayout_refrIndexes)
        self.pushButton_calibratePDs = QtWidgets.QPushButton(self.widgetRefractiveIndexCalibrateBtn)
        self.pushButton_calibratePDs.setObjectName("pushButton_calibratePDs")
        self.verticalLayout_2.addWidget(self.pushButton_calibratePDs)
        self.horizontalLayout_2.addWidget(self.widgetRefractiveIndexCalibrateBtn)
        spacerItem6 = QtWidgets.QSpacerItem(2, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem6)
        self.line_2 = QtWidgets.QFrame(self.groupBox_Calibration)
        self.line_2.setFrameShape(QtWidgets.QFrame.VLine)
        self.line_2.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_2.setObjectName("line_2")
        self.horizontalLayout_2.addWidget(self.line_2)
        self.widgetCalibrationValues = QtWidgets.QWidget(self.groupBox_Calibration)
        self.widgetCalibrationValues.setObjectName("widgetCalibrationValues")
        self.gridLayout = QtWidgets.QGridLayout(self.widgetCalibrationValues)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")
        self.label_cCA = QtWidgets.QLabel(self.widgetCalibrationValues)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_cCA.sizePolicy().hasHeightForWidth())
        self.label_cCA.setSizePolicy(sizePolicy)
        self.label_cCA.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_cCA.setObjectName("label_cCA")
        self.gridLayout.addWidget(self.label_cCA, 0, 1, 1, 1)
        self.label_cOA = QtWidgets.QLabel(self.widgetCalibrationValues)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_cOA.sizePolicy().hasHeightForWidth())
        self.label_cOA.setSizePolicy(sizePolicy)
        self.label_cOA.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_cOA.setObjectName("label_cOA")
        self.gridLayout.addWidget(self.label_cOA, 1, 1, 1, 1)
        self.label_cOAValue = QtWidgets.QLabel(self.widgetCalibrationValues)
        self.label_cOAValue.setText("")
        self.label_cOAValue.setObjectName("label_cOAValue")
        self.gridLayout.addWidget(self.label_cOAValue, 1, 2, 1, 1)
        self.label_cCAValue = QtWidgets.QLabel(self.widgetCalibrationValues)
        self.label_cCAValue.setText("")
        self.label_cCAValue.setObjectName("label_cCAValue")
        self.gridLayout.addWidget(self.label_cCAValue, 0, 2, 1, 1)
        self.label_pulseEnergy = QtWidgets.QLabel(self.widgetCalibrationValues)
        self.label_pulseEnergy.setObjectName("label_pulseEnergy")
        self.gridLayout.addWidget(self.label_pulseEnergy, 2, 1, 1, 1)
        self.label_pulseEnergyValue = QtWidgets.QLabel(self.widgetCalibrationValues)
        self.label_pulseEnergyValue.setText("")
        self.label_pulseEnergyValue.setObjectName("label_pulseEnergyValue")
        self.gridLayout.addWidget(self.label_pulseEnergyValue, 2, 2, 1, 1)
        self.horizontalLayout_2.addWidget(self.widgetCalibrationValues)
        self.verticalLayout.addWidget(self.groupBox_Calibration)
        self.groupBox_Aperture = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox_Aperture.setEnabled(False)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_Aperture.sizePolicy().hasHeightForWidth())
        self.groupBox_Aperture.setSizePolicy(sizePolicy)
        self.groupBox_Aperture.setObjectName("groupBox_Aperture")
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout(self.groupBox_Aperture)
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        self.pushButton_measureAperture = QtWidgets.QPushButton(self.groupBox_Aperture)
        self.pushButton_measureAperture.setObjectName("pushButton_measureAperture")
        self.horizontalLayout_4.addWidget(self.pushButton_measureAperture)
        self.labelApertureTransmittance = QtWidgets.QLabel(self.groupBox_Aperture)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.labelApertureTransmittance.sizePolicy().hasHeightForWidth())
        self.labelApertureTransmittance.setSizePolicy(sizePolicy)
        self.labelApertureTransmittance.setObjectName("labelApertureTransmittance")
        self.horizontalLayout_4.addWidget(self.labelApertureTransmittance)
        self.labelApertureTransmittanceValue = QtWidgets.QLabel(self.groupBox_Aperture)
        self.labelApertureTransmittanceValue.setText("")
        self.labelApertureTransmittanceValue.setObjectName("labelApertureTransmittanceValue")
        self.horizontalLayout_4.addWidget(self.labelApertureTransmittanceValue)
        spacerItem7 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_4.addItem(spacerItem7)
        self.verticalLayout.addWidget(self.groupBox_Aperture)
        self.line_3 = QtWidgets.QFrame(self.centralwidget)
        self.line_3.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_3.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_3.setObjectName("line_3")
        self.verticalLayout.addWidget(self.line_3)
        self.groupBox_Measurement = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox_Measurement.setEnabled(False)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_Measurement.sizePolicy().hasHeightForWidth())
        self.groupBox_Measurement.setSizePolicy(sizePolicy)
        self.groupBox_Measurement.setObjectName("groupBox_Measurement")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.groupBox_Measurement)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.gridLayout_startMeasurement = QtWidgets.QGridLayout()
        self.gridLayout_startMeasurement.setObjectName("gridLayout_startMeasurement")
        self.progressBarMeasurement = QtWidgets.QProgressBar(self.groupBox_Measurement)
        self.progressBarMeasurement.setProperty("value", 0)
        self.progressBarMeasurement.setObjectName("progressBarMeasurement")
        self.gridLayout_startMeasurement.addWidget(self.progressBarMeasurement, 0, 1, 1, 1)
        self.pushButton_startMeasurement = QtWidgets.QPushButton(self.groupBox_Measurement)
        self.pushButton_startMeasurement.setObjectName("pushButton_startMeasurement")
        self.gridLayout_startMeasurement.addWidget(self.pushButton_startMeasurement, 0, 0, 1, 1)
        self.checkBox_wantFit = QtWidgets.QCheckBox(self.groupBox_Measurement)
        self.checkBox_wantFit.setChecked(True)
        self.checkBox_wantFit.setObjectName("checkBox_wantFit")
        self.gridLayout_startMeasurement.addWidget(self.checkBox_wantFit, 1, 0, 1, 1)
        self.verticalLayout_3.addLayout(self.gridLayout_startMeasurement)
        self.verticalLayout.addWidget(self.groupBox_Measurement)
        spacerItem8 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem8)
        MainWindow.setCentralWidget(self.centralwidget)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "z-Scan Measurement Manager"))
        self.groupBox4Notes.setTitle(_translate("MainWindow", "Notes on the measurement: "))
        self.lineEdit_sampleMaterial.setText(_translate("MainWindow", "--"))
        self.label_solvent.setText(_translate("MainWindow", "Solvent"))
        self.label_sampleMaterial.setText(_translate("MainWindow", "Sample material"))
        self.label_laserRepRate.setText(_translate("MainWindow", "Laser rep. / Hz"))
        self.lineEdit_solvent.setText(_translate("MainWindow", "--"))
        self.label_concentration.setText(_translate("MainWindow", "Concentration / mmol/l"))
        self.label_numPositions.setText(_translate("MainWindow", "No. of positions"))
        self.labe_geomSampleLength.setText(_translate("MainWindow", "Geom. sample length / mm"))
        self.label_RefrIndexAmbient.setText(_translate("MainWindow", "Ambient refractive index"))
        self.label_RefrIndexMaterial.setText(_translate("MainWindow", "Refractive index material"))
        self.label_furtherNotes.setText(_translate("MainWindow", "Further notes"))
        self.lineEdit_furtherNotes.setText(_translate("MainWindow", "---"))
        self.groupBox_nidaqParams.setTitle(_translate("MainWindow", "Nidaq acquisition parameters"))
        self.label_samplingRate.setText(_translate("MainWindow", "Sampling rate"))
        self.label_samplesPerChannel.setText(_translate("MainWindow", "Samples/channel"))
        self.label_iterations.setText(_translate("MainWindow", "Iterations"))
        self.groupBox.setTitle(_translate("MainWindow", "1. Measure linear absorption coeficient"))
        self.pushButton_computeAlpha.setText(_translate("MainWindow", "Compute α"))
        self.label_II0.setText(_translate("MainWindow", "<html><head/><body><p>I/I<span style=\" vertical-align:sub;\">0</span></p></body></html>"))
        self.label_alphaValue.setText(_translate("MainWindow", "0.00"))
        self.label_alpha.setText(_translate("MainWindow", "<html><head/><body><p>α / mm<span style=\" vertical-align:super;\">-1</span> =</p></body></html>"))
        self.groupBox_Calibration.setTitle(_translate("MainWindow", "2. Calibration: Insert medium and open the aperture (S=1)"))
        self.label_attenuationPdRef.setText(_translate("MainWindow", "<html><head/><body><p>Attenuation PD<span style=\" vertical-align:sub;\">Ref</span></p></body></html>"))
        self.pushButton_calibratePDs.setText(_translate("MainWindow", "Calibrate photodiodes"))
        self.label_cCA.setText(_translate("MainWindow", "c_CA ="))
        self.label_cOA.setText(_translate("MainWindow", "c_OA ="))
        self.label_pulseEnergy.setText(_translate("MainWindow", "<html><head/><body><p>E<span style=\" vertical-align:sub;\">Pulse</span> / µJ =</p></body></html>"))
        self.groupBox_Aperture.setTitle(_translate("MainWindow", "3. Align aperture"))
        self.pushButton_measureAperture.setText(_translate("MainWindow", "Measure aperture transmission"))
        self.labelApertureTransmittance.setText(_translate("MainWindow", "S ="))
        self.groupBox_Measurement.setTitle(_translate("MainWindow", "4. Start measurement"))
        self.pushButton_startMeasurement.setText(_translate("MainWindow", "Start measurement"))
        self.checkBox_wantFit.setText(_translate("MainWindow", "Perform a fit"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())

