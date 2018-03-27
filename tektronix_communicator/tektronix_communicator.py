import visa 
import numpy as np
from struct import unpack
import matplotlib.pyplot as plt




class TektronixCommunicator:

    def __init__(self):
        # Channels to be queried:
        self.channels = ["CH1", "CH2", "CH4"] # needs to be ordered
        self.channelSettings = np.empty(shape=(len(self.channels), 3)) #for each channel: ymult, yoff, yzero
        self.data = None
        self.xincr = None
        self.acqlen = None
        self.waveforms = None


    def acquireCurves(self):
        rm = visa.ResourceManager()
        scope = rm.open_resource('USB0::0x0699::0x0501::C012801::INSTR')
        scope.write('DATA:ENC RPB') # encoding
        scope.write('DATA:WIDTH 1')
        scope.write('Data:Stop 1e20') # read the entire time axis.

        # get the y-axis settings of each channel
        i = 0
        for channel in self.channels:
            scope.write("SELect:" + channel + " ON") #activate channel on oscilloscope, otherwise error
            scope.write("DATA:SOU " + channel)
            self.channelSettings[i,0] = float(scope.query('WFMPRE:YMULT?'))
            self.channelSettings[i,1] = float(scope.query('WFMPRE:YOFF?'))
            self.channelSettings[i,2] = float(scope.query('WFMPRE:YZERO?'))
            i+=1

        # set "Source" to all channels that are included in self.channels
        channelsString = ""
        firstEntry = True
        for channel in self.channels:
            if firstEntry:
                channelsString += channel
                firstEntry = False
            else:
                channelsString += ", " + channel
        
        scope.write("DATA:SOU " + channelsString)

        # get the x-axis settings
        self.xincr = float(scope.query('WFMPRE:XINCR?'))
        self.acqlen = int(scope.query("HORizontal:ACQLENGTH?")) # length of a dataset

        # Acquire the curves
        scope.write('CURVENext?')   
        data = scope.read_raw()
        self.data = np.array(unpack('%sB' % len(data), data))

        scope.close()


    def processCurves(self):
        headerlength = (len(self.data)-1) % (len(self.channels)*self.acqlen) // len(self.channels)
        time = np.linspace(0, self.xincr * (self.acqlen-1), self.acqlen)

        self.waveforms = np.empty(shape=(1+len(self.channels), self.acqlen))
        self.waveforms[0] = time

        for i in range(len(self.channels)):
            start = (1+i) * headerlength + i*self.acqlen
            stop = (1+i)*(headerlength+self.acqlen)

            wave = data[start:stop]
            self.waveforms[i+1] = (wave - self.channelSettings[i,1]) * self.channelSettings[i,0]  + self.channelSettings[i,2]


    def getSignalAfterPhotodiodeOscillations(self, waveform, times):
        """ Find the position of the photodiode signal when the oscillations have stopped. """
        TIMEOFOSCILLATIONS = 1e-6

        maxPos, = np.where(waveform==waveform.max())
        maxPos = int(maxPos[0])

        requiredSteps = TIMEOFOSCILLATIONS // self.xincr

        return waveform[maxPos+requiredSteps]
        

    def getOsciValues(self):
        self.acquireCurves()
        self.processCurves()

        max_signals = np.empty(shape=(len(self.channels)), dtype=float)
        max_signals[0] = self.waveforms[1].max()                                      # Channel 1
        max_signals[1] = self.getSignalAfterPhotodiodeOscillations(self.waveforms[2]) # Channel 2
        max_signals[2] = self.getSignalAfterPhotodiodeOscillations(self.waveforms[3]) # Channel 4

        return max_signals


    def getSignalPeaks(self, iterations):
        """ Output: 
            2-dim numpy array with the first dimension denoting the channel and the second
            dimension is iterations entries long, each storing the maximum value of that
            measurement. The order of the channels is 1: ref, 2: oa, 3: ca.
        """
        max_signal_values = np.empty(shape=(len(self.channels), iterations))

        for iteration_index in iterations:
            osciVals = self.getOsciValues()
            # iterate through self.channels from right to left [CA, OA, REF], such that the result is [REF, OA, CA].
            for channel_index in range(-1,-len(self.channels)-1,-1):
                max_signal_values[channel_index, iteration_index] = osciVals[channel_index]

        return max_signal_values





if __name__ == '__main__':

    # plot waveforms
    colors = ["orange", "blue", "green", "red"]

    scopeCommunicator = TektronixCommunicator()
    max_signals = scopeCommunicator.getOsciValues()
    sc = scopeCommunicator


    plt.plot(0, max_signals[0], marker="o", color=colors[0])
    plt.plot(0, max_signals[1], marker="o", color=colors[1])
    plt.plot(0, max_signals[2], marker="o", color=colors[2])

    for i in range(len(sc.channels)):
        plt.plot(sc.waveforms[0], sc.waveforms[i+1], color=colors[i], zorder=5-i)

    plt.grid()
    plt.show()





    #scope.write("HORizontal:MODE:RECOrdlength 2000")
    #scope.write("HORizontal:MODE:SAMPLERate {0}".format(1e9))#2000 / (10*float(scope.query("HORizontal:MODE:SCAle?")))))
    #print(scope.query("HORizontal:MODE:SAMPLERate?"))


    #scope.write("HORizontal:Mode AUTO")
    #scope.write("HORizontal:MODE:AUTO:LIMITrecordlen 100000")

    #print(scope.query("HORizontal:Mode?"))
    #scope.write("HORizontal:Mode MANUAL")
    #scope.write("HORizontal:MODE:SAMPLERate {0}".format(1e9))



def findSettledPositionDEPRECATED(wave):
    """ deprecated, not in use because quite complex and not efficient"""
    for i in range(0, len(wave)):
        for j in range(i-200, i+200):
            if(wave[i] <=0 or j >= acqlen or wave[j] <=0):
                passed = False
                break
            if(j < 0 or np.abs((wave[j]-wave[i])/wave[i]) > 0.05):
                passed = False
                break
            else:
            passed = True
            
        if passed:
            break

    return i
