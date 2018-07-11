import visa 
import numpy as np
from struct import unpack
import matplotlib.pyplot as plt
import traceback
import warnings
import time



# Time of oscillations of the normal ThorLabs photodiodes.
CONSTANT_TIMEOFOSCILLATIONS = 1e-6


class TektronixCommunicator:

    def __init__(self):
        # Channels to be queried:
        self.channels = ["CH1", "CH2", "CH4"] # needs to be ordered
        self.channelSettings = np.empty(shape=(len(self.channels), 3)) #for each channel: ymult, yoff, yzero
        self.data = None
        self.xincr = None
        self.acqlen = None
        self.waveforms = None
        self.scope = None


    def allocate_scope(self):
        rm = visa.ResourceManager()
        with warnings.catch_warnings(): # Ignore the visa io warning
            warnings.simplefilter("ignore")
            self.scope = rm.open_resource('USB0::0x0699::0x0501::C012801::INSTR')

    def deallocate_scope(self):
        self.scope.close()
        self.scope=None


    def acquireCurves(self):
        #rm = visa.ResourceManager()
        #with warnings.catch_warnings(): # Ignore the visa io warning
        #    warnings.simplefilter("ignore")
        #    scope = rm.open_resource('USB0::0x0699::0x0501::C012801::INSTR')
        print("write 1")
        self.scope.write('DATA:ENC RPB') # encoding
        self.scope.write('DATA:WIDTH 1')
        self.scope.write('Data:Stop 1e20') # read the entire time axis.

        # get the y-axis settings of each channel
        i = 0
        for channel in self.channels:
            print("write 2")
            self.scope.write("SELect:" + channel + " ON") #activate channel on oscilloscope, otherwise error
            self.scope.write("DATA:SOU " + channel)
            self.channelSettings[i,0] = float(self.scope.query('WFMPRE:YMULT?'))
            self.channelSettings[i,1] = float(self.scope.query('WFMPRE:YOFF?'))
            self.channelSettings[i,2] = float(self.scope.query('WFMPRE:YZERO?'))
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
        
        print("write 3")
        self.scope.write("DATA:SOU " + channelsString)

        # get the x-axis settings
        print("x axis settings query")
        self.xincr = float(self.scope.query('WFMPRE:XINCR?'))
        self.acqlen = int(self.scope.query("HORizontal:ACQLENGTH?")) # length of a dataset

        # Acquire the curves
        print("wirte curvenext")
        self.scope.write('CURVENext?')
        time.sleep(0.2)
        print("read raw")
        data = self.scope.read_raw()
        time.sleep(0.2)
        self.data = np.array(unpack('%sB' % len(data), data))



    def processCurves(self):
        headerlength = (len(self.data)-1) % (len(self.channels)*self.acqlen) // len(self.channels)
        time = np.linspace(0, self.xincr * (self.acqlen-1), self.acqlen)

        self.waveforms = np.empty(shape=(1+len(self.channels), self.acqlen))
        self.waveforms[0] = time

        for i in range(len(self.channels)):
            start = (1+i) * headerlength + i*self.acqlen
            stop = (1+i)*(headerlength+self.acqlen)

            wave = self.data[start:stop]
            self.waveforms[i+1] = (wave - self.channelSettings[i,1]) * self.channelSettings[i,0]  + self.channelSettings[i,2]


    def getSignalAfterPhotodiodeOscillations(self, waveform):
        """ Find the position of the photodiode signal when the oscillations have stopped. """
        
        maxPos = self.getPositionOfMaximum(waveform)
        requiredSteps = CONSTANT_TIMEOFOSCILLATIONS // self.xincr
        pos = int(maxPos+requiredSteps)

        if pos < self.acqlen:
            return waveform[int(maxPos+requiredSteps)]
        else:
            return waveform[-1]


    def getOsciValues(self):

        deallocate_necessary = False
        
        if self.scope is None:
            deallocate_necessary = True
            self.allocate_scope()

        while(True):
            try:
                print("Status byte: ", self.scope.query('*STB?'))
                self.acquireCurves()
                self.processCurves()
                break
            except Exception as ex:
                print("\nTektronixCommunicator encountered an exception, will try again!!!!")
                print(ex)
                print("\n")
                time.sleep(1)
                self.scope.write('CLEAR ALL')
                self.scope.write('*CLS')
                #traceback.print_exc() #Doesn't work with Visa IO error

        if deallocate_necessary:
            self.deallocate_scope()

        # correct channel 1 and 2 for their offsets:
        self.waveforms[1] = self.correctForOffset(self.waveforms[1])
        self.waveforms[2] = self.correctForOffset(self.waveforms[2])

        max_signals = np.empty(shape=(len(self.channels)), dtype=float)
        max_signals[0] = self.computeLocalAverage(self.waveforms[1], pos="max")  # Channel 1
        max_signals[1] = self.computeLocalAverage(self.waveforms[2], pos="max")  # Channel 2
        max_signals[2] = self.getSignalAfterPhotodiodeOscillations(self.waveforms[3]) # Channel 4

        return max_signals


    def getSignalPeaks(self, iterations):
        """ Output: 
            2-dim numpy array with the first dimension denoting the channel and the second
            dimension is iterations entries long, each storing the maximum value of that
            measurement. The order of the channels is 1: ref, 2: oa, 3: ca.
        """
        max_signal_values = np.zeros(shape=(len(self.channels), iterations))

        self.allocate_scope()

        for iteration_index in range(iterations):
            osciVals = self.getOsciValues()
            # iterate through self.channels from right to left [CA, OA, REF], such that the result is [REF, OA, CA].
            for channel_index in range(-1,-len(self.channels)-1,-1):
                max_signal_values[np.abs(channel_index)-1, iteration_index] = osciVals[channel_index]

        self.deallocate_scope()

        return max_signal_values


    def computeLocalAverage(self, wave, pos):
        """ compute the local average of a waveform around the given position. For this, a certain
            amount of data points after position are taken into consideration. If position=="max",
            the local average around the maximum position of wave is returned.
        """
        AMOUNT_NEXT_DATA_POINTS = 40

        assert isinstance(pos, int) or pos == "max"

        if pos == "max":
            pos = self.getPositionOfMaximum(wave)
        return wave[pos : pos+AMOUNT_NEXT_DATA_POINTS].mean()


    def getPositionOfMaximum(self, waveform):
        maxPos, = np.where(waveform==waveform.max())
        return int(maxPos[0])


    def correctForOffset(self, wave):
        """ Correct for the offset which is assumed at the beginning of the trace."""
        return wave - self.computeLocalAverage(wave, pos=0)



if __name__ == '__main__':

    # plot waveforms
    colors = ["orange", "blue", "green", "red"]

    scopeCommunicator = TektronixCommunicator()
    sc = scopeCommunicator
    data = sc.getSignalPeaks(5)

    #print("{0:.5f}    {1:.5f}".format(data[0].mean(), data[0].std(ddof=1)))

    
    max_signals = sc.getOsciValues()

    plt.axhline(max_signals[0], marker="o", color=colors[0])
    plt.axhline(max_signals[1], marker="o", color=colors[1])
    plt.axhline(max_signals[2], marker="o", color=colors[2])

    for i in range(len(sc.channels)):
        plt.plot(sc.waveforms[0], sc.waveforms[i+1], color=colors[i], zorder=5-i)

    plt.grid()
    plt.show()

    data = sc.getSignalPeaks(5)
    data[1] = data[1] / data[0]
    data[2] = data[2] / data[0]
    print("Ref: ", data[0].mean(), data[0].std(ddof=1)/data[0].mean())
    print("OA:  ", data[1].mean(), data[1].std(ddof=1)/data[1].mean())
    print("CA:  ", data[2].mean(), data[2].std(ddof=1)/data[2].mean())





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
