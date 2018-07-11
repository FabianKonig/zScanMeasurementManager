import visa
import numpy as np
from struct import unpack
import pylab

rm = visa.ResourceManager()
scope = rm.open_resource('USB0::0x0699::0x0501::C012801::INSTR')
#scope = rm.open_resource('USB0::0x0699::0x03A3::C011056::INSTR')

scope.write('DATA:SOU CH4')
scope.write('DATA:WIDTH 1')
scope.write('DATA:ENC RPB')





ymult = float(scope.query('WFMPRE:YMULT?'))
yzero = float(scope.query('WFMPRE:YZERO?'))
yoff = float(scope.query('WFMPRE:YOFF?'))
xincr = float(scope.query('WFMPRE:XINCR?'))













scope.write('CURVE?')
data = scope.read_raw()

headerlen = 2 + int(data[1])
header = data[:headerlen]
ADC_wave = data[headerlen:-1]

ADC_wave = np.array(unpack('%sB' % len(ADC_wave),ADC_wave))
Volts = (ADC_wave - yoff) * ymult  + yzero
Time = np.arange(0, xincr * len(Volts), xincr)

pylab.plot(Time, Volts)
pylab.show()