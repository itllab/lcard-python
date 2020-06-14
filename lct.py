# -*- coding: utf-8 -*-
import ctypes
from ctypes import cdll, pointer, c_ulong, POINTER
from ctypes import *
import time
import sys
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, show
import numpy.fft as fft
from numpy import arange, sin, pi
#from numpy import *
#import numpy as np
#from scipy.signal import argrelextrema
#from scipy import fftpack


class slot(ctypes.Structure):
    _fields_ = [('Base', ctypes.c_ulong),
                ('BaseL', ctypes.c_ulong),
                ('Base1', ctypes.c_ulong),
                ('BaseL1', ctypes.c_ulong),
                ('Mem', ctypes.c_ulong),
                ('MemL', ctypes.c_ulong),
                ('Mem1', ctypes.c_ulong),
                ('MemL1', ctypes.c_ulong),
                ('Irq', ctypes.c_ulong),
                ('BoardType', ctypes.c_ulong),
                ('DSPType', ctypes.c_ulong),
                ('Dma', ctypes.c_ulong),
                ('DmaDac', ctypes.c_ulong),
                ('DTA_REG', ctypes.c_ulong),
                ('IDMA_REG', ctypes.c_ulong),
                ('CMD_REG', ctypes.c_ulong),
                ('IRQ_RST', ctypes.c_ulong),
                ('DTA_ARRAY', ctypes.c_ulong),
                ('RDY_REG', ctypes.c_ulong),
                ('CFG_REG', ctypes.c_ulong)]

class read(ctypes.Structure):
    _fields_ = [('SerNum', ctypes.c_char*9),
                ('BrdName', ctypes.c_char*7),
                ('Rev', ctypes.c_char),
                ('DspType', ctypes.c_char*5),
                ('IsDacPresent', ctypes.c_char),
                ('Quartz', ctypes.c_ulong),
                ('Reserv2', ctypes.c_char*13),
                ('KoefADC', ctypes.c_ushort*8),
                ('KoefDAC', ctypes.c_ushort*4),
                ('Custom', ctypes.c_ushort*32)]

class wadc_par_0(ctypes.Structure):
    _pack_ = 1
    _fields_ = [
                ('s_Type', ctypes.c_ulong),
                ('FIFO', ctypes.c_ulong),
                ('IrqStep', ctypes.c_ulong),
                ('Pages', ctypes.c_ulong),
                
                ('AutoInit', ctypes.c_ulong),
                
                ('dRate', ctypes.c_double),
                ('dKadr', ctypes.c_double),
                ('dScale', ctypes.c_double),
                ('Rate', ctypes.c_ulong),
                ('Kadr', ctypes.c_ulong),
                ('Scale', ctypes.c_ulong),
                ('FPDelay', ctypes.c_ulong),
                
                ('SynchroType', ctypes.c_ulong),
                ('SynchroSensitivity', ctypes.c_ulong),
                ('SynchroMode', ctypes.c_ulong),
                ('AdChannel', ctypes.c_ulong),
                ('AdPorog', ctypes.c_ulong),
                ('NCh', ctypes.c_ulong),
                ('Chn', ctypes.c_ulong * 16),
                ('IrqEna', ctypes.c_ulong),
                ('AdcEna', ctypes.c_ulong)]

wl = cdll.wlcomp
hDll = ctypes.pointer(c_ulong(wl.LoadAPIDLL('lcomp.dll')))
hErr = ctypes.pointer(c_ulong())
hIfc = ctypes.pointer(c_ulong(wl.CallCreateInstance(hDll, 0, hErr)))
print 'hDll', hDll.contents.value
print 'hIfc', hIfc.contents.value
print 'hErr', hErr.contents.value

Open = ctypes.pointer(c_ulong(wl.OpenLDevice(hIfc)))
print 'Open', Open.contents.value

Bios = ctypes.pointer(c_ulong(wl.LoadBios(hIfc, 'E440')))
print 'Bios', Bios.contents.value

Test = ctypes.pointer(c_ulong(wl.PlataTest(hIfc)))
print 'Test', Test.contents.value

sl = ctypes.pointer(slot())
wl.GetSlotParam(hIfc, sl)
print 'Slot', sl.contents.BoardType

pd = ctypes.pointer(read())
print 'ReadPlataDescr', wl.ReadPlataDescr(hIfc, pd)
print 'SerNum', pd.contents.SerNum
print 'BrdName', pd.contents.BrdName
print 'Rev', pd.contents.Rev
print 'DspType', pd.contents.DspType
print 'IsDacPresent', ord(pd.contents.IsDacPresent)
print 'Quartz', pd.contents.Quartz, 'Hz'

L_ADC_PARAM = 1
L_STREAM_ADC = 1
sp = 2
L_POINT_SIZE = 10001

pp = ctypes.pointer(wadc_par_0())

pp.contents.s_Type = L_ADC_PARAM
pp.contents.AutoInit = 1
pp.contents.dRate = 400.0
pp.contents.dKadr = 0.0
pp.contents.dScale = 0.0
pp.contents.SynchroType = 0
pp.contents.SynchroSensitivity = 0
pp.contents.SynchroMode = 0
pp.contents.AdChannel = 0
pp.contents.AdPorog = 0
pp.contents.NCh = 16
pp.contents.Chn[0] = 0
pp.contents.Chn[1] = 1
pp.contents.Chn[2] = 2
pp.contents.Chn[3] = 3
pp.contents.Chn[4] = 4
pp.contents.Chn[5] = 5
pp.contents.Chn[6] = 6
pp.contents.Chn[7] = 7
pp.contents.Chn[8] = 8
pp.contents.Chn[9] = 9
pp.contents.Chn[10] = 10
pp.contents.Chn[11] = 11
pp.contents.Chn[12] = 12
pp.contents.Chn[13] = 13
pp.contents.Chn[14] = 14
pp.contents.Chn[15] = 15
pp.contents.FIFO = 4096
pp.contents.IrqStep = 4096
pp.contents.Pages = 32
pp.contents.IrqEna = 1
pp.contents.AdcEna = 1

Size = ctypes.pointer(ctypes.c_ulong(1000000))
Data = ctypes.pointer(ctypes.c_ushort())
Sync = ctypes.pointer(ctypes.c_ulong())

poi = ctypes.pointer(c_ulong())

print 'RequestBufferStream', wl.RequestBufferStream(hIfc, Size, L_STREAM_ADC)
print 'Allocated memory size(word): ', Size[0]

print 'FillDAQparameters', wl.FillDAQparameters(hIfc, pp, sp)

print '.......... Buffer size(word):      ', Size[0]
print '.......... Pages:                  ', pp.contents.Pages
print '.......... IrqStep:                ', pp.contents.IrqStep
print '.......... FIFO:                   ', pp.contents.FIFO
print '.......... Rate:                   ', pp.contents.dRate
print '.......... Kadr:                   ', pp.contents.dKadr

print 'SetParametersStream', wl.SetParametersStream(hIfc, pp, sp, Size, 
                                                    ctypes.cast(ctypes.pointer(Data), ctypes.POINTER(c_void_p)), 
                                                    ctypes.cast(ctypes.pointer(Sync), ctypes.POINTER(c_void_p)), L_STREAM_ADC)

print '.......... Used buffer size(points):', Size[0]
print '.......... Pages:                   ', pp.contents.Pages
print '.......... IrqStep:                 ', pp.contents.IrqStep
print '.......... FIFO:                    ', pp.contents.FIFO
print '.......... Rate:                    ', pp.contents.dRate
print '.......... Kadr:                    ', pp.contents.dKadr

#print 'GetParameter', wl.GetParameter(hIfc, L_POINT_SIZE, poi)
#print '.......... Point size:', poi[0]

print 'EnableCorrection', wl.EnableCorrection(hIfc, 1)
print 'InitStartDevice', wl.InitStartLDevice(hIfc)
print 'StartDevice', wl.StartLDevice(hIfc)

fr = pp.contents.dRate * 1000   #Hz
tr = 1 / fr                     #s
x1 = []
y = []
x4 = []
y4 = []
#-----------------------ver1-------------------------------------------
'''
i = 0
while i < 65535:
    time.sleep(0.0001)
    if Data[i] < 10000:
        x2 = Data[i]*(10.0/8000)
        x1.append(x2)
        y.append(i*tr)
    elif Data[i] > 10000:
        x3 = (Data[i]-65536)*(10.0/8000)
        x1.append(x3)
        y.append(i*tr)
    i = i + 1
su = sum(x1[500:]) / 65536
sur = round(su, 4)
print sur
#print ((max(x1[500:]) - min(x1[500:])) / (1.4142135 * 2))
'''
#---------------------ver2-----------------------------------------------
'''N = pp.contents.NCh
k = 0
while k < N:
    i = 0 + k
    time.sleep(0.1)
    while i < 131072:       
        if Data[i] < 10000:
            x2 = Data[i]*(10.0/8000)
            x1.append(x2)
            y.append(i*tr)
        elif Data[i] > 10000:
            x3 = (Data[i]-65536)*(10.0/8000)
            x1.append(x3)
            y.append(i*tr)
        i = i + N
    dc = sum(x1) / 131072
    dc = round(dc, 4)
    print 'Chn', k, dc
    #ac = (max(x1) - min(x1)) / (1.4142135 * 2)
    #ac = round(ac, 4)
    #print 'Chn', k, ac
    x1 = []
    k = k + 1'''
#---------------------ver3------------------------------------------------
N = pp.contents.NCh
k = 0
while k < N:
    i = 0 + k
    time.sleep(0.4)
    while i < 131072:       
        if Data[i] < 10000:
            x2 = Data[i]*(10.0/8000)
        elif Data[i] > 10000:
            x2 = (Data[i]-65536)*(10.0/8000)
        x1.append(x2)
        y.append(i*tr)
        i = i + N
    x4.append(x1)
    y4.append(y)
    x1 = []
    y = []
    k = k + 1

k = 0
while k < N:
    dc = sum(x4[k])/(131072/N)                          # постоянное напряжение
    dc = round(dc, 4)
    print 'Chn', k, dc
    '''ac = (max(x4[k]) - min(x4[k])) / (1.4142135 * 2) # переменное напряжение
    ac = round(ac, 4)
    print 'Chn', k, ac'''
    k = k + 1
#------------------------------------------------------------------------
print 'StopDevice', wl.StopLDevice(hIfc)
print 'CloseDevice', wl.CloseLDevice(hIfc)

x5 = x4[2]
spectrum = fft.fft(x5)
freq = fft.fftfreq(len(spectrum))

threshold = 0.5 * max(abs(spectrum))
mask = abs(spectrum) > threshold
peaks = freq[mask] * fr / N
print peaks

fig = figure(1)

ax1 = fig.add_subplot(211)
ax1.plot(y4[2], x4[2])
ax1.grid(True)
ax1.set_xlim((0.002, 0.009))
ax1.set_ylabel('U, V')
l1=ax1.set_title('t, sec')
l1.set_color('g')
l1.set_fontsize('large')

ax2 = fig.add_subplot(212)
ax2.plot(freq*fr/N, abs(spectrum))
ax2.grid(True)
ax2.set_xlim((-5000, 5000))
l2 = ax2.set_xlabel('f, Hz')
l2.set_color('g')
l2.set_fontsize('large')

show()
#--------------------------------------------------