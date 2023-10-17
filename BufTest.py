# -*- coding: utf-8 -*-
"""
DAQ Buffer Measuring Test v1.0

Lukas Kostal, 16.10.2023, ICL
"""


import numpy as np
import matplotlib.pyplot as plt
from importlib import reload
import time

import SerLib as sl
sl = reload(sl)


# specify serial ports for power supply and picoammeter
psu_port = '/dev/tty.usbserial-14230'
pam_port = '/dev/tty.usbserial-14220'

# initialise measurement instruments
psu = sl.psu(psu_port, 9600)
pam = sl.pam(pam_port, 9600)

# ensure serial ports are open
psu.opn()
pam.opn()

# set timeout to 100ms
psu.timeout(100)
pam.timeout(100)

# check psu and pam connected with IDN? query
IDN_psu = psu.IDN()
IDN_pam = pam.IDN().split(',')[0]

# print connected instrument status
print()
print(f'PSU connected: {IDN_psu}')
print(f'PAM connected: {IDN_pam}')

# set psu output to 0V and enable
psu.VSET(0)
psu.OUT(True)

# reset
pam.RST()
# set trigger delay to 0 and count to 2000
pam.TRIG(0, 100)
# set integration rate to 0.01
pam.NPLC(0.01)
# set range to 2nA
pam.RANG(2e-3)
# disable zero check
pam.ZCH(False)
# disable autozero
pam.AZER(False)
# disable display
#pam.DISP(False)
# clear status model
pam.CLS()
# set buffer size to 2000
pam.POIN(100)
# clear buffer
pam.CLE()
# set storage control to start on next reading
pam.FEED()
# enable buffer full measurement event
pam.FULL()
# send operation complete query
opc = pam.OPC()
print(opc)
# start data taking
pam.INIT()
# wait for 10s
#time.sleep(10)
# enable display again
#pam.DISP(True)
# request data from buffer
data = pam.DATA()
print(data)






