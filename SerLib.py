# -*- coding: utf-8 -*-
"""
Serial Interface Library v3.0

Lukas Kostal, 20.10.2023, ICL
"""


import numpy as np
import serial
import time


# class for serial device
class dev:
    """
    Device class to represent an instrument connected over serial

    Attributes
    ----------
    port : string
        Path to the tty serial device
        
    baudrate : integer
        Baudrate in bits per second
        
    serial : object
        Instance of the PySerial serial class
        
    stat : string
        Status of the device serial port 'open' or 'closed'

    timeout : float
        Timeout in ms for recieving data

    Methods
    -------
    __init__(self, port, baudrate=9600)
        Make an instance of the device class with 8 bits, 1 stop bit, no parity

    __repr__(self)
        Representation method for the device class

    __str__(self)
        String method for the device class

    status(self)
        Return status of the device serial port
    
    opn(self)
        Check if device serial port is open if closed open it
        
    close(self)
        Check if device serial port is closed if open close it
        
    timeout(self, t),
        Set timeout for recieving data from the serial device in ms
        
    IDN(self, t),
        Send identity query to the device and return response
    """

    # initialise serial device at given port
    def __init__(self, port, baudrate=9600):
        self.port = port
        self.baudrate = baudrate
        self.serial = serial.Serial(
            port = self.port, 
            baudrate = self.baudrate,
            parity=serial.PARITY_NONE,
            stopbits=serial.STOPBITS_ONE,
            bytesize=serial.EIGHTBITS
        )
        
    # get the status of the serial channel
    def status(self):
        if self.serial.is_open() == True:
            self.stat = 'open'
        else:
            self.stat = 'closed'
            
        return self.stat
        
    # if serial channel is closed open it
    def opn(self):
        if self.serial.is_open == False:
            self.serial.open()
            self.stat = 'open'
    
    # if serial channel is open close it
    def close(self):
        if self.serial.is_open == True:
            self.serial.close()
            self.stat = 'closed'
        
    # set the timeout for recieving data from device in ms
    def timeout(self, t=100):
        self.timeout = t * 1e-3
        self.serial.timeout = t * 1e-3
    
    # send identity query to serial device
    def IDN(self):
        resp = self.Rx('*IDN?')
        
        return resp
        
    
# class for Keithley picoammeter instrument
class pam(dev):
    
    # transmit data over serial
    def Tx(self, string, t=100):
        string += '\r'
        string = bytes(string, encoding='utf-8')
        self.serial.write(string)
        time.sleep(t * 1e-3)
    
    # transmit data over serial and recieve a response
    def Rx(self, string):
        string += '\r'
        string = bytes(string, encoding='utf-8')
        self.serial.write(string)
        
        time.sleep(0.1)
        resp = self.serial.readline()
        resp = str(resp, encoding='utf-8')
        if resp[-1:] == '\r':
            resp = resp[:-2]
        
        return resp
    
    # operation completion query
    def OPC(self):
        resp = self.Rx('*OPC?')
        
        stat = False
        
        if resp == '1':
            stat = True
            
        return stat
    
    # reset the device
    def RST(self):
        self.Tx('*RST')
        
    # trigger a reading
    def INIT(self):
        self.Tx('INIT')
        
    # set range and for autorange specify lower and upper limits on range
    def RANG(self, rang='AUTO', Imin=0, Imax=0):
        self.range = rang
        
        if rang == 'AUTO':
            self.Imin = Imin
            self.Imax = Imax
            
            self.Tx('RANG:AUTO ON')
            
            if Imin != 0 and Imin != 0:
                self.Tx(f'RANG:AUTO:LLIM {Imin:1g}')
                self.Tx(f'RANG:AUTO:ULIM {Imax:1g}')
                
        else:
            self.Tx(f'RANG {rang:.1g}')
    
    # enable zero check
    def ZCH(self, enb):
        if enb == True:
            self.Tx('SYST:ZCH ON')
        if enb == False:
            self.Tx('SYST:ZCH OFF')
        
    # enable zero correct
    def ZCOR(self, enb):
        if enb == True:
            self.Tx('SYST:ZCOR ON')
        if enb == False:
            self.Tx('SYST:ZCOR OFF')
            
    # enable autozero
    def AZER(self, enb):
        if enb == True:
            self.Tx('SYST:AZER ON')
        if enb == False:
            self.Tx('SYST:AZER OFF')
            
    # perform zero offset callibration
    def ZCAL(self):
        self.ZCH(True)
        self.Tx('RANG 2e-9')
        self.Tx('INIT')
        self.Tx('SYST:ZCOR:ACQ')
        self.ZCOR(True)
        self.ZCH(False)
    
    # disable all filters
    def FOFF(self):
        self.Tx('MED OFF')
        self.Tx('AVER OFF')
        
    # trigger a reading and return value as float
    def READ(self):
        meas = self.Rx('READ?')
        
        vals = meas.split(',')
        vals[0] = vals[0][:-1]
        vals = np.array(vals, dtype=float)

        return vals[0]
    
    # set trigger number of trigger events and trigger delay
    def TRIG(self, Ntrig, Tdly=0):
        self.Tx(f'TRIG:DEL {Tdly}')
        self.Tx(f'TRIG:COUN {int(Ntrig)}')
        
    # set integration rate in number of power line cycles PLC
    def NPLC(self, Nplc):
        self.Tx(f'NPLC {Nplc}')
        
    # clear status model
    def CLS(self):
        self.Tx('*CLS')
    
    # clear internal buffer
    def CLE(self):
        self.Tx('TRAC:CLE')
        
    # set internal buffer size
    def POIN(self, Bsize):
        self.Tx(f'TRAC:POIN {int(Bsize)}')
    
    # set buffer storage control to start on next reading
    def FEED(self):
        self.Tx('TRAC:FEED:CONT NEXT')
        
    # enable buffer full measurement event
    def FULL(self):
        self.Tx('STAT:MEAS:ENAB 512')
    
    # enable SQR on buffer full measurement event
    def SRE(self):
        self.Tx('*SRE 1')
        
    # request data from buffer
    def DATA(self):
        data = self.Rx('TRAC:DATA?')
        return data
    
    # enable display
    def DISP(self, stat):
        if stat == True:
            self.Tx('DISP:ENAB ON')
        if stat == False:
            self.Tx('DISP:ENAB OFF')
            
    # function to prepare for reading into buffer
    def BPREP(self, Nmes, Tdly=0, Nplc=0.01):
        self.Nmes = Nmes
        self.Tdly = Tdly
        self.Nplc = Nplc
        
        self.TRIG(Nmes, Tdly)
        self.NPLC(Nplc)
        self.CLS()
        self.POIN(Nmes)
        
    # function to prepare for another measurement to buffer
    def BARM(self):
        self.CLE()
        self.FEED()
        self.FULL()
        
    # function to take measurement into buffer and read buffer
    def BREAD(self):
        self.INIT()
        Twait = (self.Nplc/50 + self.Tdly + 0.01) * self.Nmes + 0.01
        time.sleep(Twait)
        
        data = self.DATA()
        
        vals = data.split(',')
        vals = np.array(vals, dtype=str)
        vals = vals.reshape(-1, 3)
        for i in range(0, len(vals)):
            vals[i, 0] = vals[i, 0][:-1]
        vals = np.array(vals, dtype=float)
        
        return vals[:, 0]

            
    # function to setup and take measurements store them to buffer and read them out
    # Nmes is no of measurements, Tdly is delay between measurements, Nplc is integration time
    def BALL(self, Nmes, Tdly=0, Nplc=0.01):
        self.Nmes = Nmes
        self.Tdly = Tdly
        self.Nplc = Nplc
        
        self.TRIG(Nmes, Tdly)
        self.NPLC(Nplc)
        self.CLS()
        self.POIN(Nmes)
        self.CLE()
        self.FEED()
        self.FULL()
        self.INIT()
        
        Twait = (Nplc/50 + Tdly + 0.01) * Nmes + 0.01
        time.wait(Twait)
        
        data = self.DATA()
        
        vals = data.split(',')
        vals = np.array(vals, dtype=str)
        vals = vals.reshape(-1, 3)
        for i in range(0, len(vals)):
            vals[i, 0] = vals[i, 0][:-1]
        vals = np.array(vals, dtype=float)
        
        return vals[:, 0]
        
    
# class for TENMA power supply instrument
class psu(dev):
        
    # transmit data over serial
    def Tx(self, string):
        string = bytes(string, encoding='utf-8')
        self.serial.write(string)
        time.sleep(0.1)
        
    # transmit data over serial and recieve a response
    def Rx(self, string):
        string = bytes(string, encoding='utf-8')
        self.serial.write(string)
        
        resp = self.serial.readline()
        resp = str(resp, encoding='utf-8')
        
        return resp
    
    # set polarity
    def POL(self, V):
        if hasattr(self, 'pos') == False:
            input('Change polarity to positive')
            print()
            self.pos = True
        
        if V >= 0 and self.pos == False:
            input('Change polarity to positive')
            print()
            self.pos = True
            
        if V < 0 and self.pos == True:
            input('Change polarity to negative')
            print()
            self.pos = False
    
    # set voltage offset
    def OFST(self, V):
        
        if hasattr(self, 'ofst') == False:
            input('Set offset to 0V')
            print()
            self.ofst = 0
        
        if V > 30 and self.ofst < 30:
            input('Set offset voltage to 30V')
            print()
            self.ofst = 30
            
        if V < -30 and self.ofst > -30:
            input('Set offset voltage to -30V')
            print()
            self.ofst = -30
            
        if -30 < V < 30 and self.ofst != 0:
            input('Set offset voltage to 0V')
            print()
            self.ofst = 0
            
        return self
        
    # set the output voltage
    def VSET(self, V):

        self.OFST(V)
        V -= self.ofst
        self.POL(V)
        V_abs = np.abs(V)
        self.Tx(f'VSET1:{V_abs:.2f}')
    
    # read the set output voltage and return as float
    def getVSET(self):
        resp = self.Rx('VOUT?')
        self.V = float(resp)
        
        return self.V
    
    # enable output
    def OUT(self, enb):
        self.out = enb
        
        if enb == True:
            self.Tx('OUT1')
            
        if enb == False:
            self.Tx('OUT0')
            
    # wait at a given voltage
    def WAIT(self, V):
        if self.V == V:            
            input(f'Waiting at V={V:.2f} V')
            print()
            
    # enable beep
    def BEEP(self, enb):
        self.enb = enb
        
        if enb == True:
            self.Tx('BEEP1')
            
        if enb == False:
            self.Tx('BEEP0')