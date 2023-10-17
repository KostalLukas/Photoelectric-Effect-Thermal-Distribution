# -*- coding: utf-8 -*-
"""
Serial Communication Test v1.0

Lukas Kostal, 12.10.2023, ICL
"""


import serial
import time


# initialise serial communication with power supply
psu = serial.Serial(
    port='/dev/tty.usbserial-14230',
    baudrate=9600,
    parity=serial.PARITY_NONE,
    stopbits=serial.STOPBITS_ONE,
    bytesize=serial.EIGHTBITS
)

# initialise serial communication with picoammeter 
pam = serial.Serial(
    port='/dev/tty.usbserial-14220',
    baudrate=9600,
    parity=serial.PARITY_NONE,
    stopbits=serial.STOPBITS_ONE,
    bytesize=serial.EIGHTBITS
)

# chech serial ports are open if not open them
if psu.is_open == False:
    psu.open()
if pam.is_open == False:
    pam.open()

# set the serial read timeout to 100ms
psu.timeout = 0.1
pam.timeout = 0.1
    
# string identification query message to be sent
msg = '*IDN?'

# send the query to the power supply, read the response and print it
psu.write(bytes(msg, encoding='utf-8'))
response = psu.readline()
print(response)

# send the query to the picoammeter, read the response and print it
pam.write(bytes(msg+'\r', encoding='utf-8'))
response = pam.readline()
print(response)

# close the serial ports
psu.close()
pam.close()