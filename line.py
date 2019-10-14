
# Class of system lines
import numpy as np


class Line:

    def __init__(self, dataline):

        self.origin = int(dataline[0:4].strip())
        self.destiny = int(dataline[4:12].strip())

        if dataline[16:23].strip():
            self.R = float(dataline[16:23])/100
        else:
            self.R = 0.

        if dataline[23:29].strip():
            self.X = float(dataline[23:29])/100
        else:
            self.X = 0.

        if dataline[29:35].strip():
            self.B = float(dataline[29:35])/100
        else:
            self.B = 0.

        # Power flowing from origin to destiny and destiny to origin
        self.S_od = 0
        self.S_do = 0

    def save_flow(self, p, q, bus_origin):
        if bus_origin == self.origin:
            self.S_od = p + 1j*q
        elif bus_origin == self.destiny:
            self.S_do = p + 1j * q
        else:
            print("ERROR: The bus must be at one of the ends of the line")
