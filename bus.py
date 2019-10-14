# Class of system buses
import numpy as np


class Bus:

    def __init__(self, databus):

        bustypes = {'0': 'PQ', '1': 'PV', '2': 'Vθ'}

        self.ID = int(databus[0:4])
        self.name = databus[8:22]
        self.bustype = bustypes[databus[4:8].strip()]

        self.V = float(databus[22:26])/1000

        if databus[26:30].strip():
            self.theta = float(databus[26:30])
        else:
            self.theta = 0.

        p = []
        q = []
        for item in [databus[30:35], databus[56:60]]:
            if not item.strip():
                p.append(0.)
            else:
                p.append(float(item))

        for item in [databus[35:40], databus[60:65]]:
            if not item.strip():
                q.append(0.)
            else:
                q.append(float(item))

        self.P = (p[0] - p[1])/Sb
        self.Q = (q[0] - q[1])/Sb

    # Function to refresh values of voltage and angle on buses
    def refresh(self, ang, v):

        if self.bustype == 'PQ':
            self.theta += ang
            self.V += v
        elif self.bustype == 'PV':
            self.theta += ang

    # Function to save values of power
    def save_power(self, p, q):

        if self.bustype == 'Vθ':
            self.P = p
            self.Q = q
        elif self.bustype == 'PV':
            self.Q = q
