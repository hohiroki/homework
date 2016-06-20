import sys
import numpy as np

class Datamanage(object):
    def __init__(self, filename=None):
        self.rawdata = 'not input yet'
        if filename is not None:
            self.datRead(filename)


    def datRead(self, filename):
        f = open(filename,'rb')
        self.rawdata = f.read()


    def get_u(self):
        pass

    def get_v(self):
        pass

    def get_mach(self):
        pass

    def get_ro(self):
        pass

    def plot_colorf(self):
        pass

    def print_data(self):
        print self.rawdata[100]

    def save_npy(self):
        np.save('savenpy',self.rawdata)

    def read_npy(self, filename):
        self.uv = np.load(filename)
        self.u = self.uv[2]
        print self.u



if __name__ == "__main__":
    manage = Datamanage('1.dat')
    manage.save_npy()
    manage.read_npy('savenpy.npy')