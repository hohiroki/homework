from lbm2dq9 import *

if __name__ == '__main__':
    import sys
    lbm = LBM2DQ9(256,256)
    lbm.run(Re=5000., filename='re4000',init=True, output=True)