from lbm2dq9 import *

if __name__ == '__main__':
    import sys
    lbm = LBM2DQ9(256,256)
    lbm.run(Re=2000., filename='re2000',init=True)