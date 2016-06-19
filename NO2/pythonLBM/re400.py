from lbm2dq9 import *

if __name__ == '__main__':
    import sys
    lbm = LBM2DQ9(256,256)
    lbm.run(Re=400., filename='re400',init=True)