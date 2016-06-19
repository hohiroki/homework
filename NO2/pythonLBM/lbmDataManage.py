import numpy as np

class DataManage(object):
    def __init__(self, filename):
        self.uv = np.load(filename)
        self.u = self.uv[0]
        self.v = self.uv[1]
        self.NX = self.u.shape[0]-1
        self.NY = self.u.shape[1]-1

    def plot_u_contour(self, ff=None):
        """ """
        import matplotlib.pyplot as pl
        pl.figure(figsize=[5,5],dpi=120)
        x = np.linspace(0,1,self.NX+1)
        y = np.linspace(0,1,self.NX+1)
        Y,X=np.meshgrid(x,y)
        CS = pl.contour(X,Y,self.u)
        pl.clabel(CS, inline=1, fontsize=10)
        if ff is None:
            pl.show()
        else:
            pl.show()
            pl.savefig(ff)

    def plot_v_contour(self, ff=None):
        """ """
        import matplotlib.pyplot as pl
        pl.figure(figsize=[5,5],dpi=120)
        x = np.linspace(0,1,self.NX+1)
        y = np.linspace(0,1,self.NX+1)
        Y,X=np.meshgrid(x,y)
        CS = pl.contour(X,Y,(self.u**2+self.v**2)**0.5)
        pl.clabel(CS, inline=1, fontsize=10)
        if ff is None:
            pl.show()
        else:
            pl.show()
            pl.savefig(ff)

    def plot_streamline(self, ff=None):
        """ """
        import matplotlib.pyplot as pl
        pl.figure(figsize=[5,5],dpi=120)
        x = np.linspace(0,1,self.NX+1)
        y = np.linspace(0,1,self.NX+1)
        X,Y=np.meshgrid(x,y)
        pl.streamplot(X,Y,self.u.transpose(),self.v.transpose(),color=self.u*0, linewidth=1)
        pl.axis([0,1,0,1])
        #pl.colorbar()
        if ff is None:
            pl.show()
        else:
            pl.show()
            pl.savefig(ff)

if __name__ == "__main__":
    datamanage = DataManage("re1000.npy")
    datamanage.plot_streamline()