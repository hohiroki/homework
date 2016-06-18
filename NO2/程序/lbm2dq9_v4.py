#!/usr/bin/env python
"""
lbm 2dq9
"""
import numpy as np

class LBM2DQ9:
    """ """
    Q = 9
    #e = [[0,0],[1,0],[0,1],[-1,0],[0,-1],[1,1],[-1,1],[-1,-1],[1,-1]]
    w = [4./9, 1./9, 1./9, 1./9, 1./9, 1./36, 1./36, 1./36, 1./36]
    def __init__(self, nx, ny):
        """ """
        self.rho = np.ones((nx+1, ny+1)) 
        self.u = np.zeros((nx+1, ny+1)) 
        self.v = np.zeros((nx+1, ny+1)) 
        self.u0 = np.zeros((nx+1, ny+1)) 
        self.v0 = np.zeros((nx+1, ny+1)) 
        self.f = np.zeros((self.Q, nx+1, ny+1)) 
        self.NX = nx
        self.NY = ny

    def initialize(self, Re, U=0.1, rho0=1.):
        """ """
        dx = dy = 1.
        self.Lx = dx*self.NY
        self.Ly = dy*self.NX
        dt = dx
        c = dx/dt
        nu = U*self.Lx/Re
        self.tau_f = 3*nu + 0.5
        print("tau_f = %.5f"%self.tau_f)
        self.u[:] = 0.
        self.u[:,self.NY] = U
        self.v[:] = 0.
        self.rho[:] = rho0
        u = self.u
        v = self.v
        u2 = u*u
        v2 = v*v
        uv = 1.5*(u2 + v2)
        upv = u + v
        umv = u - v
        #e = [[0,0],[1,0],[0,1],[-1,0],[0,-1],[1,1],[-1,1],[-1,-1],[1,-1]]
        self.f[0] = self.w[0]*self.rho*(1. - uv)
        self.f[1] = self.w[1]*self.rho*(1. + 3.*u + 4.5*u2 - uv)
        self.f[2] = self.w[2]*self.rho*(1. + 3.*v + 4.5*v2 - uv)        
        self.f[3] = self.w[3]*self.rho*(1. - 3.*u + 4.5*u2 - uv)
        self.f[4] = self.w[4]*self.rho*(1. - 3.*v + 4.5*v2 - uv)
        self.f[5] = self.w[5]*self.rho*(1. + 3.*upv + 4.5*upv*upv - uv)
        self.f[6] = self.w[6]*self.rho*(1. - 3.*umv + 4.5*umv*umv - uv)
        self.f[7] = self.w[7]*self.rho*(1. - 3.*upv + 4.5*upv*upv - uv)
        self.f[8] = self.w[8]*self.rho*(1. + 3.*umv + 4.5*umv*umv - uv)

    def _feq2(self, k, rho, u, v):
        """ """
        uv = u*u + v*v
        #e = [[0,0],[1,0],[0,1],[-1,0],[0,-1],[1,1],[-1,1],[-1,-1],[1,-1]]
        if k == 0:
            return self.w[k]*rho*(1. - 1.5*uv)
        elif k == 1:
            return self.w[k]*rho*(1. + 3.*u + 4.5*u*u - 1.5*uv)
        elif k == 2:
            return self.w[k]*rho*(1. + 3.*v + 4.5*v*v - 1.5*uv)
        elif k == 3:
            return self.w[k]*rho*(1. - 3.*u + 4.5*u*u - 1.5*uv)
        elif k == 4:
            return self.w[k]*rho*(1. - 3.*v + 4.5*v*v - 1.5*uv)
        elif k == 5:
            eu = u + v
            return self.w[k]*rho*(1. + 3.*eu + 4.5*eu*eu - 1.5*uv)
        elif k == 6:
            eu = -u + v
            return self.w[k]*rho*(1. + 3.*eu + 4.5*eu*eu - 1.5*uv)
        elif k == 7:
            eu = u + v
            return self.w[k]*rho*(1. - 3.*eu + 4.5*eu*eu - 1.5*uv)
        elif k == 8:
            eu = u - v
            return self.w[k]*rho*(1. + 3.*eu + 4.5*eu*eu - 1.5*uv)

    def set_u_boundary(self):
        """ """
        self.v[0,:] = 0.
        self.v[self.NX,:] = 0.
        self.v[:,0] = 0.
        self.v[:,self.NY] = 0.
        self.u[0,:] = 0.
        self.u[self.NX,:] = 0.
        self.u[:,0] = 0.
        self.u[:,self.NY] = self.U        

    #@profile
    def evolution(self):
        """ """
        nx = self.NX
        ny = self.NY
        u = self.u
        v = self.v
        u2 = u*u
        v2 = v*v
        uv = 1.5*(u2 + v2)
        upv = u + v
        umv = u - v
        qq = 1./self.tau_f
        tt = (1. - qq)        
        #e = [[0,0],[1,0],[0,1],[-1,0],[0,-1],[1,1],[-1,1],[-1,-1],[1,-1]]
        G = (1. - uv)
        F = qq*self.w[0]*self.rho*G + tt*self.f[0]
        self.f[0,1:-1,1:-1] = F[1:-1,1:-1]
        G += 3.*u + 4.5*u2
        F = qq*self.w[1]*self.rho*G + tt*self.f[1]
        self.f[1,1:-1,1:-1] = F[0:-2,1:-1]
        G -= 6.*u
        F = qq*self.w[3]*self.rho*G + tt*self.f[3]
        self.f[3,1:-1,1:-1] = F[2:,1:-1]
        G = (1. + 3.*v + 4.5*v2 - uv)
        F = qq*self.w[2]*self.rho*G + tt*self.f[2]
        self.f[2,1:-1,1:-1] = F[1:-1,0:-2]
        G -= 6.*v
        F = qq*self.w[4]*self.rho*G + tt*self.f[4]
        self.f[4,1:-1,1:-1] = F[1:-1,2:]
        G = (1. + 3.*upv + 4.5*upv*upv - uv)
        F = qq*self.w[5]*self.rho*G + tt*self.f[5]
        self.f[5,1:-1,1:-1] = F[0:-2,0:-2]
        G -= 6.*upv
        F = qq*self.w[7]*self.rho*G + tt*self.f[7]
        self.f[7,1:-1,1:-1] = F[2:,2:]
        G = (1. - 3.*umv + 4.5*umv*umv - uv)
        F = qq*self.w[6]*self.rho*G + tt*self.f[6]
        self.f[6,1:-1,1:-1] = F[2:,0:-2]
        G += 6.*umv
        F = qq*self.w[8]*self.rho*G + tt*self.f[8]
        self.f[8,1:-1,1:-1] = F[0:-2,2:]
        
        self.u0[:] = self.u
        self.v0[:] = self.v
        self.rho = np.sum(self.f, axis=0)
        #self.u = self.f[1]-self.f[3]+self.f[5]-self.f[6]-self.f[7]+self.f[8]
        #self.v = self.f[2]-self.f[4]+self.f[5]+self.f[6]-self.f[7]-self.f[8]
        f57 = self.f[5]-self.f[7]
        f68 = self.f[6]-self.f[8]
        self.u = self.f[1]-self.f[3]+f57-f68
        self.v = self.f[2]-self.f[4]+f57+f68
        self.u /= self.rho
        self.v /= self.rho
        # boundary
        self.set_u_boundary()
        self.rho[nx,1:-1] = self.rho[nx-1,1:-1]
        self.rho[0,1:-1] = self.rho[1,1:-1]
        for k in xrange(self.Q):
            self.f[k,nx,1:-1] = self._feq2(k,self.rho[nx,1:-1],self.u[nx,1:-1], self.v[nx,1:-1]) + self.f[k,nx-1,1:-1] - self._feq2(k,self.rho[nx-1,1:-1],self.u[nx-1,1:-1], self.v[nx-1,1:-1])
            self.f[k,0,1:-1] = self._feq2(k,self.rho[0,1:-1],self.u[0,1:-1], self.v[0,1:-1]) + self.f[k,1,1:-1] - self._feq2(k,self.rho[1,1:-1],self.u[1,1:-1], self.v[1,1:-1])
        self.rho[:,0] = self.rho[:,1]
        self.rho[:,ny] = self.rho[:,ny-1]
        for k in xrange(self.Q):
            self.f[k,:,0] = self._feq2(k,self.rho[:,0],self.u[:,0], self.v[:,0]) + self.f[k,:,1] - self._feq2(k,self.rho[:,1],self.u[:,1], self.v[:,1])
            self.f[k,:,ny] = self._feq2(k,self.rho[:,ny],self.u[:,ny], self.v[:,ny]) + self.f[k,:,ny-1] - self._feq2(k,self.rho[:,ny-1],self.u[:,ny-1], self.v[:,ny-1])
            
    def error(self):
        """ """
        t1 = np.sum((self.u - self.u0)**2 + (self.v - self.v0)**2)
        t2 = np.sum(self.u*self.u0 + self.v*self.v0)
        return (t1/(t2+1e-50))**0.5

    def run(self, Re, U=0.1, rho0=1., step=np.inf, init=False, output=False):
        """ """
        self.U = U
        if init:
            self.initialize(Re, U, rho0)
        n = 0
        while n < step:
            n += 1
            self.evolution()
            if n%100 == 0:
                err = self.error()
                print('%5d step: u[NX/2][NY/2] = {%.6f, %.6f}, err = %g'%(n, self.u[self.NX/2,self.NY/2], self.v[self.NX/2,self.NY/2], err))
                if output is True and n%1000 == 0:
                    np.save('pcavity_%06d.dat'%n, (self.u, self.v))
                if err < 1e-6:
                    print('converged')
                    break
        print('done.')

    def plot_u_contour(self, ff=None):
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
            pl.savefig(ff)

    def plot_streamline(self, ff=None):
        """ """
        import matplotlib.pyplot as pl
        pl.figure(figsize=[5,5],dpi=120)
        x = np.linspace(0,1,self.NX+1)
        y = np.linspace(0,1,self.NX+1)        
        X,Y=np.meshgrid(x,y)
        pl.streamplot(X,Y,self.u.transpose(),self.v.transpose(),color=self.u, linewidth=2, cmap=pl.cm.autumn)
        pl.axis([0,1,0,1])
        pl.colorbar()
        if ff is None:
            pl.show()
        else:
            pl.savefig(ff)

if __name__ == '__main__':
    import sys
    lbm = LBM2DQ9(512,512)
    if len(sys.argv) > 1:
        #print sys.argv[1]
        lbm.run(Re=1000., step=int(sys.argv[1]), init=True)
    else:
        lbm.run(Re=1000., init=True, output=True)
    
