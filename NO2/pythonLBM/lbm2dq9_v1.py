#!/usr/bin/env python
"""
lbm 2dq9
"""
import numpy as np

class LBM2DQ9:
    """ """
    Q = 9
    #e = [[0,0],[1,0],[0,1],[-1,0],[0,-1],[1,1],[-1,1],[-1,-1],[1,-1]]
    ex = np.array([0,1,0,-1,0,1,-1,-1,1])
    ey = np.array([0,0,1,0,-1,1,1,-1,-1])
    w = [4./9, 1./9, 1./9, 1./9, 1./9, 1./36, 1./36, 1./36, 1./36]
    def __init__(self, nx, ny):
        """ """
        nx1 = nx+1
        ny1 = ny+1
        self.rho = np.ones((nx1, ny1)) 
        self.u = np.zeros((nx1, ny1)) 
        self.v = np.zeros((nx1, ny1)) 
        self.u0 = np.zeros((nx1, ny1)) 
        self.v0 = np.zeros((nx1, ny1)) 
        self.f = np.zeros((self.Q, nx1, ny1)) 
        self.F = np.zeros((self.Q, nx1, ny1))
        self.feq = np.zeros((self.Q, nx1, ny1))
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
        self._calc_feq()
        self.f[:] = self.feq[:]      

    def _calc_feq(self):
        """ """
        uv = self.u*self.u + self.v*self.v
        for k in xrange(self.Q):
            eu = self.ex[k]*self.u + self.ey[k]*self.v
            self.feq[k] = self.w[k]*self.rho*(1. + 3.*eu + 4.5*eu*eu - 1.5*uv)

    def _calc_feq2(self, k, rho, u, v):
        """ """
        uv = u*u + v*v
        eu = self.ex[k]*u + self.ey[k]*v
        return self.w[k]*rho*(1. + 3.*eu + 4.5*eu*eu - 1.5*uv)

    def evolution(self):
        """ """
        nx = self.NX
        ny = self.NY
        self._calc_feq()
        pp = 1./self.tau_f
        tt = (1. - pp)
        for k in xrange(self.Q):
            i0 = 1 - self.ex[k]
            i1 = i0 + nx - 1
            j0 = 1 - self.ey[k]
            j1 = j0 + ny - 1
            self.F[k,1:-1,1:-1] = tt*self.f[k,i0:i1,j0:j1] + pp*self.feq[k,i0:i1,j0:j1]
        self.u0[:] = self.u
        self.v0[:] = self.v
        self.u[1:-1,1:-1] = 0.
        self.v[1:-1,1:-1] = 0.
        self.f[:] = self.F[:]
        self.rho[1:-1,1:-1] = np.sum(self.f[:,1:-1,1:-1], axis=0)
        for k in xrange(self.Q):
            self.u += self.ex[k]*self.f[k]
            self.v += self.ey[k]*self.f[k]
        self.u[1:-1,1:-1] /= self.rho[1:-1,1:-1]
        self.v[1:-1,1:-1] /= self.rho[1:-1,1:-1]
        # boundary
        self.rho[nx,1:-1] = self.rho[nx-1,1:-1]
        self.rho[0,1:-1] = self.rho[1,1:-1]
        for k in xrange(self.Q):
            self.f[k,nx,1:-1] = self._calc_feq2(k,self.rho[nx,1:-1],self.u[nx,1:-1], self.v[nx,1:-1]) + self.f[k,nx-1,1:-1] - self._calc_feq2(k,self.rho[nx-1,1:-1],self.u[nx-1,1:-1], self.v[nx-1,1:-1])
            self.f[k,0,1:-1] = self._calc_feq2(k,self.rho[0,1:-1],self.u[0,1:-1], self.v[0,1:-1]) + self.f[k,1,1:-1] - self._calc_feq2(k,self.rho[1,1:-1],self.u[1,1:-1], self.v[1,1:-1])
        self.rho[:,0] = self.rho[:,1]
        self.rho[:,ny] = self.rho[:,ny-1]
        for k in xrange(self.Q):
            self.f[k,:,0] = self._calc_feq2(k,self.rho[:,0],self.u[:,0], self.v[:,0]) + self.f[k,:,1] - self._calc_feq2(k,self.rho[:,1],self.u[:,1], self.v[:,1])
            self.f[k,:,ny] = self._calc_feq2(k,self.rho[:,ny],self.u[:,ny], self.v[:,ny]) + self.f[k,:,ny-1] - self._calc_feq2(k,self.rho[:,ny-1],self.u[:,ny-1], self.v[:,ny-1])
            
    def output(self, m):
        """ """
        np.save('pcavity_%06d.dat'%m, (self.u, self.v))

    def error(self):
        """ """
        t1 = np.sum((self.u - self.u0)**2 + (self.v - self.v0)**2)
        t2 = np.sum(self.u*self.u0 + self.v*self.v0)
        return (t1/(t2+1e-30))**0.5

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
                    self.output(n)
                if err < 1e-6:
                    break
        print('done.')

if __name__ == '__main__':
    import sys
    lbm = LBM2DQ9(256,256)
    if len(sys.argv) > 1:
        #print sys.argv[1]
        lbm.run(Re=1000., step=int(sys.argv[1]), init=True)
    
