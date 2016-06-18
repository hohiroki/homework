#!/usr/bin/env python
"""
lbm 2dq9
"""
import numpy as np
import theano
import numba

@numba.jit(nopython=True)
def calc_rho(f, rho):
    for i in xrange(rho.shape[0]):
        for j in xrange(rho.shape[1]):
            rho[i][j] = f[0][i][j]+f[1][i][j]+f[2][i][j]+f[3][i][j]+f[4][i][j]+f[5][i][j]+f[6][i][j]+f[7][i][j]+f[8][i][j]


x=theano.tensor.dmatrix('x')
x0=theano.tensor.dmatrix('x0')
x1=theano.tensor.dmatrix('x1')
x2=theano.tensor.dmatrix('x2')
x3=theano.tensor.dmatrix('x3')
x4=theano.tensor.dmatrix('x4')
x5=theano.tensor.dmatrix('x5')
x6=theano.tensor.dmatrix('x6')
x7=theano.tensor.dmatrix('x7')
x8=theano.tensor.dmatrix('x8')
y=theano.tensor.dmatrix('y')
z=theano.tensor.dmatrix('z')
r1=1.+3*x+4.5*x*x-y
ff=theano.function([x,y],r1)
r2=1.5*(x*x+y*y)
ff2=theano.function([x,y],r2)
r3=(x1-x2+x3-x4-x5+x6)/y
r4=(x1-x2+x3+x4-x5-x6)/y
ff3=theano.function([x1,x2,x3,x4,x5,x6,y],r3)
ff4=theano.function([x1,x2,x3,x4,x5,x6,y],r4)
r5=x0+x1+x2+x3+x4+x5+x6+x7+x8
ff5=theano.function([x0,x1,x2,x3,x4,x5,x6,x7,x8],r5)
a = theano.tensor.scalar('a')
b = theano.tensor.scalar('b')
r6=a*x*y+b*z
gg=theano.function([a,b,x,y,z],r6)

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
    def boundary_evolution(self):
        """ """
        nx,ny = self.NX,self.NY
        w = self.w
        self.rho[nx,1:-1] = self.rho[nx-1,1:-1]
        self.rho[0,1:-1] = self.rho[1,1:-1]
        u,v = self.u[nx],self.v[nx]
        u1,v1 = self.u[nx-1],self.v[nx-1]
        r,r1 = self.rho[nx],self.rho[nx-1]
        uv,uv1 = u*u + v*v, u1*u1 + v1*v1
        self.f[0,nx] = self.f[0,nx-1]+w[0]*(r*(1-1.5*uv)-r1*(1.-1.5*uv1))
        self.f[1,nx] = self.f[1,nx-1]+w[1]*(r*(1+3*u+4.5*u*u-1.5*uv)-r1*(1+3*u1+4.5*u1*u1-1.5*uv1))
        self.f[2,nx] = self.f[2,nx-1]+w[2]*(r*(1+3*v+4.5*v*v-1.5*uv)-r1*(1+3*v1+4.5*v1*v1-1.5*uv1))
        self.f[3,nx] = self.f[3,nx-1]+w[3]*(r*(1-3*u+4.5*u*u-1.5*uv)-r1*(1-3*u1+4.5*u1*u1-1.5*uv1))
        self.f[4,nx] = self.f[4,nx-1]+w[4]*(r*(1-3*v+4.5*v*v-1.5*uv)-r1*(1-3*v1+4.5*v1*v1-1.5*uv1))
        eu, eu1 = u+v, u1+v1
        self.f[5,nx] = self.f[5,nx-1]+w[5]*(r*(1+3*eu+4.5*eu*eu-1.5*uv)-r1*(1+3*eu1+4.5*eu1*eu1-1.5*uv1))
        eu, eu1 = -u+v, -u1+v1
        self.f[6,nx] = self.f[6,nx-1]+w[6]*(r*(1+3*eu+4.5*eu*eu-1.5*uv)-r1*(1+3*eu1+4.5*eu1*eu1-1.5*uv1))
        eu, eu1 = -u-v, -u1-v1
        self.f[7,nx] = self.f[7,nx-1]+w[7]*(r*(1+3*eu+4.5*eu*eu-1.5*uv)-r1*(1+3*eu1+4.5*eu1*eu1-1.5*uv1))
        eu, eu1 = u-v, u1-v1
        self.f[8,nx] = self.f[8,nx-1]+w[8]*(r*(1+3*eu+4.5*eu*eu-1.5*uv)-r1*(1+3*eu1+4.5*eu1*eu1-1.5*uv1))
        u, v = self.u[0], self.v[0]
        u1, v1 = self.u[1], self.v[1]
        r,r1 = self.rho[0],self.rho[1]
        uv,uv1 = u*u + v*v, u1*u1 + v1*v1
        self.f[0,0] = self.f[0,1]+w[0]*(r*(1-1.5*uv)-r1*(1.-1.5*uv1))
        self.f[1,0] = self.f[1,1]+w[1]*(r*(1+3*u+4.5*u*u-1.5*uv)-r1*(1+3*u1+4.5*u1*u1-1.5*uv1))
        self.f[2,0] = self.f[2,1]+w[2]*(r*(1+3*v+4.5*v*v-1.5*uv)-r1*(1+3*v1+4.5*v1*v1-1.5*uv1))
        self.f[3,0] = self.f[3,1]+w[3]*(r*(1-3*u+4.5*u*u-1.5*uv)-r1*(1-3*u1+4.5*u1*u1-1.5*uv1))
        self.f[4,0] = self.f[4,1]+w[4]*(r*(1-3*v+4.5*v*v-1.5*uv)-r1*(1-3*v1+4.5*v1*v1-1.5*uv1))
        eu, eu1 = u+v, u1+v1
        self.f[5,0] = self.f[5,1]+w[5]*(r*(1+3*eu+4.5*eu*eu-1.5*uv)-r1*(1+3*eu1+4.5*eu1*eu1-1.5*uv1))
        eu, eu1 = -u+v, -u1+v1
        self.f[6,0] = self.f[6,1]+w[6]*(r*(1+3*eu+4.5*eu*eu-1.5*uv)-r1*(1+3*eu1+4.5*eu1*eu1-1.5*uv1))
        eu, eu1 = -u-v, -u1-v1
        self.f[7,0] = self.f[7,1]+w[7]*(r*(1+3*eu+4.5*eu*eu-1.5*uv)-r1*(1+3*eu1+4.5*eu1*eu1-1.5*uv1))
        eu, eu1 = u-v, u1-v1
        self.f[8,0] = self.f[8,1]+w[8]*(r*(1+3*eu+4.5*eu*eu-1.5*uv)-r1*(1+3*eu1+4.5*eu1*eu1-1.5*uv1))

        self.rho[:,0] = self.rho[:,1]
        self.rho[:,ny] = self.rho[:,ny-1]
        u, v = self.u[:,0], self.v[:,0]
        u1, v1 = self.u[:,1], self.v[:,1]
        r,r1 = self.rho[:,0], self.rho[:,1]
        uv,uv1 = u*u + v*v, u1*u1 + v1*v1
        self.f[0,:,0] = self.f[0,:,1]+w[0]*(r*(1-1.5*uv)-r1*(1.-1.5*uv1))
        self.f[1,:,0] = self.f[1,:,1]+w[1]*(r*(1+3*u+4.5*u*u-1.5*uv)-r1*(1+3*u1+4.5*u1*u1-1.5*uv1))
        self.f[2,:,0] = self.f[2,:,1]+w[2]*(r*(1+3*v+4.5*v*v-1.5*uv)-r1*(1+3*v1+4.5*v1*v1-1.5*uv1))
        self.f[3,:,0] = self.f[3,:,1]+w[3]*(r*(1-3*u+4.5*u*u-1.5*uv)-r1*(1-3*u1+4.5*u1*u1-1.5*uv1))
        self.f[4,:,0] = self.f[4,:,1]+w[4]*(r*(1-3*v+4.5*v*v-1.5*uv)-r1*(1-3*v1+4.5*v1*v1-1.5*uv1))
        eu, eu1 = u+v, u1+v1
        self.f[5,:,0] = self.f[5,:,1]+w[5]*(r*(1+3*eu+4.5*eu*eu-1.5*uv)-r1*(1+3*eu1+4.5*eu1*eu1-1.5*uv1))
        eu, eu1 = -u+v, -u1+v1
        self.f[6,:,0] = self.f[6,:,1]+w[6]*(r*(1+3*eu+4.5*eu*eu-1.5*uv)-r1*(1+3*eu1+4.5*eu1*eu1-1.5*uv1))
        eu, eu1 = -u-v, -u1-v1
        self.f[7,:,0] = self.f[7,:,1]+w[7]*(r*(1+3*eu+4.5*eu*eu-1.5*uv)-r1*(1+3*eu1+4.5*eu1*eu1-1.5*uv1))
        eu, eu1 = u-v, u1-v1
        self.f[8,:,0] = self.f[8,:,1]+w[8]*(r*(1+3*eu+4.5*eu*eu-1.5*uv)-r1*(1+3*eu1+4.5*eu1*eu1-1.5*uv1))

        u, v = self.u[:,ny], self.v[:,ny]
        u1, v1 = self.u[:,ny-1], self.v[:,ny-1]
        r,r1 = self.rho[:,ny], self.rho[:,ny-1]
        uv,uv1 = u*u + v*v, u1*u1 + v1*v1
        self.f[0,:,ny] = self.f[0,:,ny-1]+w[0]*(r*(1-1.5*uv)-r1*(1.-1.5*uv1))
        self.f[1,:,ny] = self.f[1,:,ny-1]+w[1]*(r*(1+3*u+4.5*u*u-1.5*uv)-r1*(1+3*u1+4.5*u1*u1-1.5*uv1))
        self.f[2,:,ny] = self.f[2,:,ny-1]+w[2]*(r*(1+3*v+4.5*v*v-1.5*uv)-r1*(1+3*v1+4.5*v1*v1-1.5*uv1))
        self.f[3,:,ny] = self.f[3,:,ny-1]+w[3]*(r*(1-3*u+4.5*u*u-1.5*uv)-r1*(1-3*u1+4.5*u1*u1-1.5*uv1))
        self.f[4,:,ny] = self.f[4,:,ny-1]+w[4]*(r*(1-3*v+4.5*v*v-1.5*uv)-r1*(1-3*v1+4.5*v1*v1-1.5*uv1))
        eu, eu1 = u+v, u1+v1
        self.f[5,:,ny] = self.f[5,:,ny-1]+w[5]*(r*(1+3*eu+4.5*eu*eu-1.5*uv)-r1*(1+3*eu1+4.5*eu1*eu1-1.5*uv1))
        eu, eu1 = -u+v, -u1+v1
        self.f[6,:,ny] = self.f[6,:,ny-1]+w[6]*(r*(1+3*eu+4.5*eu*eu-1.5*uv)-r1*(1+3*eu1+4.5*eu1*eu1-1.5*uv1))
        eu, eu1 = -u-v, -u1-v1
        self.f[7,:,ny] = self.f[7,:,ny-1]+w[7]*(r*(1+3*eu+4.5*eu*eu-1.5*uv)-r1*(1+3*eu1+4.5*eu1*eu1-1.5*uv1))
        eu, eu1 = u-v, u1-v1
        self.f[8,:,ny] = self.f[8,:,ny-1]+w[8]*(r*(1+3*eu+4.5*eu*eu-1.5*uv)-r1*(1+3*eu1+4.5*eu1*eu1-1.5*uv1))

    #@profile
    def evolution(self):
        """ """
        nx, ny = self.NX, self.NY
        u, v = self.u, self.v
        uv = ff2(u,v)
        upv = u + v
        umv = u - v
        qq = 1./self.tau_f
        tt = (1. - qq)        
        #e = [[0,0],[1,0],[0,1],[-1,0],[0,-1],[1,1],[-1,1],[-1,-1],[1,-1]]
        G = (1. - uv)
        F = gg(qq*self.w[0], tt, self.rho, G, self.f[0])
        self.f[0,1:-1,1:-1] = F[1:-1,1:-1]
        G = ff(u,uv)
        F = gg(qq*self.w[1], tt, self.rho, G, self.f[1])
        self.f[1,1:-1,1:-1] = F[0:-2,1:-1]
        G -= 6.*u
        F = gg(qq*self.w[3], tt, self.rho, G, self.f[3])
        self.f[3,1:-1,1:-1] = F[2:,1:-1]
        G = ff(v,uv)
        F = gg(qq*self.w[2], tt, self.rho, G, self.f[2])
        self.f[2,1:-1,1:-1] = F[1:-1,0:-2]
        G -= 6.*v
        F = gg(qq*self.w[4], tt, self.rho, G, self.f[4])
        self.f[4,1:-1,1:-1] = F[1:-1,2:]
        G = ff(upv,uv)
        F = gg(qq*self.w[5], tt, self.rho, G, self.f[5])
        self.f[5,1:-1,1:-1] = F[0:-2,0:-2]
        G -= 6.*upv
        F = gg(qq*self.w[7], tt, self.rho, G, self.f[7])
        self.f[7,1:-1,1:-1] = F[2:,2:]
        G = ff(umv,uv)
        F = gg(qq*self.w[8], tt, self.rho, G, self.f[8])
        self.f[8,1:-1,1:-1] = F[0:-2,2:]
        G -= 6.*umv
        F = gg(qq*self.w[6], tt, self.rho, G, self.f[6])
        self.f[6,1:-1,1:-1] = F[2:,0:-2]
        
        self.u0 = self.u
        self.v0 = self.v
        #self.rho=ff5(self.f[0],self.f[1],self.f[2],self.f[3],self.f[4],self.f[5],self.f[6],self.f[7],self.f[8])
        calc_rho(self.f, self.rho)
        self.u = ff3(self.f[1],self.f[3],self.f[5],self.f[6],self.f[7],self.f[8],self.rho)
        self.v = ff4(self.f[2],self.f[4],self.f[5],self.f[6],self.f[7],self.f[8],self.rho)        
        # boundary
        self.set_u_boundary()
        self.boundary_evolution()
            
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
    lbm = LBM2DQ9(256,256)
    if len(sys.argv) > 1:
        #print sys.argv[1]
        lbm.run(Re=1000., step=int(sys.argv[1]), init=True)
    else:
        lbm.run(Re=1000., init=True, output=True)
    
