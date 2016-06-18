import numpy as np
import matplotlib.pyplot as plt

class Topflow(object):
    def __init__(self):
        self.maxstep = 1000000
        self.step = 0
        self.Q = 9
        self.NX = 156
        self.NY = 156
        self.U = 0.1
        self.e = np.array([[0, 0],
                           [1, 0],
                           [0, 1],
                           [-1, 0],
                           [0, -1],
                           [1, 1],
                           [-1, 1],
                           [-1, -1],
                           [1, -1]])
        self.w = np.array([4./9, 1./9, 1./9, 1./9, 1./9, 1./36, 1./36, 1./36, 1./36])
        self.rho = np.ones((self.NX+1, self.NY+1))
        self.u = np.ones((self.NX+1, self.NY+1, 2))
        self.u0 = np.ones((self.NX+1, self.NY+1, 2))
        self.f = np.ones((self.NX+1, self.NY+1, self.Q))
        self.F = np.ones((self.NX+1, self.NY+1, self.Q))
        self.dx, self.dt = 1., 1.
        self.dy = 1.
        self.Lx = self.dx*self.NY
        self.Ly = self.dy*self.NX
        self.c = self.dx/self.dt
        self.rho0 = 1.
        self.Re = 1000
        self.niu = self.U*self.Lx/self.Re
        self.tau_f = 3.*self.niu+0.5

        self.massage = "not start or done yet"
        self.reError = 100
        self.erroControl = 1.E-6

    def initial(self):
        for i in range(0, self.NX+1, 1):
            for j in range(0, self.NY+1, 1):
                self.u[i][j][0] = 0.01
                self.u[i][j][1] = 0.01
                self.rho[i][j] = self.rho0
                self.u[i][self.NY][0] = self.U
                for k in range(0, self.Q, 1):
                    self.f[i][j][k] = self.feq(k, self.rho[i][j], self.u[i][j])


    def feq(self, k, rho, u):
        eu = (self.e[k][0]*u[0]+self.e[k][1]*u[1])
        uv = (u[0]**2+u[1]**2)
        return self.w[k]*rho*(1.+3.*eu+4.5*eu*eu-1.5*uv)

    def evolution(self):
        for i in range(1, self.NX, 1):
            for j in range(1,self.NY, 1):
                for k in range(0, self.Q, 1):
                    ip = i - self.e[k][0]
                    jp = j - self.e[k][1]
                    self.F[i][j][k] = self.f[ip][jp][k] + (self.feq(k, self.rho[ip][jp], self.u[ip][jp])-self.f[ip][jp][k])/self.tau_f

        for i in range(1, self.NX, 1):
            for j in range(1, self.NY, 1):
                self.u0[i][j][0] = self.u[i][j][0]
                self.u0[i][j][1] = self.u[i][j][1]
                self.rho[i][j] = 0
                self.u[i][j][0] = 0
                self.u[i][j][1] = 0
                for k in range(0, self.Q, 1):
                    self.f[i][j][k] = self.F[i][j][k]
                    self.rho[i][j] += self.f[i][j][k]
                    self.u[i][j][0] += self.e[k][0]*self.f[i][j][k]
                    self.u[i][j][1] += self.e[k][1]*self.f[i][j][k]
                self.u[i][j][0] = self.u[i][j][0]/self.rho[i][j]
                self.u[i][j][1] = self.u[i][j][1]/self.rho[i][j]

        #lef and right boundary
        for j in range(1, self.NY, 1):
            for k in range(0, self.Q, 1):
                self.rho[self.NX][j] = self.rho[self.NX-1][j]
                self.f[self.NX][j][k] = self.feq(k, self.rho[self.NX][j], self.u[self.NX][j]) + self.f[self.NX-1][j][k]- self.feq(k, self.rho[self.NX-1][j], self.u[self.NX-1][j])
                self.rho[0][j] = self.rho[1][j]
                self.f[0][j][k] = self.feq(k, self.rho[0][j], self.u[0][j]) + self.f[1][j][k]- self.feq(k, self.rho[1][j], self.u[1][j])

        #top and bottom boundary
        for i in range(0, self.NX+1, 1):
            for k in range(0, self.Q, 1):
                self.rho[i][0] = self.rho[i][1]
                self.f[i][0][k] = self.feq(k, self.rho[i][0], self.u[i][0]) + self.f[i][1][k]- self.feq(k, self.rho[i][1], self.u[i][1])
                self.rho[i][self.NY] = self.rho[1][self.NY-1]
                self.u[i][self.NY][0] = self.U
                self.f[i][self.NY][k] = self.feq(k, self.rho[i][self.NY], self.u[i][self.NY]) + self.f[i][self.NY-1][k]- self.feq(k, self.rho[i][self.NY-1], self.u[i][self.NY-1])

    #output
    def output(self):
        print "The residual error :"
        print self.reError
        #draw figure

    #residual error
    def reErroCacu(self):
        du = self.u0 - self.u
        total_erro = (du[:][:][0]*du[:][:][0]).sum()+(du[:][:][1]*du[:][:][1]).sum()
        uu2 = (self.u[:][:][0]*self.u[:][:][0]).sum() + (self.u[:][:][1]*self.u[:][:][1]).sum()
        self.reError = np.sqrt(total_erro/uu2)
    #index
    def main(self):
        self.initial()
        for n in range(1, self.maxstep, 1):
            self.evolution()
            if(n%10 == 0):
                self.reErroCacu()
                self.output()
                if(self.reError < self.erroControl):
                    break
        if(self.step >=(self.maxstep-2)):
            self.massage = "caculate step over maxstep"
        else:
            self.massage = "caculate done"

if __name__ == '__main__':
    sim = Topflow()
    sim.main()


'''
Y, X = np.mgrid[-3:3:100j, -3:3:100j]
U = -1 - X**2 + Y
V = 1 + X - Y**2
speed = np.sqrt(U*U + V*V)

fig0, ax0 = plt.subplots()
strm = ax0.streamplot(X, Y, U, V, color=U, linewidth=2, cmap=plt.cm.autumn)
fig0.colorbar(strm.lines)

fig1, (ax1, ax2) = plt.subplots(ncols=2)
ax1.streamplot(X, Y, U, V, density=[0.5, 1])

lw = 5*speed / speed.max()
ax2.streamplot(X, Y, U, V, density=0.6, color='k', linewidth=lw)

plt.show()
'''