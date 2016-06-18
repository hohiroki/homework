%物性与，初始条件
global gama T0 P0 R Pref Rouref Uref Lref eps epsmax dt

gama = 1.26;
T0 = 3200;
P0 = 6*1e6;
R = 287.314;
Pref = 1e5;
Rouref = 1;
Uref = (Pref/Rouref);
Lref = 5e-3;
eps = 1;
epsmax = 1e-6;
dt = 1/(xnum*ynum*10);