from mayavi import mlab
import pandas as pd
import numpy as np
import math as math
import os
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from matplotlib import animation
import quaternion
from scipy.spatial.transform import Rotation as R
q1 = np.quaternion(4.0599e-01,7.5234e-01, -5.0865e-01, -1.0219e-01)
file='axtest_omega_0.1_F_0.001.txt'
ncol=13
omega=-0.1
lag=0.0
N=1
N_virtual=44
ac=1.98
folder=file.split('.txt')[0]
isExist = os.path.exists(folder)
if not isExist:
    os.mkdir(folder)
df = pd.read_csv(file, delimiter = " ", index_col=False, header=None)
data=df.to_numpy()
it=0;
s1=r't='
sfull=f'{s1:2s}{data[it,0]:.4f}'
centrs = np.zeros((N, 3))
moments = np.zeros((N, 3))
momentsx = np.zeros((N, 3))
quat=[]
t=0
for i in range(N):
    quat.append(q1)
def read_data(it):
    for k in range(N):
        t=data[it,0]
        centrs[k]=(data[it,1+ncol*k],data[it,2+ncol*k],data[it,3+ncol*k])
        quat[k]=np.quaternion(data[it,7+ncol*k],data[it,4+ncol*k],data[it,5+ncol*k],data[it,6+ncol*k])
        moments[k]=(data[it,8+ncol*k],data[it,9+ncol*k],data[it,10+ncol*k])
read_data(it)
x, y, z = np.mgrid[-2:2:100j, -2:2:100j, -2:2:100j]
x1,y1=np.mgrid[-4:4:1,-4:4:1]
z1=0*x1-0.5*math.sqrt(2.0)
q1 = np.quaternion(4.0599e-01,7.5234e-01, -5.0865e-01, -1.0219e-01)
q2 = np.quaternion(4.0599e-01,-7.5234e-01, 5.0865e-01, 1.0219e-01)
tz = -0.125 * 2. * math.pi;
ty = 0.1312899719849809 * 2. ** math.pi;
tx = 0.25 * 2. *  math.pi;
m1=(R.from_rotvec(tx * np.array([1, 0, 0]))*R.from_rotvec(ty * np.array([0, 1, 0]))*R.from_rotvec(tz * np.array([0, 0, 1]))).as_matrix()
mu1=np.array([1, 0, 0])
m0=np.linalg.inv(m1)
#m=m1.as_matrix()
mu0=mu1
#print(m0)
B0=(math.cos(t*omega+lag),math.sin(omega*t+lag),0);
mlab.mesh(m0[0,0]*x1+m0[0,1]*y1+m0[0,2]*z1,m0[1,0]*x1+m0[1,1]*y1+m0[1,2]*z1,m0[2,0]*x1+m0[2,1]*y1+m0[2,2]*z1, color=(0.8,1.0,0.7)) #0.5*math.sqrt(2.0)
B=m1.dot(B0)
Bx=m0.dot(B0)
print(m1.dot(B0))
print(m0.dot(B0))
for k in range(N):
    momentsx[k]=m0.dot(moments[k])
    X0=centrs[k]
    X0=np.zeros(3)#coment out
    NX0=(m0.dot(X0))
    ax1=quaternion.as_rotation_matrix(quat[k])
    m=ax1
    mx=ax1.transpose()
    mu=m0.dot(ax1.dot(mu0))
    F = (ac*(m[0,0]*(x-NX0[0])+m[0,1]*(y-NX0[1])+m[0,2]*(z-NX0[2])))**4 + (ac*(m[1,0]*(x-NX0[0])+m[1,1]*(y-NX0[1])+m[1,2]*(z-NX0[2])))**4 + (ac*(m[2,0]*(x-NX0[0])+m[2,1]*(y-NX0[1])+m[2,2]*(z-NX0[2])))**4 -1
    F1 = (ac*(mx[0,0]*(x-NX0[0])+mx[0,1]*(y-NX0[1])+mx[0,2]*(z-NX0[2])))**4 + (ac*(mx[1,0]*(x-NX0[0])+mx[1,1]*(y-NX0[1])+mx[1,2]*(z-NX0[2])))**4 + (ac*(mx[2,0]*(x-NX0[0])+mx[2,1]*(y-NX0[1])+mx[2,2]*(z-NX0[2])))**4 -1
    mlab.contour3d(x,y,z,F, contours = [0],color=(1.0,0.647,0.0))
    mlab.contour3d(x,y,z,F1, contours = [0],color=(1.0,0.0,0.647))
    mlab.quiver3d(NX0[0],NX0[1],NX0[2],mu[0],mu[1],mu[2],color=(0.0,1.0,0.0))
    mlab.quiver3d(NX0[0],NX0[1],NX0[2],momentsx[k][0],momentsx[k][1],momentsx[k][2],color=(0.0,0.0,1.0))
    #mlab.quiver3d(NX0[0],NX0[1],NX0[2],B[0],B[1],B[2],color=(1,0.0,0.0))
    mlab.quiver3d(NX0[0],NX0[1],NX0[2],Bx[0],Bx[1],Bx[2],color=(1.0,0.0,0.0))
mlab.text(0.05,0.95,sfull, width=0.1, color=(0.0,0.0,0.0))
#exec(open("a.py").read())
view =mlab.view(83.35643132532017, 90.82451424273874, 10.753908743193113, (-0.31768616,  0.68231384, -0.65842488), 19.727559240933086)
mlab.view(reset_roll=False)
for i in range(0,20):
    print("i=%d",i)
    mlab.clf()
    t=data[it,0]
    s1=r't='
    sfull=f'{s1:2s}{data[it,0]:.4f}'
    read_data(it)
    B0=(math.cos(t*omega+lag),math.sin(omega*t+lag),0);
    mlab.mesh(m0[0,0]*x1+m0[0,1]*y1+m0[0,2]*z1,m0[1,0]*x1+m0[1,1]*y1+m0[1,2]*z1,m0[2,0]*x1+m0[2,1]*y1+m0[2,2]*z1, color=(0.8,1.0,0.7)) #0.5*math.sqrt(2.0)
    B=m1.dot(B0)
    Bx=m0.dot(B0)
    for k in range(N):
        momentsx[k]=m0.dot(moments[k])
        X0=centrs[k]
        X0=np.zeros(3)#coment out
        NX0=(m0.dot(X0))
        ax1=quaternion.as_rotation_matrix(quat[k])
        #m=ax1.transpose()
        m=ax1
        #m=np.identity(3)
        mu=m0.dot(ax1.dot(mu0))
        F = (ac*(m[0,0]*(x-NX0[0])+m[0,1]*(y-NX0[1])+m[0,2]*(z-NX0[2])))**4 + (ac*(m[1,0]*(x-NX0[0])+m[1,1]*(y-NX0[1])+m[1,2]*(z-NX0[2])))**4 + (ac*(m[2,0]*(x-NX0[0])+m[2,1]*(y-NX0[1])+m[2,2]*(z-NX0[2])))**4 -1
        mlab.contour3d(x,y,z,F, contours = [0],color=(1.0,0.647,0.0))
        mlab.quiver3d(NX0[0],NX0[1],NX0[2],mu[0],mu[1],mu[2],color=(0.0,1.0,0.0))
        mlab.quiver3d(NX0[0],NX0[1],NX0[2],momentsx[k][0],momentsx[k][1],momentsx[k][2],color=(0.0,0.0,1.0))
        #mlab.quiver3d(NX0[0],NX0[1],NX0[2],B[0],B[1],B[2],color=(1,0.0,0.0))
        mlab.quiver3d(NX0[0],NX0[1],NX0[2],Bx[0],Bx[1],Bx[2],color=(1.0,0.0,0.0))
    mlab.text(0.05,0.95,sfull, width=0.1, color=(0.0,0.0,0.0))
    mlab.view(83.35643132532017, 90.82451424273874, 10.753908743193113, (-0.31768616,  0.68231384, -0.65842488), 19.727559240933086)
    file="%s/anum%05d.png"% (folder, it)
    mlab.savefig(file, magnification=2.0)
    it=it+5
#mlab.close() 

