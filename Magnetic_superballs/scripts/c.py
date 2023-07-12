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
file='axtest_omega_0.1_F_0.001.txt'
omega=-0.1
N=2
N_virtual=44
lag=0.0
ncol=13
folder=file.split('.txt')[0]
isExist = os.path.exists(folder)
if not isExist:
    os.mkdir(folder)
df = pd.read_csv(file, delimiter = " ", index_col=False, header=None)
data=df.to_numpy()
quat=np.quaternion(1,0,0,0)
quat1=np.quaternion(1,0,0,0)
#mlab.init_notebook(backend='x3d')
#print(data[:,1])
mlab.figure(size=(1920/2,1080/2))
tz = -0.125 * 2. * math.pi;
ty = 0.1312899719849809 * 2. ** math.pi;
tx = 0.25 * 2. *  math.pi;
radius_big_particle=0.5;
radius_midpoints = 2 ** (1 / 4.0) / 6 * 2 * radius_big_particle
radius_coners = 3 ** (1 / 4.0) / 6 * 2 * radius_big_particle
x_midpoints = 1.0 / 2 ** (1 / 4.0) / 3 * 2 * radius_big_particle
x_coner = 1.0 / 3 ** (1 / 4.0) / 3 * 2 * radius_big_particle
radius_red = 0.19782
x_red = 0.27953
y_red = 0.153645
radius_red1 = 0.257166
x_red1 = 0.15215
y_red1 = 0.240383
radius_red2 = 0.237841
x_red2 = 0.178519
y_red2 = 0.2598

Radius = (radius_big_particle,radius_coners,radius_midpoints,radius_red,radius_red1,radius_red1)
colors = np.zeros((6, 3))
colors=((0.0,0.0,1.0),(1.0,0.647,0.0),(1.0,1.0,1.0),(1.0,0.0,1.0),(1.0,0.0,1.0),(1.0,0.0,1.0))
Points = np.zeros((92, 3))
Ptype = np.int_(np.ones(92))
for j in range(8):
    Ptype[j] = 1
for j in range(8, 20):
    Ptype[j] = 2
for j in range(20, 44):
    Ptype[j] = 3
for j in range(44, 68):
    Ptype[j] = 4
for j in range(68, 92):
    Ptype[j] = 5
Points[0] = (x_coner, x_coner, x_coner)
Points[1] = (x_coner, x_coner, -x_coner)
Points[2] = (x_coner, -x_coner, -x_coner)
Points[3] = (x_coner, -x_coner, x_coner)
Points[4] = (-x_coner, -x_coner, -x_coner)
Points[5] = (-x_coner, -x_coner, x_coner)
Points[6] = (-x_coner, x_coner, x_coner)
Points[7] = (-x_coner, x_coner, -x_coner)
Points[8] = (-x_midpoints, 0, -x_midpoints)
Points[9] = (-x_midpoints, 0, x_midpoints)
Points[10] = (x_midpoints, 0, x_midpoints)
Points[11] = (x_midpoints, 0, -x_midpoints)
Points[12] = (0, x_midpoints, x_midpoints)
Points[13] = (0, x_midpoints, -x_midpoints)
Points[14] = (0, -x_midpoints, -x_midpoints)
Points[15] = (0, -x_midpoints, x_midpoints)
Points[16] = (-x_midpoints, x_midpoints, 0)
Points[17] = (x_midpoints, x_midpoints, 0)
Points[18] = (x_midpoints, -x_midpoints, 0)
Points[19] = (-x_midpoints, -x_midpoints, 0)
Points[20] = (x_red, x_red, y_red)
Points[21] = (x_red, x_red, -y_red)
Points[22] = (x_red, -x_red, y_red)
Points[23] = (-x_red, x_red, y_red)
Points[24] = (-x_red, -x_red, y_red)
Points[25] = (x_red, -x_red, -y_red)
Points[26] = (-x_red, x_red, -y_red)
Points[27] = (-x_red, -x_red, -y_red)
Points[28] = (x_red, y_red, x_red)
Points[29] = (x_red, y_red, -x_red)
Points[30] = (x_red, -y_red, x_red)
Points[31] = (-x_red, y_red, x_red)
Points[32] = (-x_red, -y_red, x_red)
Points[33] = (x_red, -y_red, -x_red)
Points[34] = (-x_red, y_red, -x_red)
Points[35] = (-x_red, -y_red, -x_red)
Points[36] = (y_red, x_red, x_red)
Points[37] = (y_red, x_red, -x_red)
Points[38] = (y_red, -x_red, x_red)
Points[39] = (-y_red, x_red, x_red)
Points[40] = (-y_red, -x_red, x_red)
Points[41] = (y_red, -x_red, -x_red)
Points[42] = (-y_red, x_red, -x_red)
Points[43] = (-y_red, -x_red, -x_red)
Points[44] = (x_red1, x_red1, y_red1)
Points[45] = (x_red1, x_red1, -y_red1)
Points[46] = (x_red1, -x_red1, y_red1)
Points[47] = (-x_red1, x_red1, y_red1)
Points[48] = (-x_red1, -x_red1, y_red1)
Points[49] = (x_red1, -x_red1, -y_red1)
Points[50] = (-x_red1, x_red1, -y_red1)
Points[51] = (-x_red1, -x_red1, -y_red1)
Points[52] = (x_red1, y_red1, x_red1)
Points[53] = (x_red1, y_red1, -x_red1)
Points[54] = (x_red1, -y_red1, x_red1)
Points[55] = (-x_red1, y_red1, x_red1)
Points[56] = (-x_red1, -y_red1, x_red1)
Points[57] = (x_red1, -y_red1, -x_red1)
Points[58] = (-x_red1, y_red1, -x_red1)
Points[59] = (-x_red1, -y_red1, -x_red1)
Points[60] = (y_red1, x_red1, x_red1)
Points[61] = (y_red1, x_red1, -x_red1)
Points[62] = (y_red1, -x_red1, x_red1)
Points[63] = (-y_red1, x_red1, x_red1)
Points[64] = (-y_red1, -x_red1, x_red1)
Points[65] = (y_red1, -x_red1, -x_red1)
Points[66] = (-y_red1, x_red1, -x_red1)
Points[67] = (-y_red1, -x_red1, -x_red1)
Points[68] = (x_red2, y_red2, 0)
Points[69] = (y_red2, x_red2, 0)
Points[70] = (x_red2, 0, y_red2)
Points[71] = (y_red2, 0, x_red2)
Points[72] = (0, x_red2, y_red2)
Points[73] = (0, y_red2, x_red2)
Points[74] = (-x_red2, -y_red2, 0)
Points[75] = (-y_red2, -x_red2, 0)
Points[76] = (-x_red2, 0, -y_red2)
Points[77] = (-y_red2, 0, -x_red2)
Points[78] = (0, -x_red2, -y_red2)
Points[79] = (0, -y_red2, -x_red2)
Points[80] = (-x_red2, y_red2, 0)
Points[81] = (-y_red2, x_red2, 0)
Points[82] = (-x_red2, 0, y_red2)
Points[83] = (-y_red2, 0, x_red2)
Points[84] = (0, -x_red2, y_red2)
Points[85] = (0, -y_red2, x_red2)
Points[86] = (x_red2, -y_red2, 0)
Points[87] = (y_red2, -x_red2, 0)
Points[88] = (x_red2, 0, -y_red2)
Points[89] = (y_red2, 0, -x_red2)
Points[90] = (0, x_red2, -y_red2)
Points[91] = (0, y_red2, -x_red2)

m1=R.from_rotvec(tx * np.array([1, 0, 0]))*R.from_rotvec(ty * np.array([0, 1, 0]))*R.from_rotvec(tz * np.array([0, 0, 1]))
mu1=np.array([1, 0, 0])
m2=m1.as_matrix()
for i in range(N_virtual):
    Points[i]=m2.dot(Points[i])
[phi,theta] = np.mgrid[0:2*np.pi:40j,0:np.pi:40j]
x = np.cos(phi)*np.sin(theta)
y = np.sin(phi)*np.sin(theta)
z = np.cos(theta)
#x1,y1=np.mgrid[-2:2:1,-2:2:1]
def plot_sphere(p, col):    
    r,a,b,c = p
    #r=1
    return mlab.mesh(r*x+a, r*y+b, r*z+c,color=col )  
mlab.figure(size=(1920/2,1080/2))
x1,y1=np.mgrid[-2:2:1,-2:2:1]
it=0;
t=data[it,0]
s1=r't='
sfull=f'{s1:2s}{data[it,0]:.4f}'
for k in range(N):
    centrs=(data[it,1+ncol*k],data[it,2+ncol*k],data[it,3+ncol*k])
    q1 = np.quaternion(data[it,7+ncol*k],data[it,4+ncol*k],data[it,5+ncol*k],data[it,6+ncol*k])
    mx1=quaternion.as_rotation_matrix(q1)
    c = (0.5, centrs[0],centrs[1],centrs[2] )
    plot_sphere(c,colors[0])
    mu=mx1.dot(mu1)
    mlab.quiver3d(centrs[0],centrs[1],centrs[2],mu[0],mu[1],mu[2],color=(0.0,1.0,0.0))
    mlab.quiver3d(centrs[0],centrs[1],centrs[2],math.cos(t*omega+lag),math.sin(omega*t+lag),0,color=(1.0,0.0,0.0))
    
    for i in range (N_virtual):
        apoint=mx1.dot(Points[i])+centrs
        c = (Radius[Ptype[i]], apoint[0],apoint[1],apoint[2] )
        plot_sphere(c,colors[Ptype[i]])        

mlab.mesh(x1,y1,0*x1-0.5*math.sqrt(2.0), color=(0.8,1.0,0.7)) #0.5*math.sqrt(2.0)         
mlab.text(0.05,0.05,sfull, width=0.1, color=(0.0,0.0,0.0))
mlab.view(reset_roll=False)
view = mlab.view()
roll = mlab.roll()
for i in range(0,400):
    print("i=%d",i)
    mlab.clf()
    t=data[it,0]
    s1=r't='
    sfull=f'{s1:2s}{data[it,0]:.4f}'
    for k in range(N):
        centrs=(data[it,1+ncol*k],data[it,2+ncol*k],data[it,3+ncol*k])
        q1 = np.quaternion(data[it,7+ncol*k],data[it,4+ncol*k],data[it,5+ncol*k],data[it,6+ncol*k])
        mx1=quaternion.as_rotation_matrix(q1)
        c = (0.5, centrs[0],centrs[1],centrs[2] )
        plot_sphere(c,colors[0])
        mu=mx1.dot(mu1)
        mlab.quiver3d(centrs[0],centrs[1],centrs[2],mu[0],mu[1],mu[2],color=(0.0,1.0,0.0))
        mlab.quiver3d(centrs[0],centrs[1],centrs[2],math.cos(t*omega+lag),math.sin(omega*t+lag),0,color=(1.0,0.0,0.0))
        for i in range (N_virtual):
            apoint=mx1.dot(Points[i])+centrs
            c = (Radius[Ptype[i]], apoint[0],apoint[1],apoint[2] )
            plot_sphere(c,colors[Ptype[i]])        
    mlab.mesh(x1,y1,0*x1-0.5*math.sqrt(2.0), color=(0.8,1.0,0.7)) #0.5*math.sqrt(2.0)         
    mlab.text(0.5,0.85,sfull, width=0.1, color=(0.0,0.0,0.0))
    mlab.view(*view)
    mlab.roll(roll)
    file="%s/anum%05d.png"% (folder, it)
    mlab.savefig(file, magnification=2.0)
    it=it+5
#mlab.close()
