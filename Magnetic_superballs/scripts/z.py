import pyvista as pv
from pyvista import examples
import pandas as pd
import numpy as np
import math as math
import os
import quaternion
from scipy.spatial.transform import Rotation as R
file='axtest_N_2_Nv_92_omega_-0.95_F_0.001_alpha_0.94.txt'
folder=file.split('.txt')[0]
ncol=13
N=2
omega=-0.1
df = pd.read_csv(file, delimiter = " ", index_col=False, header=None)
data=df.to_numpy()
print(data.shape)
it=0;
m0=np.zeros(3)
m0[0]=1
centrs = np.zeros((N, 3))
angle = np.zeros(N)
axis = np.zeros((N, 3))
moments = np.zeros((N, 3))
#momentsx = np.zeros((N, 3))
#momentsx1 = np.zeros((N, 3))
B=np.zeros(3)
quat=[]
def read_data(it):
    global B
    B=(data[it,11],data[it,12],data[it,13]) 
    for k in range(N):
        t=data[it,0]
        centrs[k]=(data[it,1+ncol*k],data[it,2+ncol*k],data[it,3+ncol*k])
        moments[k]=(data[it,8+ncol*k],data[it,9+ncol*k],data[it,10+ncol*k])
        r = R.from_quat((data[it,4+ncol*k],data[it,5+ncol*k],data[it,6+ncol*k],data[it,7+ncol*k]))
        #q1=np.quaternion(data[it,7+ncol*k],data[it,4+ncol*k],data[it,5+ncol*k],data[it,6+ncol*k])
        #ax1=quaternion.as_rotation_matrix(q1)
        #x=r.as_matrix()
        #momentsx[k]=x.dot(m0)
        #momentsx1[k]=ax1.dot(m0)
        temp=r.as_rotvec()
        norm=np.linalg.norm(temp)
        #print(norm)
        angle[k]=np.rad2deg(norm)
        axis[k]=temp*(1./norm)
        #print(moments[k])
        #print(momentsx[k])
        #print(momentsx1[k])
read_data(it)
print(B)
print(data[0,11])
print(data[0,12])
print(data[0,10])
center=np.zeros(3)
center1=np.zeros(3)
#center1=(1.0, 1.0, -0.5)
tz = -0.125 * 2. * 180;
ty = 0.1312899719849809 * 2. * 180;
tx = 0.25 * 2. *  180;

light1 = pv.Light(
    position=(0, 0.0, 100.0),
    focal_point=(0, 0, 0))
a=pv.Superquadric(center=center, scale=(1, 1, 1), size=0.5, theta_roundness=0.5, phi_roundness=0.5, theta_resolution=160, phi_resolution=160, toroidal=False, thickness=0.333333333333)
arr=pv.Arrow(start=center, direction=(1.0, 0.0, 0.0), tip_length=0.1, tip_radius=0.02, tip_resolution=20, shaft_radius=0.005, shaft_resolution=20, scale=1.0)
plane=pv.Plane(center=(0, 0, -0.5*math.sqrt(2.0)), direction=(0, 0, 1), i_size=4, j_size=4, i_resolution=10, j_resolution=10)
a.rotate_z(tz, point=center1, inplace=True).rotate_y(ty, point=center1, inplace=True).rotate_x(tx, point=center1, inplace=True)
alist=[]
arrlist=[]
#arrlist1=[]
Blist=[]
for i in range(N): 
    alist.append(a.rotate_vector(axis[i], angle[i], point=center, inplace=False).translate(centrs[i], inplace=False))
    arrlist.append(pv.Arrow(start=centrs[i], direction=moments[i], tip_length=0.1, tip_radius=0.02, tip_resolution=20, shaft_radius=0.005, shaft_resolution=20, scale=1.0))
    #arrlist1.append(arr.rotate_vector(axis[i], angle[i], point=center, inplace=False).translate(centrs[i], inplace=False))
    Blist.append(pv.Arrow(start=centrs[i], direction=B, tip_length=0.1, tip_radius=0.02, tip_resolution=20, shaft_radius=0.005, shaft_resolution=20, scale=1.0))

#p = pv.Plotter(lighting=None, off_screen=True, window_size=[1000, 1000])
#p = pv.Plotter(lighting='three lights', off_screen=True, window_size=[1000, 1000])
p = pv.Plotter( off_screen=True, window_size=[1504, 1504])
mfile=folder+"_.mp4"
p.open_movie(mfile, framerate=5, quality=10)
#p.add_light(light1)
#p.add_mesh(a.rotate_vector((1, 1, 1), 30, inplace=True), color=[1.0,0.647,0.0], specular=1.0, specular_power=10)
#p.add_mesh(a, color=[1.0,0.647,0.0], specular=1.0, specular_power=10)
for i in range(N):
    p.add_mesh(alist[i], color=[1.0,0.647,0.0], specular=1.0, specular_power=10)
    p.add_mesh(arrlist[i], color=[0.0,1.0,0.0], specular=1.0, specular_power=10)
    #p.add_mesh(arrlist1[i], color=[0.0,0.0,0.1], specular=1.0, specular_power=10)
    p.add_mesh(Blist[i], color=[1.0,0.0,0.0], specular=1.0, specular_power=10)
       
p.add_mesh(pv.Arrow(start=(1.5,1.5,-0.5*math.sqrt(2.0)+0.05), direction=B, tip_length=0.1, tip_radius=0.06, tip_resolution=20, shaft_radius=0.015, shaft_resolution=20, scale=0.3), color=[1.0,0.0,0.0], specular=1.0, specular_power=10)
points = np.array([[1.5,1.5,-0.5*math.sqrt(2.0)+0.05]])
labels = ['B']
p.add_point_labels(points, labels, italic=True, font_size=60, text_color='red', shape_color=[0.8,1.0,0.7], fill_shape=True, margin=0,
                            point_color='red', point_size=2,  shape=None,
                            render_points_as_spheres=True,
                            always_visible=True, shadow=True)     
    
p.view_vector((5.0, 2, 3))
p.camera.zoom(0.75)
p.camera.position=(5.37735183808318, 2.086683310163463, 3.0073711232728435)
p.camera.clipping_range=(1.3886312789079893, 13.079698586438862)
p.camera.focal_point=(0.17591583728790283, 0.006108909845352173, -0.11349047720432281)
print("i=%d",it)
p.add_mesh(plane, lighting=True, color=[0.8,1.0,0.7])
#p.add_floor('-z', lighting=True, color='tan', pad=2.0)
p.enable_shadows()
s1="\n     t="
sfull=f'{s1}{data[it,0]:.4f}'
p.add_text(sfull, position='upper_left', color=[0.8,1.0,0.7] ,shadow=False, font_size=36)

#p.screenshot('airplane.png')
#p.show(screenshot='airplane.png')
p.show(auto_close=False)  # only necessary for an off-screen movie
# Run through each frame
print("!!!!")
print(p.camera.clipping_range)
print(p.camera.focal_point)

p.write_frame() 
it=it+1
for k in range(0,800):
    read_data(it)
    print("i=%d",it)
    for i in range(N): 
        alist[i]=a.rotate_vector(axis[i], angle[i], point=center1, inplace=False).translate(centrs[i])
        arrlist[i]=pv.Arrow(start=centrs[i], direction=moments[i], tip_length=0.1, tip_radius=0.02, tip_resolution=20, shaft_radius=0.005, shaft_resolution=20, scale=1.0)
        Blist[i]=pv.Arrow(start=centrs[i], direction=B, tip_length=0.1, tip_radius=0.02, tip_resolution=20, shaft_radius=0.005, shaft_resolution=20, scale=1.0)
    p.clear_actors()
    for i in range(N):
        p.add_mesh(alist[i], color=[1.0,0.647,0.0], specular=1.0, specular_power=10)
        p.add_mesh(arrlist[i], color=[0.0,1.0,0.0], specular=1.0, specular_power=10)
        p.add_mesh(Blist[i], color=[1.0,0.0,0.0], specular=1.0, specular_power=10)    
    p.add_mesh(pv.Arrow(start=(1.5,1.5,-0.5*math.sqrt(2.0)+0.05), direction=B, tip_length=0.1, tip_radius=0.06, tip_resolution=20, shaft_radius=0.015, shaft_resolution=20, scale=0.3), color=[1.0,0.0,0.0], specular=1.0, specular_power=10)
    p.add_point_labels(points, labels, italic=True, font_size=60, text_color='red', shape_color=[0.8,1.0,0.7], fill_shape=True, margin=0,
                            point_color='red', point_size=2,  shape=None,
                            render_points_as_spheres=True,
                            always_visible=True, shadow=True) 
    #p.view_vector((5.0, 2, 3))
    #p.camera.zoom(0.75)
    p.add_mesh(plane, lighting=True, color=[0.8,1.0,0.7])
    #p.add_floor('-z', lighting=True, color='tan', pad=2.0)
    #p.enable_shadows()
    p.view_vector((5.0, 2, 3))
    p.camera.zoom(0.75)
    p.camera.position=(5.37735183808318, 2.086683310163463, 3.0073711232728435)
    p.camera.clipping_range=(1.3886312789079893, 13.079698586438862)
    p.camera.focal_point=(0.17591583728790283, 0.006108909845352173, -0.11349047720432281)
    s1="\n     t="
    sfull=f'{s1}{data[it,0]:.4f}'
    p.add_text(sfull, position='upper_left', color=[0.8,1.0,0.7],shadow=False, font_size=36)
    p.update()
    p.write_frame()
    it=it+1
    #p.show(screenshot='airplanea.png')
p.close()
