
from __future__ import print_function
import espressomd
from espressomd.magnetostatics import *
from espressomd.galilei import*
from espressomd.analyze import *
from espressomd.virtual_sites import VirtualSitesRelative
from espressomd.io.writer import h5md
from espressomd import visualization
from threading import Thread
from espressomd.shapes import Wall
from espressomd.io.writer import h5md

print(espressomd.features())
required_features = ["LENNARD_JONES"]
espressomd.assert_features(required_features)
import numpy as np
from scipy.spatial.transform import Rotation as Rot
import math
import sys

# Implortant parameter

n_part = 2
n_virtual=0 #n_virtual=92
#H_ext=10.0
phi=12
H_ext=0.500
H_fiel=(H_ext,0.0, 0.0)
g_c=(0.0,-2.0,0.0)
#H_fiel=(H_ext/math.sqrt(2.0),0.0, H_ext/math.sqrt(2.0))
temperature_kT=0.00728
#temperature_kT=0.213843
time_step = 0.0001
dipole_prefactor=1.0
radius_big_particle=0.5
visc=1

#q=2
radius_midpoints=2**(1/4.0)/6*2*radius_big_particle
radius_coners=3**(1/4.0)/6*2*radius_big_particle
x_midpoints=1.0/2**(1/4.0)/3*2*radius_big_particle
x_coner=1.0/3**(1/4.0)/3*2*radius_big_particle
radius_red=0.19782
x_red=0.27953
y_red=0.153645
radius_red1=0.257166
x_red1=0.15215
y_red1=0.240383
radius_red2=0.237841
x_red2=0.178519
y_red2=0.2598
# radius_red2=0.218021
# x_red2=0.179663
# y_red2=0.279588

## Functions

# Function definition is here
def printangle( filename, sys, n_part, H_field):
    f = open(filename, 'a')
    radius_vector=np.zeros((92, 3))
    for i in range(n_part-1):
        radius_vector[i]=np.copy(sys.part.by_id(i+1) .pos)-np.copy(sys.part.by_id(i) .pos)
    radius_vector[n_part-1] = np.copy(sys.part.by_id(n_part-1) .pos) - np.copy(sys.part.by_id(0) .pos)
    mult=1.0
    for i in range(n_part):
        if(np.cross(radius_vector[i], H_field)[2]<0.0):
            mult=-1.0
        f.write(
            "%15g " % (mult*math.acos(np.dot(radius_vector[i], H_field)/(np.linalg.norm(radius_vector[i])*np.linalg.norm(H_field)))*180.0/math.pi)  )
    f.write("\n")
    f.close()
    return;

def printpos( filename, sys, n_part, nvirtual):
    cur_pos = sys.part.all().pos
    if nvirtual==92:
        f = open(filename, "w")
    for i in range(n_part):
        f.write("%15g %15g %15g %15g %i\n" % (cur_pos[i][0], cur_pos[i][1], cur_pos[i][2],radius_big_particle, i))
        for j in range(8):
            f.write("%15g %15g %15g %15g %i\n" % ( cur_pos[i*nvirtual+n_part+j][0], cur_pos[i*nvirtual+n_part+j][1], cur_pos[i*nvirtual+n_part+j][2], radius_coners, i))
        for j in range(8, 20):
            f.write(
                "%15g %15g %15g %15g %i\n" % (cur_pos[i * nvirtual + n_part + j][0], cur_pos[i * nvirtual + n_part + j][1],
                cur_pos[i * nvirtual + n_part + j][2], radius_midpoints, i))
        for j in range(20, 44):
            f.write(
                "%15g %15g %15g %15g %i\n" % (cur_pos[i * nvirtual + n_part + j][0], cur_pos[i * nvirtual + n_part + j][1],
                    cur_pos[i * nvirtual + n_part + j][2], radius_red, i))
        for j in range(44, 68):
            f.write(
                "%15g %15g %15g %15g %i\n" % (cur_pos[i * nvirtual + n_part + j][0], cur_pos[i * nvirtual + n_part + j][1],
                    cur_pos[i * nvirtual + n_part + j][2], radius_red1, i))
        for j in range(68, 92):
            f.write(
                "%15g %15g %15g %15g %i\n" % (cur_pos[i * nvirtual + n_part + j][0], cur_pos[i * nvirtual + n_part + j][1],
                    cur_pos[i * nvirtual + n_part + j][2], radius_red, i))
    f.close()
    return;

# Function definition is here
def printcords( filename, sys, n_part, nvirtual):
    cur_pos = sys.part.all().pos
    f = open(filename, "w")
    for i in range(n_part):
        f.write("%i %15g %15g %15g %15g\n" % (i,cur_pos[i][0], cur_pos[i][1],cur_pos[i][2],i))
    f.write("\n\n")
    for j in range(n_part):
        for i in range(n_part,n_part+4):
            f.write("%i %15g %15g %15g %15g\n" % (i+j*nvirtual,cur_pos[i+j*nvirtual][0], cur_pos[i+j*nvirtual][1],cur_pos[i+j*nvirtual][2],j))
        for i in range(n_part,n_part+1):
            f.write("%i %15g %15g %15g %15g\n" % (i+j*nvirtual,cur_pos[j*nvirtual+i][0], cur_pos[j*nvirtual+i][1],cur_pos[j*nvirtual+i][2],j))
        f.write("\n\n")
        for i in range(n_part+4,n_part+8):
            f.write("%i %15g %15g %15g %15g\n" % (i+j*nvirtual,cur_pos[j*nvirtual+i][0], cur_pos[j*nvirtual+i][1],cur_pos[j*nvirtual+i][2],j))
        for i in range(n_part+4,n_part+5):
            f.write("%i %15g %15g %15g %15g\n" % (i+j*nvirtual,cur_pos[j*nvirtual+i][0], cur_pos[j*nvirtual+i][1],cur_pos[j*nvirtual+i][2],j))
        f.write("\n\n")
        i=n_part
        f.write("%i %15g %15g %15g %15g\n" % (i+j*nvirtual,cur_pos[j*nvirtual+i][0], cur_pos[j*nvirtual+i][1],cur_pos[j*nvirtual+i][2],j))
        i=n_part+1
        f.write("%i %15g %15g %15g %15g\n" % (i+j*nvirtual,cur_pos[j*nvirtual+i][0], cur_pos[j*nvirtual+i][1],cur_pos[j*nvirtual+i][2],j))
        i=n_part+7
        f.write("%i %15g %15g %15g %15g\n" % (i+j*nvirtual,cur_pos[j*nvirtual+i][0], cur_pos[j*nvirtual+i][1],cur_pos[j*nvirtual+i][2],j))
        i=n_part+6
        f.write("%i %15g %15g %15g %15g\n" % (i+j*nvirtual,cur_pos[j*nvirtual+i][0], cur_pos[j*nvirtual+i][1],cur_pos[j*nvirtual+i][2],j))
        i=n_part
        f.write("%i %15g %15g %15g %15g\n" % (i+j*nvirtual,cur_pos[j*nvirtual+i][0], cur_pos[j*nvirtual+i][1],cur_pos[j*nvirtual+i][2],j))
        f.write("\n\n")
        i=n_part+2
        f.write("%i %15g %15g %15g %15g\n" % (i+j*nvirtual,cur_pos[j*nvirtual+i][0], cur_pos[j*nvirtual+i][1],cur_pos[j*nvirtual+i][2],j))
        i=n_part+3
        f.write("%i %15g %15g %15g %15g\n" % (i+j*nvirtual,cur_pos[j*nvirtual+i][0], cur_pos[j*nvirtual+i][1],cur_pos[j*nvirtual+i][2],j))
        i=n_part+5
        f.write("%i %15g %15g %15g %15g\n" % (i+j*nvirtual,cur_pos[j*nvirtual+i][0], cur_pos[j*nvirtual+i][1],cur_pos[j*nvirtual+i][2],j))
        i=n_part+4
        f.write("%i %15g %15g %15g %15g\n" % (i+j*nvirtual,cur_pos[j*nvirtual+i][0], cur_pos[j*nvirtual+i][1],cur_pos[j*nvirtual+i][2],j))
        i=n_part+2
        f.write("%i %15g %15g %15g %15g\n" % (i+j*nvirtual,cur_pos[j*nvirtual+i][0], cur_pos[j*nvirtual+i][1],cur_pos[j*nvirtual+i][2],j))
        f.write("\n\n")
        for i in range(n_part+8,n_part+12):
            f.write("%i %15g %15g %15g %15g\n" % (i+j*nvirtual,cur_pos[j*nvirtual+i][0], cur_pos[j*nvirtual+i][1],cur_pos[j*nvirtual+i][2],j))
        for i in range(n_part+8,n_part+9):
            f.write("%i %15g %15g %15g %15g\n" % (i+j*nvirtual,cur_pos[j*nvirtual+i][0], cur_pos[j*nvirtual+i][1],cur_pos[j*nvirtual+i][2],j))
        f.write("\n\n")
        for i in range(n_part+12,n_part+16):
            f.write("%i %15g %15g %15g %15g\n" % (i+j*nvirtual,cur_pos[j*nvirtual+i][0], cur_pos[j*nvirtual+i][1],cur_pos[j*nvirtual+i][2],j))
        for i in range(n_part+12,n_part+13):
            f.write("%i %15g %15g %15g %15g\n" % (i+j*nvirtual,cur_pos[j*nvirtual+i][0], cur_pos[j*nvirtual+i][1],cur_pos[j*nvirtual+i][2],j))
        f.write("\n\n")
        for i in range(n_part+16,n_part+20):
            f.write("%i %15g %15g %15g %15g\n" % (i+j*nvirtual,cur_pos[j*nvirtual+i][0], cur_pos[j*nvirtual+i][1],cur_pos[j*nvirtual+i][2],j))
        for i in range(n_part+16,n_part+17):
            f.write("%i %15g %15g %15g %15g\n" % (i+j*nvirtual,cur_pos[j*nvirtual+i][0], cur_pos[j*nvirtual+i][1],cur_pos[j*nvirtual+i][2],j))
        f.write("\n\n")
    f.close()
    return;
### printdipoles
def printdip( filename, sys, n_part,radius_big_particle):
    cur_pos = sys.part.by_id([0,n_part]).pos
    dipoles = sys.part.by_id([0,n_part]).dip
    print(dipoles)
    dipoles = dipoles * radius_big_particle
    print(dipoles)    
    f = open(filename, "w")
    for i in range(n_part):
        f.write("%i %15g %15g %15g %15g %15g %15g %i\n" % (i,cur_pos[i][0], cur_pos[i][1],cur_pos[i][2], dipoles[i][0],dipoles[i][1],dipoles[i][2] ,i))
    f.close()
    return


def print_mathematica( filename):
    fout = open(filename, 'w')
    for i in range( n_part):
        fout.write("part%d={{Blue,Ball[{%.15f, %.15f, %.15f}, %.15f]}" % (i,system.part.by_id(i) .pos[0],system.part.by_id(i) .pos[1], system.part.by_id(i) .pos[2], radius_big_particle ))
        for j in range(8):
            fout.write(", {Orange,Ball[{%.15f, %.15f, %.15f}, %.15f]}" % (system.part.by_id(i*92+n_part+j) .pos[0],system.part.by_id(i*92+n_part+j) .pos[1],system.part.by_id(i*92+n_part+j) .pos[2],radius_coners ))
        for j in range(8,20):
            fout.write(", {White,Ball[{%.15f, %.15f, %.15f}, %.15f]}" % (system.part.by_id(i*92+n_part+j) .pos[0],system.part.by_id(i*92+n_part+j) .pos[1],system.part.by_id(i*92+n_part+j) .pos[2], radius_midpoints ))
        for j in range(20, 44):
            fout.write(", {Magenta,Ball[{%.15f, %.15f, %.15f}, %.15f]}" % (system.part.by_id(i*92+n_part+j) .pos[0], system.part.by_id(i*92+n_part+j) .pos[1],system.part.by_id(i*92+n_part+j) .pos[2], radius_red))
        for j in range(44, 68):
            fout.write(", {Red,Ball[{%.15f, %.15f, %.15f}, %.15f]}" % (system.part.by_id(i*92+n_part+j) .pos[0], system.part.by_id(i*92+n_part+j) .pos[1],system.part.by_id(i*92+n_part+j) .pos[2], radius_red1))
        for j in range(68, 92):
            fout.write(", {Green,Ball[{%.15f, %.15f, %.15f}, %.15f]}" % (system.part.by_id(i*92+n_part+j) .pos[0], system.part.by_id(i*92+n_part+j) .pos[1],system.part.by_id(i*92+n_part+j) .pos[2], radius_red2))
#            print("i={0}, j={1}".format(i*92+n_part+j, j))
        fout.write(",{Green, Arrow[{{%.15f, %.15f, %.15f},am*{%.15f, %.15f, %.15f }+{%.15f, %.15f, %.15f }}]}};\n" % (system.part.by_id(i) .pos[0],system.part.by_id(i) .pos[1], system.part.by_id(i) .pos[2],
                                                                                                                      system.part.by_id(i) .dip[0], system.part.by_id(i) .dip[1], system.part.by_id(i) .dip[2],
                                                                                                                      system.part.by_id(i) .pos[0],system.part.by_id(i) .pos[1], system.part.by_id(i) .pos[2]))
    fout.close()
    return



#




#q=1.5
# radius_midpoints=2**(1/6.0)/4*2*radius_big_particle
# radius_coners=3**(1/6.0)/4*2*radius_big_particle
# x_midpoints=1.0/2**(1/3.0)/4*2*radius_big_particle
# x_coner=1.0/3**(1/3.0)/4*2*radius_big_particle
#radius_red=0.342124
#x_red=0.0919776
#y_red=0.155246



####  end important parameters

print('{0} {1} {2} {3}'.format(radius_midpoints, x_midpoints, radius_coners, x_coner))
box_l=np.ones(3)*15
system = espressomd.System(box_l=box_l)
system.periodicity=np.zeros(3)
#system.seed = 123
system.cell_system.skin = 4.0


system.time_step = time_step


system.virtual_sites = espressomd.virtual_sites.VirtualSitesRelative(have_quaternion=False)
#system.virtual_sites = espressomd.virtual_sites.VirtualSitesRelative(have_velocity=False, have_quaternion=False)
# Add particles to the simulation box at random positions
phi0=math.atan(1.0/math.sqrt(2.))*180/math.pi;

#dipole=(1, 1, math.sqrt(2.)*math.tan((phi+phi0)*math.pi/180.));
#dipole1=(-1,- 1, math.sqrt(2.)*math.tan((phi+phi0)*math.pi/180.));

#dipole1=(1, 1, math.sqrt(2.)*math.tan((phi+phi0)*math.pi/180.));
#dipole=(0.5, 0, 0.5);
#dipole=(1, 0, 0);
dipole=(0, 0, 1);
dipole1=(0, 0, 1);
#dipole=(0.524834, 0.279258, 0.804092);
dipole=dipole/np.linalg.norm(dipole)
dipole1=dipole1/np.linalg.norm(dipole1)

if(n_virtual>1):
    Points = np.zeros((92, 3))
    Ptype=np.int_(np.ones(92))
    for j in range(8):
        Ptype[j]=2
    for j in range(8,20):
        Ptype[j]=3
    for j in range(44, 68):
        Ptype[j] = 4
    for j in range(68, 92):
        Ptype[j] = 5
    Points[0]=(x_coner,x_coner,x_coner)
    Points[1]=(x_coner,x_coner,-x_coner)
    Points[2]=(x_coner,-x_coner,-x_coner)
    Points[3]=(x_coner,-x_coner,x_coner)
    Points[4]=(-x_coner,-x_coner,-x_coner)
    Points[5]=(-x_coner,-x_coner,x_coner)
    Points[6]=(-x_coner,x_coner,x_coner)
    Points[7]=(-x_coner,x_coner,-x_coner)
    Points[8]=(-x_midpoints,0,-x_midpoints)
    Points[9] =(-x_midpoints,0,x_midpoints)
    Points[10] =(x_midpoints,0,x_midpoints)
    Points[11] =(x_midpoints,0,-x_midpoints)
    Points[12] =(0,x_midpoints,x_midpoints)
    Points[13] =(0,x_midpoints,-x_midpoints)
    Points[14] =(0,-x_midpoints,-x_midpoints)
    Points[15] =(0,-x_midpoints,x_midpoints)
    Points[16] =(-x_midpoints,x_midpoints,0)
    Points[17] =(x_midpoints,x_midpoints,0)
    Points[18] =(x_midpoints,-x_midpoints,0)
    Points[19] =(-x_midpoints,-x_midpoints,0)
    Points[20] =(x_red, x_red, y_red)
    Points[21] =(x_red, x_red, -y_red)
    Points[22] =(x_red, -x_red, y_red)
    Points[23] =(-x_red, x_red, y_red)
    Points[24] =(-x_red, -x_red, y_red)
    Points[25] =(x_red, -x_red, -y_red)
    Points[26] =(-x_red, x_red, -y_red)
    Points[27] =(-x_red, -x_red, -y_red)
    Points[28] =(x_red, y_red, x_red)
    Points[29] =(x_red, y_red, -x_red)
    Points[30] = (x_red, -y_red, x_red)
    Points[31] =(-x_red, y_red, x_red)
    Points[32] =(-x_red, -y_red, x_red)
    Points[33] =(x_red, -y_red, -x_red)
    Points[34] =(-x_red, y_red, -x_red)
    Points[35] =(-x_red, -y_red, -x_red)
    Points[36] =(y_red, x_red,x_red)
    Points[37] =(y_red, x_red, -x_red)
    Points[38] =(y_red, -x_red,x_red)
    Points[39] =(-y_red, x_red,x_red)
    Points[40] =(-y_red, -x_red,x_red)
    Points[41] =(y_red, -x_red, -x_red)
    Points[42] =(-y_red, x_red, -x_red)
    Points[43] =(-y_red, -x_red, -x_red)
    Points[44] =(x_red1, x_red1, y_red1)
    Points[45] =(x_red1, x_red1, -y_red1)
    Points[46] =(x_red1, -x_red1, y_red1)
    Points[47] =(-x_red1, x_red1, y_red1)
    Points[48] =(-x_red1, -x_red1, y_red1)
    Points[49] =(x_red1, -x_red1, -y_red1)
    Points[50] =(-x_red1, x_red1, -y_red1)
    Points[51] =(-x_red1, -x_red1, -y_red1)
    Points[52] =(x_red1, y_red1, x_red1)
    Points[53] =(x_red1, y_red1, -x_red1)
    Points[54] =(x_red1, -y_red1, x_red1)
    Points[55] =(-x_red1, y_red1, x_red1)
    Points[56] =(-x_red1, -y_red1, x_red1)
    Points[57] =(x_red1, -y_red1, -x_red1)
    Points[58] =(-x_red1, y_red1, -x_red1)
    Points[59] =(-x_red1, -y_red1, -x_red1)
    Points[60] =(y_red1, x_red1,x_red1)
    Points[61] =(y_red1, x_red1, -x_red1)
    Points[62] =(y_red1, -x_red1,x_red1)
    Points[63] =(-y_red1, x_red1,x_red1)
    Points[64] =(-y_red1, -x_red1,x_red1)
    Points[65] =(y_red1, -x_red1, -x_red1)
    Points[66] =(-y_red1, x_red1, -x_red1)
    Points[67] =(-y_red1, -x_red1, -x_red1)
    Points[68] =(x_red2, y_red2, 0)
    Points[69] =(y_red2, x_red2, 0)
    Points[70] =(x_red2, 0, y_red2)
    Points[71] =(y_red2, 0, x_red2)
    Points[72] =(0, x_red2, y_red2)
    Points[73] =(0, y_red2, x_red2)
    Points[74] =(-x_red2, -y_red2, 0)
    Points[75] =(-y_red2, -x_red2, 0)
    Points[76] =(-x_red2, 0, -y_red2)
    Points[77] =(-y_red2, 0, -x_red2)
    Points[78] =(0, -x_red2, -y_red2)
    Points[79] =(0, -y_red2, -x_red2)
    Points[80] =(-x_red2, y_red2, 0)
    Points[81] =(-y_red2, x_red2, 0)
    Points[82] =(-x_red2, 0, y_red2)
    Points[83] =(-y_red2, 0, x_red2)
    Points[84] =(0, -x_red2, y_red2)
    Points[85] =(0, -y_red2, x_red2)
    Points[86] =(x_red2, -y_red2, 0)
    Points[87] =(y_red2, -x_red2, 0)
    Points[88] =(x_red2, 0, -y_red2)
    Points[89] =(y_red2, 0, -x_red2)
    Points[90] =(0, x_red2, -y_red2)
    Points[91] =(0, y_red2, -x_red2)
ty = 0.25*2*np.pi;
ry = (Rot.from_rotvec(ty * np.array([0, 1, 0]))).as_matrix()
dp=ry.dot(dipole);
dp1=ry.dot(dipole1);
print(dp)
print(dp1)


for i in range(n_part):
    rpos1 = np.array([0,  0, 1*i])
    rpos = ry.dot(rpos1)
    if i % 2 == 0:
        system.part.add(id=i, type=0, pos=rpos, dip=dp, gamma=visc, gamma_rot=1, rotation=[1, 1, 1], fix=[1, 1, 1])
    else :
        system.part.add(id=i, type=0, pos=rpos, dip=dp1, gamma=visc, gamma_rot=1, rotation=[1, 1, 1], fix=[1, 1, 1])
    #rpos1=np.array([0.2868926887846606*i, 0.2868926887846606*i, 1*i])
    #rpos = rot.dot(rpos1)
    #print(rpos)
    #system.part.add(id=i,type=0, pos=rpos, dip=dipole, gamma=visc, gamma_rot=1, rotation=[1,1,1], fix=[1,1,1])
    #system.part.add(id=i,type=0, pos=(6*i, 6*i, 6*i), dip=dipole, gamma=1, gamma_rot=1, rotation=[1,1,1])

for i in range(n_part):
    for j in range(n_virtual):
        system.part.add(type=Ptype[j], pos=np.add(system.part.by_id(i).pos, ry.dot(Points[j]))).vs_auto_relate_to(system.part.by_id(i).id)
if(n_virtual>1):
    print_mathematica("cords1a_test.txt")
#sys.exit(0)
#for i in range(n_part):
#    system.part.by_id(i) .quat=(math.cos(math.pi/3),math.sin(math.pi/3),0,0)
#part_pos = system.part[0:n_part].pos
#printcords( "cords1.txt", system, n_part)
#printdip( "dip1.txt", system, n_part,radius_big_particle)
for i in range(n_part):
    print(system.part.by_id(i))

# for i in range(22,46):
#     for j in range(i+1,45):
#         print('i={0} j={1} d={2} '.format(i, j, np.linalg.norm(system.part.by_id(i) .pos-system.part[j].pos)), )
# input()

#Rotate particles
#system.thermostat.set_langevin(kT=100, gamma=visc, gamma_rotation=1.0, act_on_virtual=False, seed=42)
#system.integrator.run(10000)
for i in range(n_part):
  system.part.by_id(i).fix=(0,0,0)
#system.thermostat.set_langevin(kT=0, gamma=visc, gamma_rotation=1.0, act_on_virtual=False)
#system.integrator.run(10)

#printcords( "cords1a.txt", system, n_part)
#printdip( "dip1a.txt", system, n_part,radius_big_particle)

cutoff_par=2**(1.0/6.0)
lj_eps = 1.0
lj_sig00 = 2.0*radius_big_particle/cutoff_par
lj_cut00= cutoff_par*lj_sig00
system.non_bonded_inter[0, 0].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig00,cutoff=lj_cut00, shift='auto')

if(n_virtual>1):
    lj_sig01 = (radius_big_particle+radius_red)/cutoff_par
    lj_sig02 = (radius_big_particle+radius_coners)/cutoff_par
    lj_sig03 = (radius_big_particle+radius_midpoints)/cutoff_par
    lj_sig11 = 2.0*radius_red/cutoff_par
    lj_sig12 = (radius_red+radius_coners)/cutoff_par
    lj_sig13 = (radius_red+radius_midpoints)/cutoff_par
    lj_sig22 = 2.0*radius_coners/cutoff_par
    lj_sig23 = (radius_coners+radius_midpoints)/cutoff_par
    lj_sig33 = 2.0*radius_midpoints/cutoff_par
    lj_sig04 = (radius_big_particle+radius_red1)/cutoff_par
    lj_sig14 = (radius_red+radius_red1)/cutoff_par
    lj_sig24 = (radius_coners+radius_red1)/cutoff_par
    lj_sig34 = (radius_midpoints+radius_red1)/cutoff_par
    lj_sig44 = 2*radius_red1/cutoff_par
    lj_sig05 = (radius_big_particle+radius_red2)/cutoff_par
    lj_sig15 = (radius_red+radius_red2)/cutoff_par
    lj_sig25 = (radius_coners+radius_red2)/cutoff_par
    lj_sig35 = (radius_midpoints+radius_red2)/cutoff_par
    lj_sig45 = (radius_red1+radius_red2)/cutoff_par
    lj_sig55 = 2*radius_red2/cutoff_par


    lj_cut01= cutoff_par*lj_sig01
    lj_cut02= cutoff_par*lj_sig02
    lj_cut03= cutoff_par*lj_sig03
    lj_cut11= cutoff_par*lj_sig11
    lj_cut12= cutoff_par*lj_sig12
    lj_cut13= cutoff_par*lj_sig13
    lj_cut22= cutoff_par*lj_sig22
    lj_cut23= cutoff_par*lj_sig23
    lj_cut33= cutoff_par*lj_sig33
    lj_cut04= cutoff_par*lj_sig04
    lj_cut14= cutoff_par*lj_sig14
    lj_cut24= cutoff_par*lj_sig24
    lj_cut34= cutoff_par*lj_sig34
    lj_cut44= cutoff_par*lj_sig44
    lj_cut05= cutoff_par*lj_sig05
    lj_cut15= cutoff_par*lj_sig15
    lj_cut25= cutoff_par*lj_sig25
    lj_cut35= cutoff_par*lj_sig35
    lj_cut45= cutoff_par*lj_sig45
    lj_cut55= cutoff_par*lj_sig55

    system.non_bonded_inter[0, 1].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig01,cutoff=lj_cut01, shift='auto')
    system.non_bonded_inter[0, 2].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig02,cutoff=lj_cut02, shift='auto')
    system.non_bonded_inter[0, 3].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig03,cutoff=lj_cut03, shift='auto')
    system.non_bonded_inter[1, 1].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig11,cutoff=lj_cut11, shift='auto')
    system.non_bonded_inter[1, 2].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig12,cutoff=lj_cut12, shift='auto')
    system.non_bonded_inter[1, 3].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig13,cutoff=lj_cut13, shift='auto')
    system.non_bonded_inter[2, 2].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig22,cutoff=lj_cut22, shift='auto')
    system.non_bonded_inter[3, 3].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig33,cutoff=lj_cut33, shift='auto')
    system.non_bonded_inter[2, 3].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig23,cutoff=lj_cut23, shift='auto')
    system.non_bonded_inter[0, 4].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig04,cutoff=lj_cut04, shift='auto')
    system.non_bonded_inter[1, 4].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig14,cutoff=lj_cut14, shift='auto')
    system.non_bonded_inter[2, 4].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig24,cutoff=lj_cut24, shift='auto')
    system.non_bonded_inter[3, 4].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig34,cutoff=lj_cut34, shift='auto')
    system.non_bonded_inter[4, 4].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig44,cutoff=lj_cut44, shift='auto')
    system.non_bonded_inter[0, 5].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig05,cutoff=lj_cut05, shift='auto')
    system.non_bonded_inter[1, 5].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig15,cutoff=lj_cut15, shift='auto')
    system.non_bonded_inter[2, 5].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig25,cutoff=lj_cut25, shift='auto')
    system.non_bonded_inter[3, 5].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig35,cutoff=lj_cut35, shift='auto')
    system.non_bonded_inter[4, 5].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig45,cutoff=lj_cut45, shift='auto')
    system.non_bonded_inter[5, 5].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig55,cutoff=lj_cut55, shift='auto')

lj_cap = 10
system.force_cap=lj_cap


system.thermostat.turn_off()
# ### Warmup

# In many cases, including this tutorial, particles are initially placed randomly in the simulation box. It is therefore possible that particles overlap, resulting in a huge repulsive force between them. In this case, integrating the equations of motion would not be numerically stable. Hence, it is necessary to remove this overlap. This is done by limiting the maximum force between two particles, integrating the equations of motion, and increasing the force limit step by step as follows:

# In[11]:


warm_steps  = 1000
warm_n_time = 20 #warm_n_time = 20
min_dist    = 5

print(system.analysis.energy())
print(system.analysis.energy()["total"])
#sys.exit()
system.integrator.set_steepest_descent(f_max=20.0, gamma=visc, max_steps=2000, max_displacement=0.01)
system.integrator.run(20)
print("minimize energy")
print(system.analysis.energy()["total"])

#printcords( "cords1b.txt", system, n_part)
#printdip( "dip1b.txt", system, n_part,radius_big_particle)



#system.thermostat.set_langevin(kT=temperature_kT, gamma=visc, gamma_rotation=1.0, act_on_virtual=False, seed=42)
#system.integrator.set_brownian_dynamics()
system.integrator.set_vv()
omega=0.5
Bmax=10.5#1.05
system.integrator.External_magnetic_Field_par=(Bmax,0.0, 0.0, 0.0, 1.0, 0.1, 6.0, 7.0 ,8.9)
#system.integrator.External_magnetic_Field_par=(1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 6.0, 7.0 ,8.9)
system.integrator.PARAMAGNETIC_DIPOLES=False
system.integrator.EXTERNAL_MAGNETIC_ROTATING_FIELD=True
#system.integrator.External_magnetic_Field_par[1]=2.0;
print("!!!!!!!!!!!!!!!!!!!!!!Look here\n")
print(system.integrator.External_magnetic_Field_par)
print(system.integrator.EXTERNAL_MAGNETIC_ROTATING_FIELD)
print(system.integrator.PARAMAGNETIC_DIPOLES)
#quit()
print(system.analysis.energy()["total"])
print("Warm faze\n")
i = 0
nonbonded_energy = system.analysis.energy()["non_bonded"]
system.force_cap=1
while i < warm_n_time:
   print('{0} {1} {2} {3} {4} {5}'.format(i,np.linalg.norm(system.part.by_id(1) .pos-system.part.by_id(0) .pos),system.analysis.energy()["total"], system.analysis.energy()["dipolar"], system.analysis.energy()["kinetic"],system.analysis.energy()["non_bonded"]))
   system.integrator.run(warm_steps)
   system.galilei.kill_particle_motion(rotation =1)
   nonbonded_energy = system.analysis.energy()["non_bonded"]
   print('{0} {1} {2} {3} {4} {5}'.format(i,np.linalg.norm(system.part.by_id(1) .pos-system.part.by_id(0) .pos),system.analysis.energy()["total"], system.analysis.energy()["dipolar"], system.analysis.energy()["kinetic"],system.analysis.energy()["non_bonded"]))
   #print(system.analysis.energy())
   i+=1
   lj_cap += 1.0
   system.force_cap=lj_cap

print(system.analysis.energy())

system.time_step = time_step
system.time=0.0;
#system.force_cap=20

#system.thermostat.set_langevin(kT=temperature_kT, gamma=visc, gamma_rotation=1.0, act_on_virtual=False, seed=42)
dds_cpu = DipolarDirectSumCpu(prefactor=dipole_prefactor)
system.actors.add(dds_cpu)
#system.integrator.run(steps=0, recalc_forces=True)

#printcords( "cords.txt", system, n_part)
#printdip( "dip.txt", system, n_part,radius_big_particle)

#H_constraint = espressomd.constraints.HomogeneousMagneticField(H=H_fiel)
#gravity = espressomd.constraints.Gravity(g=g_c)


#system.constraints.add(H_constraint)
#system.constraints.add(gravity)
#system.constraints.add(shape=Wall(dist=-4, normal=[0, 1, 0]),  particle_type=4)

# ### Integrating equations of motion and taking measurements

# Once warmup is done, the force capping is switched off by setting it to zero.

# In[12]:

#h5_file = h5md.H5md(filename="sample.h5", write_pos=True, write_vel=True,write_force=True, write_species=True, write_mass=False,write_charge=False, write_ordered=True)
local_temperature_kT=temperature_kT*0.0;
system.thermostat.set_langevin(kT=local_temperature_kT, gamma=visc, gamma_rotation=1.0, act_on_virtual=False, seed=42)
print('{0} {1} {2} {3} {4}'.format(i, system.analysis.energy()["dipolar"], system.analysis.energy()["kinetic"],np.linalg.norm(system.part.by_id(1) .pos-system.part.by_id(0) .pos),np.dot(system.part.by_id(1) .dip,system.part.by_id(0) .dip)))
act_min_dist = system.analysis.energy()["dipolar"]
#while act_min_dist > -20 :
print(system.analysis.energy())
print('{0} {1} {2} {3} {4}'.format(i, system.analysis.energy()["dipolar"], system.analysis.energy()["kinetic"],
                                   np.linalg.norm(system.part.by_id(1) .pos - system.part.by_id(0) .pos),
                                   np.dot(system.part.by_id(1) .dip, system.part.by_id(0) .dip)))
#input()


print("Here!!!")
dip_=system.part.by_id(0).dip; #dip_=dip_/np.linalg.norm(dip_)
dip1_=system.part.by_id(1).dip; #dip1_=dip_/np.linalg.norm(dip1_)
print(dip_)
print(dip1_)



for i in range(10):
    system.integrator.run(1)
    dip_=system.part.by_id(0).dip; #dip_=dip_/np.linalg.norm(dip_)
    dip1_=system.part.by_id(1).dip; #dip1_=dip_/np.linalg.norm(dip1_)
    d1=system.part.by_id(1).pos-system.part.by_id(0).pos
    print(dip_)
    print(dip1_)
    print(d1)
    print('{0} {1} {2} {3}'.format(i, np.linalg.norm(dip_),np.linalg.norm(dip1_), np.linalg.norm(d1)))

#quit()
system.time=0.0;
system.integrator.External_magnetic_Field_par=(Bmax,omega, 0.0, 0.0, 1.0, 0.1, 6.0, 0.1 ,8.9)

#visualizer = visualization.openGLLive(system)
visualizer = visualization.openGLLive(system,
                        director_arrows=True,
                        director_arrows_type_scale=[1.5, 1.0],
                        director_arrows_type_radii=[0.1, 0.4],
                        director_arrows_type_colors=[[1.0, 0, 0], [0, 1.0, 0]])
def main_thread():
    fname = 'angles.txt'
    f = open(fname, "w")
    #f.close()
    for i in range(1000):  #for i in range(100000):
        print('{0} {1} {2} {3} {4}'.format(i, system.analysis.energy()["dipolar"], system.analysis.energy()["kinetic"],np.linalg.norm(system.part.by_id(1) .pos-system.part.by_id(0) .pos),np.dot(system.part.by_id(1) .dip,system.part.by_id(0) .dip)))
        for i in range(100): #for i in range(10):
            system.integrator.run(300) #system.integrator.run(300)
            dip_=system.part.by_id(0).dip; dip_=dip_/np.linalg.norm(dip_)
            dip1_=system.part.by_id(1).dip; dip1_=dip1_/np.linalg.norm(dip1_)
            d1=system.part.by_id(1).pos-system.part.by_id(0).pos
            print('{0} {1} {2} {3} {4} {5} {6} {7} {8} {9}'.format(system.time, dip_[0], dip_[1], dip_[2], dip1_[0] , math.cos(system.time*omega), math.sin(system.time*omega), np.dot(d1, dip_),np.linalg.norm(d1),np.linalg.norm(system.part.by_id(1).v)), file=f)
            #if(n_virtual>1):
 
                #printangle(fname,system, n_part, H_fiel)
        visualizer.update()
        #h5_file.write()
        #print(system.analysis.energy()["dipolar"])
        print("Forces")
        print('{0} {1} {2}'.format(system.part.by_id(0).f[0], system.part.by_id(0).f[1], system.part.by_id(0).f[2]))
        print('{0} {1} {2}'.format(system.part.by_id(1).f[0], system.part.by_id(1).f[1], system.part.by_id(1).f[2]))
        print("dipole")
        print(system.part.by_id(0).dip)
        print(system.part.by_id(1).dip)
        print(np.linalg.norm(system.part.by_id(0).dip))
        print('{0} {1} {2}'.format(system.part.by_id(1).f[0], system.part.by_id(1).f[1], system.part.by_id(1).f[2]))
        d=system.part.by_id(1).pos-system.part.by_id(0).pos
        field=system.part.by_id(0).f
        torque=system.part.by_id(0).torque_lab
        print(d)
        print(np.dot(field/np.linalg.norm(field), d/np.linalg.norm(d)))
        B_=( math.cos(system.time*omega),  math.sin(system.time*omega),0.)
        print(np.dot(B_/np.linalg.norm(B_), d/np.linalg.norm(d)))
        print("AAA")
        print(np.dot(B_, system.part.by_id(0).dip))
        print(torque)
        
        #print('{0} {1} {2}'.format(i, system.analysis.energy()["dipolar"], system.analysis.energy()["kinetic"]),system.time)
        print(system.part.by_id(0).dip)
        print(system.part.by_id(1).dip)
        #print('{0} {1} {2} {3} {4} {5} {6}'.format(i, np.dot(H_fiel,system.part.by_id(0) .dip), np.dot(H_fiel,system.part.by_id(1) .dip), (system.part.by_id(1) .pos[2]-system.part.by_id(0) .pos[2]), (system.part.by_id(1) .pos[0]-system.part.by_id(0) .pos[0]), (system.part.by_id(1) .pos[1]-system.part.by_id(0) .pos[1]), np.dot(system.part.by_id(1) .dip,system.part.by_id(0) .dip)))
        #print('{0} {1} {2}'.format(i, np.dot(H_fiel, system.part.by_id(0) .dip), system.part.by_id(0) .pos[1]))
#,  (system.part.by_id(1) .pos[2]-system.part.by_id(0) .pos[2]), np.linalg.norm(system.part.by_id(1) .pos[0]-system.part.by_id(0) .pos[0]) ,np.dot(system.part.by_id(1) .dip,system.part.by_id(0) .dip)))  
        #string='data/cords_'
        #string +=str(i)
        #string +='.txt'
        #printpos( string, system, n_part, n_virtual)
        #string='data/dip_2d_'
        #string +=str(i)
        #string +='.txt'
        #printdip( string, system, n_part,radius_big_particle)
        global local_temperature_kT 
        #local_temperature_kT/=1.02
        #system.thermostat.set_langevin(kT=local_temperature_kT, gamma=visc, gamma_rotation=1.0, act_on_virtual=False)
        #print('{0} {1} {2} {3}'.format(system.analysis.energy()["dipolar"], system.analysis.energy()["kinetic"]),system.part.by_id(0) .dip,system.part.by_id(1) .dip)    
#    with open('final_pos.txt', 'wb') as f:
#        for i in range(21*n_part):
#            np.savetxt(f, np.transpose(system.part.by_id(i) .pos), fmt='%.15e')
    #print_mathematica("final_pos_mat.txt")
    #h5_file.flush()
    #h5_file.close()
    f.close()
t = Thread(target=main_thread)
t.daemon = True
t.start()
visualizer.start()

print(np.linalg.norm(system.part.by_id(1) .pos-system.part.by_id(0) .pos))
for i in range(n_part):
    print(system.part.by_id(i).pos)
print(system.analysis.energy())

system.force_cap=0


for i in range(n_part+1):
    print(system.part.by_id(i) )
print(system.analysis.energy())
#printcords( "cords2.txt", system, n_part)
#printdip( "dip2.txt", system, n_part,radius_big_particle)
system.integrator.set_steepest_descent(f_max=20000.0, gamma=0.01, max_steps=2000, max_displacement=0.0005)
system.integrator.run(2000)
print("minimize energy")
print(system.analysis.energy()["dipolar"])
