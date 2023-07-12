
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

print(espressomd.features())
required_features = ["LENNARD_JONES"]
espressomd.assert_features(required_features)
import numpy as np
import math

## Functions

# Function definition is here
def printcords( filename, sys, n_part):
    cur_pos = sys.part[:].pos
    f = open(filename, "w")
    for i in range(n_part):
        f.write("%i %15g %15g %15g %15g\n" % (i,cur_pos[i][0], cur_pos[i][1],cur_pos[i][2],i))
    f.write("\n\n")
    for j in range(n_part):
        for i in range(n_part,n_part+4):
            f.write("%i %15g %15g %15g %15g\n" % (i+j*20,cur_pos[i+j*20][0], cur_pos[i+j*20][1],cur_pos[i+j*20][2],j))
        for i in range(n_part,n_part+1):
            f.write("%i %15g %15g %15g %15g\n" % (i+j*20,cur_pos[j*20+i][0], cur_pos[j*20+i][1],cur_pos[j*20+i][2],j))
        f.write("\n\n")
        for i in range(n_part+4,n_part+8):
            f.write("%i %15g %15g %15g %15g\n" % (i+j*20,cur_pos[j*20+i][0], cur_pos[j*20+i][1],cur_pos[j*20+i][2],j))
        for i in range(n_part+4,n_part+5):
            f.write("%i %15g %15g %15g %15g\n" % (i+j*20,cur_pos[j*20+i][0], cur_pos[j*20+i][1],cur_pos[j*20+i][2],j))
        f.write("\n\n")
        i=n_part
        f.write("%i %15g %15g %15g %15g\n" % (i+j*20,cur_pos[j*20+i][0], cur_pos[j*20+i][1],cur_pos[j*20+i][2],j))
        i=n_part+1
        f.write("%i %15g %15g %15g %15g\n" % (i+j*20,cur_pos[j*20+i][0], cur_pos[j*20+i][1],cur_pos[j*20+i][2],j))
        i=n_part+7
        f.write("%i %15g %15g %15g %15g\n" % (i+j*20,cur_pos[j*20+i][0], cur_pos[j*20+i][1],cur_pos[j*20+i][2],j))
        i=n_part+6
        f.write("%i %15g %15g %15g %15g\n" % (i+j*20,cur_pos[j*20+i][0], cur_pos[j*20+i][1],cur_pos[j*20+i][2],j))
        i=n_part
        f.write("%i %15g %15g %15g %15g\n" % (i+j*20,cur_pos[j*20+i][0], cur_pos[j*20+i][1],cur_pos[j*20+i][2],j))
        f.write("\n\n")
        i=n_part+2
        f.write("%i %15g %15g %15g %15g\n" % (i+j*20,cur_pos[j*20+i][0], cur_pos[j*20+i][1],cur_pos[j*20+i][2],j))
        i=n_part+3
        f.write("%i %15g %15g %15g %15g\n" % (i+j*20,cur_pos[j*20+i][0], cur_pos[j*20+i][1],cur_pos[j*20+i][2],j))
        i=n_part+5
        f.write("%i %15g %15g %15g %15g\n" % (i+j*20,cur_pos[j*20+i][0], cur_pos[j*20+i][1],cur_pos[j*20+i][2],j))
        i=n_part+4
        f.write("%i %15g %15g %15g %15g\n" % (i+j*20,cur_pos[j*20+i][0], cur_pos[j*20+i][1],cur_pos[j*20+i][2],j))
        i=n_part+2
        f.write("%i %15g %15g %15g %15g\n" % (i+j*20,cur_pos[j*20+i][0], cur_pos[j*20+i][1],cur_pos[j*20+i][2],j))
        f.write("\n\n")
        for i in range(n_part+8,n_part+12):
            f.write("%i %15g %15g %15g %15g\n" % (i+j*20,cur_pos[j*20+i][0], cur_pos[j*20+i][1],cur_pos[j*20+i][2],j))
        for i in range(n_part+8,n_part+9):
            f.write("%i %15g %15g %15g %15g\n" % (i+j*20,cur_pos[j*20+i][0], cur_pos[j*20+i][1],cur_pos[j*20+i][2],j))
        f.write("\n\n")
        for i in range(n_part+12,n_part+16):
            f.write("%i %15g %15g %15g %15g\n" % (i+j*20,cur_pos[j*20+i][0], cur_pos[j*20+i][1],cur_pos[j*20+i][2],j))
        for i in range(n_part+12,n_part+13):
            f.write("%i %15g %15g %15g %15g\n" % (i+j*20,cur_pos[j*20+i][0], cur_pos[j*20+i][1],cur_pos[j*20+i][2],j))
        f.write("\n\n")
        for i in range(n_part+16,n_part+20):
            f.write("%i %15g %15g %15g %15g\n" % (i+j*20,cur_pos[j*20+i][0], cur_pos[j*20+i][1],cur_pos[j*20+i][2],j))
        for i in range(n_part+16,n_part+17):
            f.write("%i %15g %15g %15g %15g\n" % (i+j*20,cur_pos[j*20+i][0], cur_pos[j*20+i][1],cur_pos[j*20+i][2],j))
        f.write("\n\n")
    f.close()
    return;
### printdipoles
def printdip( filename, sys, n_part,radius_big_particle):
    cur_pos = sys.part[0:n_part].pos
    dipoles = sys.part[0:n_part].dip
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
        fout.write("part%d={{Blue,Ball[{%.15f, %.15f, %.15f}, %.15f]}" % (i,system.part[i].pos[0],system.part[i].pos[1], system.part[i].pos[2], radius_big_particle ))
        for j in range(8):
            fout.write(", {Orange,Ball[{%.15f, %.15f, %.15f}, %.15f]}" % (system.part[i*44+n_part+j].pos[0],system.part[i*44+n_part+j].pos[1],system.part[i*44+n_part+j].pos[2],radius_coners ))
        for j in range(8,20):
            fout.write(", {White,Ball[{%.15f, %.15f, %.15f}, %.15f]}" % (system.part[i*44+n_part+j].pos[0],system.part[i*44+n_part+j].pos[1],system.part[i*44+n_part+j].pos[2], radius_midpoints ))
        for j in range(20, 44):
            fout.write(", {Magenta,Ball[{%.15f, %.15f, %.15f}, %.15f]}" % (system.part[i * 44 + n_part + j].pos[0], system.part[i * 44 + n_part + j].pos[1],system.part[i * 44 + n_part + j].pos[2], radius_red))
            print("i={0}, j={1}".format(i*44+n_part+j, j))
        fout.write(",{Green, Arrow[{{%.15f, %.15f, %.15f},am*{%.15f, %.15f, %.15f }+{%.15f, %.15f, %.15f }}]}};\n" % (system.part[i].pos[0],system.part[i].pos[1], system.part[i].pos[2],
                                                                                                                      system.part[i].dip[0], system.part[i].dip[1], system.part[i].dip[2],
                                                                                                                      system.part[i].pos[0],system.part[i].pos[1], system.part[i].pos[2]))
    fout.close()
    return



## Implortant parameter

n_part = 20
#H_ext=10.0
#phi=12
phi=19
H_ext=10.00
H_fiel=(0.0,0.0, H_ext)
g_c=(0.0,-2.0,0.0)
#H_fiel=(H_ext/math.sqrt(2.0),0.0, H_ext/math.sqrt(2.0))
#temperature_kT=0.00728
temperature_kT=0.1
time_step = 0.0001
dipole_prefactor=5.0
radius_big_particle=0.5
visc=40

#q=2
radius_midpoints=2**(1/4.0)/6*2*radius_big_particle
radius_coners=3**(1/4.0)/6*2*radius_big_particle
x_midpoints=1.0/2**(1/4.0)/3*2*radius_big_particle
x_coner=1.0/3**(1/4.0)/3*2*radius_big_particle
radius_red=0.197964
x_red=0.279866
y_red=0.134857




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
system.seed = 42
system.cell_system.skin = 4.0


system.time_step = time_step


system.virtual_sites = espressomd.virtual_sites.VirtualSitesRelative(have_velocity=True, have_quaternion=True)
#system.virtual_sites = espressomd.virtual_sites.VirtualSitesRelative(have_velocity=False, have_quaternion=False)
# Add particles to the simulation box at random positions
phi0=math.atan(1.0/math.sqrt(2.))*180/math.pi;
dipole=(1, 1, math.sqrt(2.)*math.tan((phi+phi0)*math.pi/180.));
#dipole=(0.5, 0, 0.5);
#dipole=(0, 1,0 );
dipole=dipole/np.linalg.norm(dipole)
for i in range(n_part):
    #rpos=(np.random.random(1) * 10,0, np.random.random(1) * 10)
    rpos=np.random.random(3)* 1.5*math.sqrt(n_part)
    rpos[1]=0
    if i==0 :
        rpos[0]=0
        rpos[2]=0
    else:
        rpos=np.copy(system.part[i-1].pos)
        rpos[0]=-0.42*np.random.random(1)
        rpos[2]=rpos[2]+2.5*radius_big_particle
    system.part.add(id=i,type=0, pos=rpos, dip=dipole, gamma=visc, gamma_rot=1, rotation=[1,1,1], fix=[1,1,1])
    #system.part.add(id=i,type=0, pos=(6*i, 6*i, 6*i), dip=dipole, gamma=1, gamma_rot=1, rotation=[1,1,1])
for i in range(n_part):
    system.part.add(type=2, pos=np.add(system.part[i].pos,(x_coner,x_coner,x_coner))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=2, pos=np.add(system.part[i].pos,(x_coner,x_coner,-x_coner))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=2, pos=np.add(system.part[i].pos,(x_coner,-x_coner,-x_coner))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=2, pos=np.add(system.part[i].pos,(x_coner,-x_coner,x_coner))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=2, pos=np.add(system.part[i].pos,(-x_coner,-x_coner,-x_coner))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=2, pos=np.add(system.part[i].pos,(-x_coner,-x_coner,x_coner))).vs_auto_relate_to(system.part[i].id) 
    system.part.add(type=2, pos=np.add(system.part[i].pos,(-x_coner,x_coner,x_coner))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=2, pos=np.add(system.part[i].pos,(-x_coner,x_coner,-x_coner))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=3, pos=np.add(system.part[i].pos,(-x_midpoints,0,-x_midpoints))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=3, pos=np.add(system.part[i].pos,(-x_midpoints,0,x_midpoints))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=3, pos=np.add(system.part[i].pos,(x_midpoints,0,x_midpoints))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=3, pos=np.add(system.part[i].pos,(x_midpoints,0,-x_midpoints))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=3, pos=np.add(system.part[i].pos,(0,x_midpoints,x_midpoints))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=3, pos=np.add(system.part[i].pos,(0,x_midpoints,-x_midpoints))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=3, pos=np.add(system.part[i].pos,(0,-x_midpoints,-x_midpoints))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=3, pos=np.add(system.part[i].pos,(0,-x_midpoints,x_midpoints))).vs_auto_relate_to(system.part[i].id)        
    system.part.add(type=3, pos=np.add(system.part[i].pos,(-x_midpoints,x_midpoints,0))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=3, pos=np.add(system.part[i].pos,(x_midpoints,x_midpoints,0))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=3, pos=np.add(system.part[i].pos,(x_midpoints,-x_midpoints,0))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=3, pos=np.add(system.part[i].pos,(-x_midpoints,-x_midpoints,0))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=1, pos=np.add(system.part[i].pos, (x_red, x_red, y_red))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=1, pos=np.add(system.part[i].pos, (x_red, x_red, -y_red))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=1, pos=np.add(system.part[i].pos, (x_red, -x_red, y_red))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=1, pos=np.add(system.part[i].pos, (-x_red, x_red, y_red))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=1, pos=np.add(system.part[i].pos, (-x_red, -x_red, y_red))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=1, pos=np.add(system.part[i].pos, (x_red, -x_red, -y_red))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=1, pos=np.add(system.part[i].pos, (-x_red, x_red, -y_red))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=1, pos=np.add(system.part[i].pos, (-x_red, -x_red, -y_red))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=1, pos=np.add(system.part[i].pos, (x_red, y_red, x_red))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=1, pos=np.add(system.part[i].pos, (x_red, y_red, -x_red))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=1, pos=np.add(system.part[i].pos, (x_red, -y_red, x_red))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=1, pos=np.add(system.part[i].pos, (-x_red, y_red, x_red))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=1, pos=np.add(system.part[i].pos, (-x_red, -y_red, x_red))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=1, pos=np.add(system.part[i].pos, (x_red, -y_red, -x_red))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=1, pos=np.add(system.part[i].pos, (-x_red, y_red, -x_red))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=1, pos=np.add(system.part[i].pos, (-x_red, -y_red, -x_red))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=1, pos=np.add(system.part[i].pos, (y_red, x_red,x_red))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=1, pos=np.add(system.part[i].pos, (y_red, x_red, -x_red))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=1, pos=np.add(system.part[i].pos, (y_red, -x_red,x_red))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=1, pos=np.add(system.part[i].pos, (-y_red, x_red,x_red))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=1, pos=np.add(system.part[i].pos, (-y_red, -x_red,x_red))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=1, pos=np.add(system.part[i].pos, (y_red, -x_red, -x_red))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=1, pos=np.add(system.part[i].pos, (-y_red, x_red, -x_red))).vs_auto_relate_to(system.part[i].id)
    system.part.add(type=1, pos=np.add(system.part[i].pos, (-y_red, -x_red, -x_red))).vs_auto_relate_to(system.part[i].id)

#for i in range(n_part):
#    system.part[i].quat=(math.cos(math.pi/3),math.sin(math.pi/3),0,0)
part_pos = system.part[0:n_part].pos
printcords( "cords1.txt", system, n_part)
printdip( "dip1.txt", system, n_part,radius_big_particle)
for i in range(n_part):
    print(system.part[i])

# for i in range(22,46):
#     for j in range(i+1,45):
#         print('i={0} j={1} d={2} '.format(i, j, np.linalg.norm(system.part[i].pos-system.part[j].pos)), )
# input()

#Rotate particles
system.thermostat.set_langevin(kT=100, gamma=visc, gamma_rotation=1.0, act_on_virtual=False)
system.integrator.run(10000)
for i in range(n_part):
  system.part[i].fix=(0,0,0)
system.thermostat.set_langevin(kT=0, gamma=visc, gamma_rotation=1.0, act_on_virtual=False)
system.integrator.run(10)

printcords( "cords1a.txt", system, n_part)
printdip( "dip1a.txt", system, n_part,radius_big_particle)
print_mathematica("cords1a_test.txt")

cutoff_par=2**(1.0/6.0)
lj_eps = 1.0
lj_sig00 = 2.0*radius_big_particle/cutoff_par
lj_sig01 = (radius_big_particle+radius_red)/cutoff_par
lj_sig02 = (radius_big_particle+radius_coners)/cutoff_par
lj_sig03 = (radius_big_particle+radius_midpoints)/cutoff_par
lj_sig11 = 2.0*radius_red/cutoff_par
lj_sig12 = (radius_red+radius_coners)/cutoff_par
lj_sig13 = (radius_red+radius_midpoints)/cutoff_par
lj_sig22 = 2.0*radius_coners/cutoff_par
lj_sig23 = (radius_coners+radius_midpoints)/cutoff_par
lj_sig33 = 2.0*radius_midpoints/cutoff_par
lj_sig40 = radius_big_particle/cutoff_par
lj_sig41 = radius_red/cutoff_par
lj_sig42 = radius_coners/cutoff_par
lj_sig43 = radius_midpoints/cutoff_par
lj_cut00= cutoff_par*lj_sig00
lj_cut01= cutoff_par*lj_sig01
lj_cut02= cutoff_par*lj_sig02
lj_cut03= cutoff_par*lj_sig03
lj_cut11= cutoff_par*lj_sig11
lj_cut12= cutoff_par*lj_sig12
lj_cut13= cutoff_par*lj_sig13
lj_cut22= cutoff_par*lj_sig22
lj_cut23= cutoff_par*lj_sig23
lj_cut33= cutoff_par*lj_sig33
lj_cut40= cutoff_par*lj_sig40
lj_cut41= cutoff_par*lj_sig41
lj_cut42= cutoff_par*lj_sig42
lj_cut43= cutoff_par*lj_sig43
lj_cap = 10
system.non_bonded_inter[0, 0].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig00,cutoff=lj_cut00, shift='auto')
system.non_bonded_inter[0, 1].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig01,cutoff=lj_cut01, shift='auto')
system.non_bonded_inter[0, 2].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig02,cutoff=lj_cut02, shift='auto')
system.non_bonded_inter[0, 3].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig03,cutoff=lj_cut03, shift='auto')
system.non_bonded_inter[4, 0].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig40,cutoff=lj_cut40, shift='auto')
system.non_bonded_inter[4, 1].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig41,cutoff=lj_cut41, shift='auto')
system.non_bonded_inter[4, 2].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig42,cutoff=lj_cut42, shift='auto')
system.non_bonded_inter[4, 3].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig43,cutoff=lj_cut43, shift='auto')
system.non_bonded_inter[1, 1].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig11,cutoff=lj_cut11, shift='auto')
system.non_bonded_inter[1, 2].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig12,cutoff=lj_cut12, shift='auto')
system.non_bonded_inter[1, 3].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig13,cutoff=lj_cut13, shift='auto')
system.non_bonded_inter[2, 2].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig22,cutoff=lj_cut22, shift='auto')
system.non_bonded_inter[3, 3].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig33,cutoff=lj_cut33, shift='auto')
system.non_bonded_inter[2, 3].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig23,cutoff=lj_cut23, shift='auto')
system.force_cap=lj_cap

system.thermostat.turn_off()
# ### Warmup

# In many cases, including this tutorial, particles are initially placed randomly in the simulation box. It is therefore possible that particles overlap, resulting in a huge repulsive force between them. In this case, integrating the equations of motion would not be numerically stable. Hence, it is necessary to remove this overlap. This is done by limiting the maximum force between two particles, integrating the equations of motion, and increasing the force limit step by step as follows:

# In[11]:


warm_steps  = 1000
warm_n_time = 20
min_dist    = 5

system.minimize_energy.init(f_max=20.0, gamma=visc, max_steps=2000, max_displacement=0.01)
system.minimize_energy.minimize()
print("minimize energy")
print(system.analysis.energy()["total"])

printcords( "cords1b.txt", system, n_part)
printdip( "dip1b.txt", system, n_part,radius_big_particle)
system.thermostat.set_langevin(kT=temperature_kT, gamma=visc, gamma_rotation=1.0, act_on_virtual=False)
print(system.analysis.energy()["total"])
print("Warm faze\n")
i = 0
nonbonded_energy = system.analysis.energy()["non_bonded"]
system.force_cap=1
while i < warm_n_time:
    print('{0} {1} {2} {3} {4} {5}'.format(i,np.linalg.norm(system.part[1].pos-system.part[0].pos),system.analysis.energy()["total"], system.analysis.energy()["dipolar"], system.analysis.energy()["kinetic"],system.analysis.energy()["non_bonded"]))
    system.integrator.run(warm_steps)
    system.galilei.kill_particle_motion(rotation =1)
    nonbonded_energy = system.analysis.energy()["non_bonded"]
    print('{0} {1} {2} {3} {4} {5}'.format(i,np.linalg.norm(system.part[1].pos-system.part[0].pos),system.analysis.energy()["total"], system.analysis.energy()["dipolar"], system.analysis.energy()["kinetic"],system.analysis.energy()["non_bonded"]))
    #print(system.analysis.energy())
    i+=1
    lj_cap += 1.0
    system.force_cap=lj_cap

print(system.analysis.energy())

system.thermostat.set_langevin(kT=10*temperature_kT, gamma=visc, gamma_rotation=1.0, act_on_virtual=False)
dds_cpu = DipolarDirectSumCpu(prefactor=dipole_prefactor)
system.actors.add(dds_cpu)
system.integrator.run(steps=0, recalc_forces=True)

printcords( "cords.txt", system, n_part)
printdip( "dip.txt", system, n_part,radius_big_particle)

H_constraint = espressomd.constraints.HomogeneousMagneticField(H=H_fiel)
gravity = espressomd.constraints.Gravity(g=g_c)
system.constraints.add(H_constraint)
system.constraints.add(gravity)
system.constraints.add(shape=Wall(
    dist=-4, normal=[0, 1, 0]),  particle_type=4)

# ### Integrating equations of motion and taking measurements

# Once warmup is done, the force capping is switched off by setting it to zero.

# In[12]:
local_temperature_kT=10*temperature_kT;
system.thermostat.set_langevin(kT=local_temperature_kT, gamma=visc, gamma_rotation=1.0, act_on_virtual=False)
print('{0} {1} {2} {3} {4}'.format(i, system.analysis.energy()["dipolar"], system.analysis.energy()["kinetic"],np.linalg.norm(system.part[1].pos-system.part[0].pos),np.dot(system.part[1].dip,system.part[0].dip)))
act_min_dist = system.analysis.energy()["dipolar"]
#while act_min_dist > -20 :
print(system.analysis.energy())
print('{0} {1} {2} {3} {4}'.format(i, system.analysis.energy()["dipolar"], system.analysis.energy()["kinetic"],
                                   np.linalg.norm(system.part[1].pos - system.part[0].pos),
                                   np.dot(system.part[1].dip, system.part[0].dip)))
#input()
visualizer = visualization.openGLLive(system)
def main_thread():
    for i in range(2000):
        print('{0} {1} {2} {3} {4}'.format(i, system.analysis.energy()["dipolar"], system.analysis.energy()["kinetic"],np.linalg.norm(system.part[1].pos-system.part[0].pos),np.dot(system.part[1].dip,system.part[0].dip)))
        system.integrator.run(3000)
        visualizer.update()
        #print(system.analysis.energy()["dipolar"])
        print('{0} {1} {2} {3} {4}'.format(i, system.analysis.energy()["dipolar"], system.analysis.energy()["kinetic"],np.linalg.norm(system.part[1].pos-system.part[0].pos),np.dot(system.part[1].dip,system.part[0].dip)))
        print('{0} {1} {2} {3} {4} {5} {6}'.format(i, np.dot(H_fiel,system.part[0].dip), np.dot(H_fiel,system.part[1].dip), (system.part[1].pos[2]-system.part[0].pos[2]), (system.part[1].pos[0]-system.part[0].pos[0]), (system.part[1].pos[1]-system.part[0].pos[1]), np.dot(system.part[1].dip,system.part[0].dip)))
        print('{0} {1} {2}'.format(i, np.dot(H_fiel, system.part[0].dip), system.part[0].pos[1]))
#,  (system.part[1].pos[2]-system.part[0].pos[2]), np.linalg.norm(system.part[1].pos[0]-system.part[0].pos[0]) ,np.dot(system.part[1].dip,system.part[0].dip)))  
        string='data/cords_2d_'
        string +=`i`
        string +='.txt'
        printcords( string, system, n_part)
        string='data/dip_2d_'
        string +=`i`
        string +='.txt'
        printdip( string, system, n_part,radius_big_particle)
        global local_temperature_kT 
        local_temperature_kT/=1.005
        system.thermostat.set_langevin(kT=local_temperature_kT, gamma=visc, gamma_rotation=1.0, act_on_virtual=False)
        #print('{0} {1} {2} {3}'.format(system.analysis.energy()["dipolar"], system.analysis.energy()["kinetic"]),system.part[0].dip,system.part[1].dip)    
    with open('final_pos.txt', 'wb') as f:
        for i in range(21*n_part):
            np.savetxt(f, np.transpose(system.part[i].pos), fmt='%.15e')
    print_mathematica("final_pos_mat.txt")
t = Thread(target=main_thread)
t.daemon = True
t.start()
visualizer.start()

print(np.linalg.norm(system.part[1].pos-system.part[0].pos))
for i in range(n_part):
    print(system.part[i].pos)
print(system.analysis.energy())

system.force_cap=0


for i in range(n_part+1):
    print(system.part[i])
print(system.analysis.energy())
printcords( "cords2.txt", system, n_part)
printdip( "dip2.txt", system, n_part,radius_big_particle)
system.minimize_energy.init(f_max=20000.0, gamma=0.01, max_steps=2000, max_displacement=0.0005)
system.minimize_energy.minimize()
print("minimize energy")
print(system.analysis.energy()["dipolar"])
