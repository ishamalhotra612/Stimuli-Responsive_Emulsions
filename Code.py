import sys
from hoomd import *
from hoomd import md
from hoomd import group
from hoomd import deprecated
from pprint import pprint
import freud
import random
import numpy 

context.initialize('--mode=gpu')


import math
sigma=1   # we assume this as our unit

def various_calc(step, file2):
    file2.write(str(wallForce2.get_net_force(group.all())))
    file2.write('\t')
    file2.write(str(wallForce3.get_net_force(group.all())))
    file2.write("\n")
    file2.flush()
    


NAME='Initial_Configuration.gsd'
system= init.read_gsd(filename=NAME)
nl = md.nlist.cell()
all = group.all()

########
#####   INTERFACE
########

zcm=0
xcm=0
ycm=0
n=0
## center of mass zcm
for p in system.particles:
    zcm = zcm + p.position[2]
    xcm = xcm + p.position[0]
    ycm = ycm + p.position[1]
    n=n+1   
zcm=zcm/n
xcm=xcm/n
ycm=ycm/n
print(zcm)

########
#### WALL
########
#float(sys.argv[4]) = d that defines the distance between the two interfaces.
mywall=md.wall.group(md.wall.plane(origin=(0,0,zcm+float(sys.argv[4])),normal=(0,0,-1))) 
mywall3=md.wall.group(md.wall.plane(origin=(0,0,zcm-float(sys.argv[4])),normal=(0,0,1)))
my_rcut=float(sys.argv[3])*math.pow(2,1./6.)+sigma*float(sys.argv[5])
wallForce2=md.wall.lj(mywall, r_cut=my_rcut)
wallForce2.force_coeff.set(['A','B','C','M','N','K','P'],epsilon=float(sys.argv[2]),sigma=float(sys.argv[3]), r_extrap=float(sys.argv[3])*1.11)
wallForce3=md.wall.lj(----mywall3, r_cut=my_rcut)
wallForce3.force_coeff.set(['A','B','C','M','N','K','P'],epsilon=float(sys.argv[2]),sigma=float(sys.argv[3]), r_extrap=float(sys.argv[3])*1.11)

sigmacylinder=20.*sigma
# float(sys.argv[11]) gives the value of R_c defined in the paper
mycylinder=md.wall.group(md.wall.cylinder(r=float(sys.argv[11]), origin=(xcm,ycm,0), axis=(0,0,1),inside=True))
my_rcut_cylinder=sigmacylinder
cylinderForce=md.wall.force_shifted_lj(mycylinder,r_cut=my_rcut_cylinder)
cylinderForce.force_coeff.set(['A','B','C','M','N','K','P'],epsilon=5.0,sigma=sigmacylinder)


########
######  THERMORESPONSIVE FIELD
########

def longT(r,rmin,rmax,epsilon,alpha,delta,beta,sigma):
    # function defined by Rovigatti and co.
    V=0.5*alpha*epsilon*( math.cos(delta*(r/sigma)*(r/sigma)+beta) -1 )
    F=alpha*epsilon* math.sin(delta*(r/sigma)*(r/sigma)+beta) *delta*r/(sigma*sigma) 
    return (V,F);

# short range
sigma_respons=1.
mycut_s=math.pow(2.,1./6.)*sigma
#long range
R0=1.5
mycut_l=R0*sigma
mydelta=math.pi/( 2.25-math.pow(2.,1./3.) )
mybeta=2.*math.pi -2.25*mydelta
thermoPot_l = md.pair.table(width=1000,nlist=nl)
#alpha=float(sys.argv[9]) is an input parameter, alpha=0 for swollen state and alpha=1 for collapsed state
thermoPot_l.pair_coeff.set(['P','K','M','C'],['P','K','M','C'],func=longT,rmin=mycut_s,rmax=mycut_l,coeff=dict(epsilon=1.,alpha=float(sys.argv[9]),delta=mydelta,beta=mybeta,sigma=sigma_respons))

########
#### output
########

NAME='energy_EQ.log'
dump_run=1e3
log1 = analyze.log(filename=NAME, quantities=['time', 'potential_energy', 'temperature','pair_lj_energy', 'volume', 'kinetic_energy', 'pressure'], period=dump_run, overwrite=True)


########
#### WCA repulsion
########

my_rcut=math.pow(2,1./6.)    
potRep = md.pair.lj(r_cut=my_rcut,nlist=nl)
potRep.pair_coeff.set(['M','C','P','K'],['M','C','P','K'],epsilon=1.,sigma=1.)
potRep.set_params(mode="shift")

########
#### define bond potential
########

CCfene=md.bond.fene()
CCfene.bond_coeff.set('CoreCore', k=7.5, r0=1.5, sigma=1.,epsilon=1.)


#########
####   actual RUN
#########

totrun=float(sys.argv[6])
del_t=0.00009
md.integrate.mode_standard(dt=del_t)
b=md.integrate.brownian(group=all,kT=1.0,seed=int(sys.argv[1]))


NAME='Trajectory.gsd'
dump.gsd(filename=NAME, group=group.all(), period=1e5, overwrite=True);


NAME3='forces.txt'
file2 = open(NAME3, "w+")
run(int(sys.argv[6]), callback_period=1000, callback= lambda step: various_calc(step, file2), quiet=True )
file2.close()


#########
#### output the snapshot
#########
NAME='Final_Configuration.gsd'
dump.gsd(filename=NAME, group=group.all(), truncate=True, period=None, phase=0,static=['attribute', 'property' , 'momentum', 'topology'])
