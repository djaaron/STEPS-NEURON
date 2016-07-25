import mesh_coarse_diffdemo
import math
import numpy
import pylab
import random
import time

import steps.model as smodel
import steps.solver as solvmod
import steps.geom as stetmesh
import steps.rng as srng

import steps.utilities.meshio as smeshio

from neuron import h, gui
hoc_file = "coarse_diffdemo_cells.hoc"
h.load_file('stdrun.hoc')
h.load_file(hoc_file)


# simulation run control
tstop = 30
dt = 0.25 
INT = tstop
DT = dt
h("tstop = {0}".format(tstop))
h("dt = {0}".format(dt))
h("steps_per_ms = {0}".format(1.0/dt))
NINJECT = 1.0e5

# The diffusion constant for our diffusing species (m^2/s)
DCST= 60.0 #20.0e-12

beg_time = time.time()

########################################################################

def printtime(end_time):
    
	totsecs = int(end_time-beg_time)
	sec = totsecs%60
	totmin = totsecs/60
	min = totmin%60
	hours = totmin/60
	
	print 'Simulation time: %d h, %d min, %d sec' %(hours, min, sec)

########################################################################

def gen_model():         
    mdl = smodel.Model()        
    A = smodel.Spec('A', mdl)                        
    vsys = smodel.Volsys('cytosolv', mdl)
    diff_A = smodel.Diff('diff_A', vsys, A)
    diff_A.setDcst(DCST)
    return mdl

########################################################################

def gen_geom():
    # Generate STEPS mesh and associate tets
    coarse_results = mesh_coarse_diffdemo.coarse(hoc_file,"gcd","novel_name")
    t_h = coarse_results[0]
    mesh = coarse_results[1]
    # Find the total number of tetrahedrons in the mesh	
    ntets = mesh.countTets()
    print "Mesh has # Tets = ", ntets
    # Create a compartment containing all tetrahedron
    comp = stetmesh.TmComp('cyto', mesh, range(ntets))
    comp.addVolsys('cytosolv')
    return t_h, mesh

########################################################################

geom_results = gen_geom()
t_h = geom_results[0]
tmgeom = geom_results[1]
model = gen_model()

rng = srng.create('mt19937', 512) 
rng.initialize(2903)

tpnts = numpy.arange(0.0, INT, DT)
# count "time points"
ntpnts = tpnts.shape[0]

sim = solvmod.Tetexact(model, tmgeom, rng)

## RUN SIMULATION ##
sim.reset()
h('finitialize()')
print "NEURON tstop: ", h.tstop, " and dt: ", h.dt
print "STEPS tstop: ", tstop, " and dt: ", dt
sim.setTetCount(39, 'A', NINJECT)
for j in range(ntpnts):
	# Advance STEPS in time
    sim.run(tpnts[j])
    # Update NEURON tet sections with new concentrations
    for tet in t_h.keys():
    	tetconc = sim.getTetConc(int(tet), 'A')
    	print "tet: ", tet, " conc[A]: ", tetconc
    	for s in t_h[tet]:
    		s(0.5).molA_cdiffdemo = tetconc
    		s(0.5).rand_cdiffdemo = random.uniform(0.8,1.2)
    # Advance NEURON in time
    if (j > 0.0):
    	h("fadvance()")
    	# Update STEPS tets with results from NEURON sections
    	for tet in t_h.keys():
	    	for s in t_h[tet]:
	    		sim.setTetCount(tet, 'A', s(0.5).Ai)

printtime(time.time())
print "NEURON t: ", h.t, " STEPS t: ", tpnts[j]
