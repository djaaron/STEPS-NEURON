import fractions
import sys
# import Meshpy
from meshpy.tet import MeshInfo, build
# import STEPS
import steps.utilities.meshio as smeshio
import steps.geom as sgeom
import steps
# import NEURON
from neuron import h, gui
h.load_file('stdrun.hoc')


###=========================================== ABOUT ===========================================###
# DESCRIPTION:
#	Reads section 3D info from a given .hoc (NEURON) file, then
#	generates a coarse cube-voxel mesh of user-defined granularity with Meshpy/Tetgen
#	saves that mesh in STEPS formats (xml & ASCII) with user-inputed file name, loads
#	the format into STEPS and associates all tets in the STEPS mesh with the .hoc sections
#	which exist inside those tets by a tet_hoc dictionary; tet_hoc[tet_index] = [encapsulated .hoc sections]
#	and concludes by returning this dictionary. 
# AUTHOR:
#	Daniel Aaron, Blue Brain Project, 25 July 2016. 
###=========================================== ABOUT ===========================================###

def coarse_gen(cd,bdim,bmin,file_name):
	# INPUT: voxel cube side length (fl), box dimensions (fl), box minimum values (fl), name of file to create for mesh output (str)
	# This function generates a coarse mesh given info from NEURON sections.
	# Returns this mesh as a STEPS tuple, of which entry 0 is a steps.geom.Tetmesh object,
	# also saves STEPS format files for the Tetmesh to current directory under given file name (xml & ASCII file formats).

	## CALCULATE BOX INFO ##
	box_w = bdim[0]
	box_l = bdim[1]
	box_h = bdim[2]

	min_x = bmin[0]
	min_y = bmin[1]
	min_z = bmin[2]

	cl = cd

	if not (box_w % cl == box_l % cl == box_h % cl == 0):
		print 'ERRROR: voxel cube side length not a common factor of all box dimensions'
		sys.exit()
	
	wpoints = 1+(box_w)/cl
	lpoints = 1+(box_l)/cl 
	hpoints = 1+(box_h)/cl

	print "cube side length: ", cl
	print "box w,l,h: ", box_w,box_l,box_h
	print "w,l,h # points: ", wpoints,lpoints,hpoints

	cpoints = []
	hfacets = []
	vfacets = []

	## GENERATE POINTS AND FACETS ##
	for hp in range(int(hpoints)):
		for lp in range(int(lpoints)):
			for wp in range(int(wpoints)):
				cpoints.append((min_x+wp*cl,min_y+lp*cl,min_z+hp*cl))
				pindex = (hp*lpoints*wpoints)+(lp*wpoints)+wp
				# horizontal facets
				if (wp < int(wpoints)-1 and lp < int(lpoints)-1):
					hfacets.append([int(pindex), int(pindex+1), int(pindex+1+wpoints), int(pindex+wpoints)])
				# vertical facets
				if (hp > 0):
					if (wp > 0):
						vfacets.append([int(pindex),int(pindex-1),int(pindex-1-lpoints*wpoints),int(pindex-lpoints*wpoints)])
					if (lp > 0):
						vfacets.append([int(pindex),int(pindex-wpoints),int(pindex-wpoints-lpoints*wpoints),int(pindex-lpoints*wpoints)])
	all_facets = hfacets+vfacets

	## PASS MESH TO STEPS ##
	mesh_info = MeshInfo()
	mesh_info.set_points(cpoints)
	mesh_info.set_facets(all_facets)
	m = build(mesh_info)

	# write mesh proxies
	nodeproxy = smeshio.ElementProxy('node',3)
	for i, p in enumerate(m.points):
	    nodeproxy.insert(i,p)

	tetproxy = smeshio.ElementProxy('tet',4)
	newtet = [0,0,0,0]
	for i, t in enumerate(m.elements):
	    newtet[0] = nodeproxy.getSTEPSID(int(t[0]))
	    newtet[1] = nodeproxy.getSTEPSID(int(t[1]))
	    newtet[2] = nodeproxy.getSTEPSID(int(t[2]))
	    newtet[3] = nodeproxy.getSTEPSID(int(t[3]))
	    tetproxy.insert(i, newtet)
	 
	# build mesh from proxies and save in STEPS format (xml & ASCII)
	nodedata = nodeproxy.getAllData()
	tetdata = tetproxy.getAllData()
	newmesh = steps.geom.Tetmesh(nodedata, tetdata)
	smeshio.saveMesh(file_name,newmesh)

	# Load mesh into STEPS
	steps_mesh = smeshio.loadMesh(file_name)
	print "STEPS loaded Tetmesh successfully."
	print "# Vertices: ", steps_mesh[0].countVertices()
	print "# Tets: ", steps_mesh[0].countTets()
	print "# Faces/Tris: ", steps_mesh[0].countTris()
	return steps_mesh

def tet_associate(sm):
	# INPUT: Tetmesh object
	# This function associates .hoc sections to the tets which contain them, and returns the association dictionary. 
	# It's packaged separately from the coarse() function so it can be used independently in scripts which load a Tetmesh from file
	# and don't need any Tetmesh to be generated from a .hoc file. 
	
	tet_hoc = {}
	for s in h.allsec():
		for j in range(int(h.n3d())):
			containing_tet = sm.findTetByPoint((h.x3d(j),h.y3d(j),h.z3d(j)))
			if (containing_tet) not in tet_hoc.keys():
				tet_hoc[containing_tet] = [s]
			elif (containing_tet) in tet_hoc.keys():
				if (s) not in tet_hoc[containing_tet]:
					tet_hoc[containing_tet].append(s)

	return tet_hoc

def coarse(hoc_filename,cube_length,save_filename):
	# INPUT: NEURON .hoc filename to import (str), voxel cube side length (fl/str for gcd), name of file to create for mesh output (str)
	# 		>> cube_length: 'gcd' (gcd of box dims) OR floating-point (must be common factor of all box dims)
	# This function reads in NEURON data and passes their info to
	# coarse_gen(), then associates the tets of the STEPS Tetmesh object returned
	# to the NEURON sections which exist inside of them. 
	# Returns a tet_hoc dictionary -- tet_hoc[tet_index] = [encapsulated hoc section references] -- as well as the Tetmesh object

	## GET HOC SECTION INFO ##
	h.load_file(hoc_filename)

	allp = [[],[],[]]
	for s in h.allsec():
		for j in range(int(h.n3d())):
			allp[0].append(h.x3d(j))
			allp[1].append(h.y3d(j))
			allp[2].append(h.z3d(j))

	maxl = [max(allp[0]),max(allp[1]),max(allp[2])]
	minl = [min(allp[0]),min(allp[1]),min(allp[2])]
	bdim = [ maxl[0] - minl[0], maxl[1] - minl[1], maxl[2] - minl[2] ]
	print "dims: ", bdim
	print "mins: ", minl

	## CREATE COARSE MESH ##
	if (cube_length == 'gcd'):
		gcd = fractions.gcd(fractions.gcd(bdim[0],bdim[1]),fractions.gcd(bdim[2],bdim[1]))
		print "GCD: ", gcd
		cube_length = gcd

	sm = coarse_gen(cube_length,bdim,minl,save_filename)

	## ASSOCIATE HOC SECTIONS WITH THEIR TETS ##
	tet_hoc = tet_associate(sm[0])
	return tet_hoc, sm[0]





