#!/usr/bin/env python
# coding: utf-8

#read water configurations in NX4A format,
#calculate the pairwise products of local order parameter for ice 7,
#and output them with pair distances.

#2012-10-19 tested with ice VII crystal structure.

import numpy
import numpy.linalg
import math
import sys



############################################################################## Quaternion
def quat2rotmat(q):
    a,b,c,d = q
    sp11=(a*a+b*b-(c*c+d*d))
    sp12=-2.0*(a*d+b*c)
    sp13=2.0*(b*d-a*c)
    sp21=2.0*(a*d-b*c)
    sp22=a*a+c*c-(b*b+d*d)
    sp23=-2.0*(a*b+c*d)
    sp31=2.0*(a*c+b*d)
    sp32=2.0*(a*b-c*d)
    sp33=a*a+d*d-(b*b+c*c)
    return numpy.matrix([[sp11,sp12,sp13], [sp21,sp22,sp23], [sp31,sp32,sp33]])

############################################################################## Clustering
#The cluster connectivity is saved in the "group dictionary" as a tree structure.
#The value of the group dictionary, group[i], indicates:
#group[i]>=0 ==> group[i] is the parent node for node i
#group[i]<0  ==> i is the root node; -group[i] is number of nodes in the same tree.

#x: node label
#group[]: group dictionary
#return value: root node of the group.
def MyGroup(x,group):
    while group[x] >= 0:
        x = group[x]
    return x


#register the nodes x and y as connected
def BindNodes(x,y,group):
    xg = MyGroup(x,group)
    yg = MyGroup(y,group)
    if xg != yg:
		#merge the trees
        group[yg] += group[xg]
        group[xg] =  yg


##################################################### Adjacency by distance and direction
#determine whether vector points to the direction
def IsAdjacent( vec, direc, thres ):
    siz = numpy.linalg.norm(vec)
    if siz < thres:
        #distance is short enough
	ip  = numpy.dot(vec, direc)
	if ip < -siz*0.95 or siz*0.95 < ip:
            #the vector directs
            return True
    return False


############################################################### Periodic boundary utility
#relative position vector at the periodic boundary condition
def Wrap( vector, box ):
	for dim in range(len(vector)):
		vector[dim] -= math.floor( vector[dim] / box[dim] + 0.5 ) * box[dim]
	return vector


################################################################################## Loader
#load a specific file type of molecular coordinate
def LoadNX4A(file):
	h1v = numpy.array((0., +7.56950327263661182e-01,  5.20784245882928820e-01))
	h2v = numpy.array((0., -7.56950327263661182e-01,  5.20784245882928820e-01))
	line = file.readline()
	nmol = int(line)
	waters = []
	for i in range(nmol):
		line = file.readline()
		columns = line.split()
		x,y,z,a,b,c,d = map(float, columns)
		mat = quat2rotmat([a,b,c,d])
		h1 = numpy.dot(mat, h1v)
		h2 = numpy.dot(mat, h2v)
		com = numpy.array([x,y,z])
		waters.append((com,h1,h2))
	return waters

##########################################################################Order parameter
#1 in ice VII, smaller in plastics
#IT IS COMPLETELY WRONG!!!
def OrderParameter( site, neighbors, latticeGroup ):
	myg = MyGroup( site, latticeGroup )
	op = 0.0
	for n in neighbors:
		g = MyGroup( n, latticeGroup )
		if g == myg:
			#increase if the neighbor is on the same sublattice.
			op += 1
		else:
			#decrease otherwise
			op -= 1
	#normalize so that op of ice 7 becomes 1.
	return op / 4.0



def OrderParameter2( site, neighbors, waters, box ):
    pos = waters[site][0]
    op = 0.0
    for n in neighbors:
        npos = waters[n][0]
        dir = Wrap(npos - pos, box)
        parity = 1
        for dim in range(3):
            if dir[dim] < 0:
                parity *= -1
        op += parity
    return op / 4.0


#return value: array of 7-element lists.
#7 elements = 3 for coordinates, 4 for quaternions.
def Configure(file):
	while True:
		line = file.readline()
		if len(line) == 0:
			break
		columns = line.split()
		if len(columns) > 0:
			if columns[0] == "@BOX3":
				line = file.readline()
				box = numpy.array(map(float,line.split()))
			elif columns[0] =="@NX4A":
				waters = LoadNX4A(file)
	return waters,box


def Configure_gro(file):
    while True:
        line = file.readline()
        if len(line) == 0:
            break
        if line.find("frame") >= 0:
            line = file.readline()
            nsite = int(line)
            waters = []
            water = []
            for i in range(nsite):
                line = file.readline()
                molnum = line[0:5]
                molname = line[5:10]
                elemname = line[10:15]
                sitenum = line[15:20]
                x = float(line[20:28])*10
                y = float(line[28:36])*10
                z = float(line[36:44])*10
                xyz = numpy.array([x,y,z])
                if elemname == "   OW":
                    water.append(xyz)
                elif elemname == "  HW1":
                    water.append(xyz)
                elif elemname == "  HW2":
                    water.append(xyz)
                    waters.append(water)
                    water = []
            line = file.readline()
            box = line.split()
            box[0] = float(box[0])*10
            box[1] = float(box[1])*10
            box[2] = float(box[2])*10
            box = numpy.array(box)
            break
    #convert hydrogen positions from absolute to relative
    for water in waters:
        com = water[0]
        h1  = water[1]
        h2  = water[2]
        water[1] = Wrap(h1-com,box)
        water[2] = Wrap(h2-com,box)
    return waters, box
            



def main(debug=False):
	#crystal axes
	ex = numpy.array([1.0,0.0,0.0])
	ey = numpy.array([0.0,1.0,0.0])
	ez = numpy.array([0.0,0.0,1.0])

	#stage 1: load the structure
	waters,box = Configure_gro(sys.stdin)
        if len(waters) ==0:
            sys.exit(0)
	if debug:
		print waters
	#stage 2: determine the double lattices
	thres = 3.7 #Å; threshold for O-O adjacency on the same sublattice
	latticeGroup = dict()
	for i in range(len(waters)):
		#-1 means an isolated node.
		latticeGroup[i] = -1
	for i in range(len(waters)):
		for j in range(i+1, len(waters)):
			delta = Wrap(waters[i][0]-waters[j][0],box)
			#find the nearest neighbor molecule along the x,y, and z axes.
			if IsAdjacent(delta, ex, thres):
				BindNodes(i,j,latticeGroup)
			if IsAdjacent(delta, ey, thres):
				BindNodes(i,j,latticeGroup)
			if IsAdjacent(delta, ez, thres):
				BindNodes(i,j,latticeGroup)
	if debug:
		for i in range(len(waters)):
			print i,MyGroup(i, latticeGroup)
	#check the number of groups (should be 2)
	#but you may have to loosen the condition.
	groups = set()
	for i in range(len(waters)):
		groups.add(MyGroup(i, latticeGroup))
	if len(groups) != 2:
		sizes = [-latticeGroup[i] for i in groups]
		print "Illegal number of groups: ", sizes
		sys.exit(1)
        else:
            sys.stderr.write("No problem. Go ahead!\n")
	#stage 3: determine the HB networks
	hb = dict()
	for i in range(len(waters)):
		for j in range(i+1, len(waters)):
			delta = Wrap(waters[i][0]-waters[j][0],box)
			if numpy.linalg.norm(delta) < 3.2:#Å; threshold for O-O adjacency
				#molecule pair is close enough
				#find the shortest OH pair distance among four possible combinations
				d1 = numpy.linalg.norm( delta + waters[i][1] )
				d2 = numpy.linalg.norm( delta + waters[i][2] )
				d3 = numpy.linalg.norm( delta - waters[j][1] )
				d4 = numpy.linalg.norm( delta - waters[j][2] )
				if min(d1,d2,d3,d4) < 2.2:#Å; threshold for O-H adjacency
					#register as a hydrogen bond.
					#bond direction is ignored.
					if not hb.has_key(i):
						hb[i] = []
					if not hb.has_key(j):
						hb[j] = []
					hb[i].append(j)  #hb[i] contains the labels of the HB partners
					hb[j].append(i)
        if debug:
            print hb[0]
	#stage 4: determine the order parameters and output it with distance
        op = []
	for i in range(len(waters)):
            op.append(OrderParameter2(i,hb[i],waters,box))
            if debug:
                print i, op[i]
	for i in range(len(waters)):
		gi = MyGroup(i,latticeGroup)
		opi = op[i]
		for j in range(i, len(waters)):
			gj = MyGroup(j,latticeGroup)
			opj = op[j]
			if gi == gj:
				#i and j are on the same sublattice
				delta = Wrap(waters[i][0]-waters[j][0],box)
				dist  = numpy.linalg.norm( delta )
				if dist < box[0]/2:
                                    print i,j, dist, opi*opj
while True:	
    main()
