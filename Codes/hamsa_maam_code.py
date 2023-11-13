import time
start_time = time.time()
#USE PYTHON 3.* and Networkx 2.3
import math
import numpy as np
import mdtraj as md
from mdtraj.formats import XTCTrajectoryFile
import networkx as nx
import matplotlib.pyplot as plt

print("importing done")

#input filenames
xtcname = 'water_traj/bulk_tip3p_all_part1.xtc'
groname = 'water_traj/bulk_tip3p_all.gro'

print("filenames read")

#output filenames
extnname='.dat'
file1name='hamsa_maam_outputs/Network_Stats_'+extnname
file2name='hamsa_maam_outputs/molhb_histogram_'+extnname
file3name='hamsa_maam_outputs/hbweight_histogram_'+extnname
file5name='hamsa_maam_outputs/Pathlength_histogram_'+extnname
file6name='hamsa_maam_outputs/PerFrame_HistCount_'+extnname
file7name='hamsa_maam_outputs/Centrality_histogram_'+extnname
#
traj = md.load(xtcname,top=groname)   
info = traj.xyz.shape


n_frame = info[0]
n_atom = info[1]
n_coord = info[2]
coord = traj.xyz
topol = traj.topology
print (topol)
nres =  topol.n_residues

traj_file = XTCTrajectoryFile(xtcname)
xyz,tot_time,step,box = traj_file.read()

box_len =[]
for i in range(0,n_frame):
	box_len.append(box[i][0][0])

#print len(box_len)
#Donors
donors = topol.select('name O')
ndon = len(donors)
#Acceptors
acceptors = topol.select('name O')
nacc = len(acceptors)


def AINT(val):
	return int(val)

def anint(value):
	return round(value)

def distance_sq(frame,n1,n2):
	coor1 = coord[frame][n1]
	coor2 = coord[frame][n2]
	boxl = box_len[frame]
	dx = coor2[0] - coor1[0]
	dy = coor2[1] - coor1[1]
	dz = coor2[2] - coor1[2]
	dx = dx - (boxl * anint(dx/(boxl)))
	dy = dy - (boxl * anint(dy/(boxl)))
	dz = dz - (boxl * anint(dz/(boxl)))
	new_dist_sq = (dx)**2 + (dy)**2 + (dz)**2
	return new_dist_sq

def get_pbcacc(frame,n1,n2):
	coor1 = coord[frame][n1]
	coor2 = coord[frame][n2]
	boxl = box_len[frame]
	dx = coor2[0] - coor1[0]
	dy = coor2[1] - coor1[1]
	dz = coor2[2] - coor1[2]
	dx = dx - (boxl * anint(dx/(boxl)))
	dy = dy - (boxl * anint(dy/(boxl)))
	dz = dz - (boxl * anint(dz/(boxl)))
	nxa = coor1[0] + dx
	nya = coor1[1] + dy
	nza = coor1[2] + dz
	newcoor = np.array((nxa,nya,nza))
	return newcoor

def vecsub(A,B):
	xx = A[0] - B[0]
	yy = A[1] - B[1]
	zz = A[2] - B[2]
	sub = np.array((xx,yy,zz))
	return sub
#

def vecadd(A,B):
	xx = A[0] + B[0]
	yy = A[1] + B[1]
	zz = A[2] + B[2]
	add = np.array((xx,yy,zz))
	return add
#

def mag(A):
	mag_A = math.sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2])
	return mag_A

def dot(A,B):
	dot = A[0]*B[0] + A[1]*B[1] + A[2]*B[2]
	return dot

def angle(frame,n1,n2,n3):
	coor1 = coord[frame][n1]
	coor2 = coord[frame][n2]
	coor3 = get_pbcacc(frame,n2,n3)
	vec1 = vecsub(coor1,coor2)
	vec2 = vecsub(coor3,coor2)
	x = (dot(vec1,vec2))/(mag(vec1)*mag(vec2))
	angle_rad = (np.arccos(x))  
	angle = np.rad2deg(angle_rad)      # angle in degrees
	return angle

file1 = open(file1name,"w")
file1.write("#Frame \t #Hbonds \t #HBNetNodes \t #HBNetEdges \n")
file2 = open(file2name,"w")
file2.write("#hbond \t WatAcceptorHist \t WatDonateHist \t WatTotalHist \n")
file3 = open(file3name,"w")
file3.write("#weight \t wwhist\n")
file5 = open(file5name,"w")
file5.write("#pathlength \t LinePath \n")
file6 = open(file6name, "w")
file6.write("#Frame\t #WatHBAcc\t #WatHBDon\t #WatTotHB\t #HBPaths \t #TwoEdges \t #MultiEdges \t #HBNodesNoinpath \n")
file7 = open(file7name, "w")
file7.write("#Frame \t WatClCen \t WatBtCen\n")
hbnet = nx.MultiDiGraph()
watnet = nx.MultiDiGraph() 
wmin = 0.4 # 0.14/0.35 as 3.5 A is the h-bond cut-off distance
wmax = 0.6 # 0.14/0.23 as 2.3 A is the absolute hard sphere radius because gr of wat-wat oxygens is absolutely 0.0 upto 2.3 A
cosmin = np.cos(30*np.pi/180)
cosnorm = 1.0 - cosmin
print(cosmin,cosnorm)
#
nhb = 40
nhb_hist = np.zeros((nhb,3),int)
mpl = 100
path_hist = np.zeros((mpl),int)
wbin = 22
wbin_width = 0.05
weight_hist = np.zeros((wbin),int)
mcen = 100
cen_hist = np.zeros((mcen,2),int)
btbin_width = 0.005
clbin_width = 0.01
#
print("No. of frames: %d"%(n_frame))
for frame in range(0,min(5, n_frame)):
	print("Working on frame " + str(frame))
	count = 0
	for i in range (0,nres):
		rname = topol.residue(i).name
		hbnet.add_node(i,name=rname)
		if(rname == "HOH"):
			watnet.add_node(i)
	#
	for index1 in donors:
		rid1 = topol.atom(index1).residue.index
		for index2 in acceptors:
			rid2 = topol.atom(index2).residue.index
			#if index1 != index2:
			if rid1 != rid2:
				DAdist_sq = distance_sq(frame,index1,index2)
				if DAdist_sq <= 0.1225:
					angle_h1 = angle(frame,index1+1,index1,index2)
					if angle_h1 <= 30:
						count = count + 1
						#frac = 0.14/math.sqrt(DAdist_sq)
						#rratio = (frac - wmin)/(wmax - wmin)
						#cosang = np.cos(angle_h1*np.pi/180)
						#angwt = (cosang - cosmin)/cosnorm
						wt = 1.0
						hbnet.add_edge(rid1,rid2,weight=wt)
					angle_h2 = angle(frame,index1+2,index1,index2)
					if angle_h2 <= 30:
						count = count + 1
						#frac = 0.14/math.sqrt(DAdist_sq)
						#rratio = (frac - wmin)/(wmax - wmin)
						#cosang = np.cos(angle_h2*np.pi/180)
						#angwt = (cosang - cosmin)/cosnorm
						wt = 1.0 
						hbnet.add_edge(rid1,rid2,weight=wt)
					#Angle loops
				#heavy atom distance loop
		#end of acceptor loop
	#end of donor loop
	#in_degree = number of h-bonds accepted by the molecule
	#out_degree = number of h-bonds donated by the molecule
	print("Network created for frame %d"%(frame))
	print("--- %s seconds ---" % (time.time() - start_time))                              
	wwcount = 0
	#weight distribution
	for n, nbrsdict in hbnet.adj.items():
		for nbr,keydict in nbrsdict.items():
			for key,eattr in keydict.items():
				wt = eattr['weight']
				bin = int(wt/wbin_width)
				weight_hist[bin] = weight_hist[bin] + 1
				wwcount = wwcount + 1
				#
	#path length properties
	path = dict(nx.all_pairs_dijkstra_path(hbnet))
	print("after dijkstar path for 1st network")
	print("--- %s seconds ---" % (time.time() - start_time))                              
	#
	waccsum = 0
	wdonsum = 0
	whbsum = 0
	pcount = 0
	mecount=0
	mecount2=0
	file1.write("%d\t\t%d\t\t%d\t\t%d\n"%(frame,count,hbnet.number_of_nodes(), hbnet.number_of_edges()))
	for nd in list(hbnet.nodes):
		acc = hbnet.in_degree(nd)
		don = hbnet.out_degree(nd)
		thb = hbnet.degree(nd)
		waccsum = waccsum + acc
		wdonsum = wdonsum + don
		whbsum  = whbsum + thb
		if(acc < nhb): nhb_hist[acc][0] = nhb_hist[acc][0] + 1
		else: print ("Number of Hbonds exceeds %d"%(nhb))
		if(don < nhb): nhb_hist[don][1] = nhb_hist[don][1] + 1
		else: print ("Number of Hbonds exceeds %d"%(nhb))
		if(thb < nhb): nhb_hist[thb][2] = nhb_hist[thb][2] + 1
		else: print ("Number of Hbonds exceeds %d"%(nhb))
		for nd2 in list(hbnet.nodes):
			if(nd != nd2):
				nedges = hbnet.number_of_edges(nd,nd2)
				if(nedges > 1):
					mecount = mecount + 1
					#print(mecount,nd,nd2,nedges) being directed graph one edge is counted once 
				if(nedges > 2):
					mecount2 = mecount2 + 1
		#
	#hbnet node loop
	print("after loop network for 1st network")
	print("--- %s seconds ---" % (time.time() - start_time))                              
	for ss in path.keys():
		Atar = path[ss]
		for tt in Atar.keys():
			nconn = len(path[ss][tt])-1
			if(nconn > 0):
				pcount = pcount + 1
				if(nconn < mpl): path_hist[nconn] = path_hist[nconn] + 1
				else: print("Total pathlength %d exceeds maximum pathlength of %d"%(nconn,mpl))
				#print frame,ss,tt,nconn,path[ss][tt]
	
	print("after path histogram for 1st network")
	print("--- %s seconds ---" % (time.time() - start_time))                              
	#Centrality Analysis
	#clcen - closeness centrality gives the product of ratio of nodes that can reach a given node (inward paths) and the inverse of the average distance to the given node. closeness distance is based on the incoming distance 
	#btcen - between centrality
	clcen1 = nx.closeness_centrality(hbnet,distance="weight")
	btcen1 = nx.betweenness_centrality(hbnet,weight="weight")
	hnodes = hbnet.nodes()
	nrcount1 = 0
	for nd in hnodes:
		if(hbnet.in_degree(nd) == 0): nrcount1 = nrcount1 + 1
		clbin = int(clcen1[nd]/clbin_width)
		btbin = int(btcen1[nd]/btbin_width)
		if(clbin < mcen):
			cen_hist[clbin][0] = cen_hist[clbin][0] + 1
		else:
			print ("Water %d Closeness %f binid %d is greater than %d"%(nd,clcen1[nd],clbin,mcen))
		if(btbin < mcen):
			cen_hist[btbin][1] = cen_hist[btbin][1] + 1
		else:
			print ("Water %d Betweenness %f binid %d is greater than %d"%(nd,btcen1[nd],btbin,mcen))
	#
	#plotting
	#if(frame==0):
	#	nx.draw(hbnet,nx.circular_layout(hbnet),nodecolor='r', edge_color='b')
	#	plt.show()
	#Clearing networks
	hbnet.clear()
	file6.write("%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\n"%(frame,waccsum,wdonsum,whbsum,pcount,mecount,mecount2,nrcount1))
	#
	print("after centrality histogram for networks")
	print("--- %s seconds ---" % (time.time() - start_time))                              
	#print count,len(nodes),minwt,maxwt,x
#frame loop
#Normalization could use 0,1,2,6,7,8 - #water nodes * nframes; 3,4,5 - #sugar nodes*nframes
for i in range (0,nhb):
	file2.write("%d\t\t%d\t\t%d\t\t%d\n"%(i,nhb_hist[i][0],nhb_hist[i][1],nhb_hist[i][2]))
#
#Normalization 
#Linepaths: if all nodes were connected, no of possible paths would be (n-1)*(n-2) as suggested by betweenness centrality in networkx tutorial
#Linepath normalization: (n-1)*(n-2)*nframes
for i in range (0,mpl):
	file5.write("%d\t\t%d\n"%(i,path_hist[i]))
#
#Normaliztion could use ttcount,twcount,wtcount,wwcount listed in file6
for i in range (0,wbin):
	wt = i*wbin_width + wbin_width/2.0
	file3.write("%f\t\t%d\n"%(wt,weight_hist[i]))
#
#Normalization: 0,1 - #sugar nodes*nframes; 2,3,4,5 - #water nodes*nframes
for i in range (0,mcen):
	cen1 = i*clbin_width + clbin_width/2.0
	cen2 = i*btbin_width + btbin_width/2.0
	file7.write("%f\t\t%d\t\t%f\t\t%d\n"%(cen1,cen_hist[i][0],cen2,cen_hist[i][1]))
#
#closing files
file1.close()
file2.close()
file3.close()
#file4.close()
file5.close()
file6.close()
file7.close()
print("--- %s seconds ---" % (time.time() - start_time))                              
