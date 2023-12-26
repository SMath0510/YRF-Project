import os
import sys
import string
import math
import time
import numpy as np
import extract_pdb
import matplotlib.pyplot as plt
import networkx as nx
from mdtraj.formats import XTCTrajectoryFile
import copy

print("importing done")

start_time = time.time()

#input filenames
xtcname='bulk_tip3p_all_part1.xtc' 
groname='bulk_tip3p_all.gro' #centered gro file given as input

print("filenames read")

#output filenames
extnname='.dat'
file1name='shaun_outputs/Network_Stats_'+extnname
file2name='shaun_outputs/molhb_histogram_'+extnname
file3name='shaun_outputs/hbweight_histogram_'+extnname
file5name='shaun_outputs/Pathlength_histogram_'+extnname
file6name='shaun_outputs/PerFrame_HistCount_'+extnname
file7name='shaun_outputs/Centrality_histogram_'+extnname

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

deltat = 1 #changed from 1.0
nframes = 5 #changed from 2
frame_start=2000.0
coord=[]
atomname=[]
resid=[]
R_a=[]

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

## reading the traj file
traj_file = XTCTrajectoryFile(xtcname)
xyz,tot_time,step,box = traj_file.read()

box_len =[]
for i in range(0,nframes):
	box_len.append(box[i][0][0])
#########################

'''
Time starts at 2000 and ends at 3000
just 1000 frames. Thats all.
'''

for nf in range(0, nframes):
	ftime = frame_start + nf * deltat
	frame_file_name = f"frame_{nf}.pdb"
	cmd = f'echo 1 | gmx -quiet trjconv -f {xtcname} -s {groname} -o {frame_file_name} -b {ftime} -e {ftime + deltat}'
	# Add print statement to check the command
	print(cmd)
	os.system(cmd)

	ffile=open(frame_file_name,'r')
	lcount = 0
	natoms = 0
	coord=[]
    
	x_main=[]
	y_main=[]
	z_main=[]   
	
	C=extract_pdb.extractor(frame_file_name)

	index=C[1,:]
	chain=C[2,:]
	res=C[3,:]
	res_num=C[4,:]
	x_coord=C[5,:]
	y_coord=C[6,:]
	z_coord=C[7,:]
	
	x_coord=extract_pdb.str_to_float(x_coord)
	y_coord=extract_pdb.str_to_float(y_coord)
	z_coord=extract_pdb.str_to_float(z_coord)
	res_num=extract_pdb.str_to_int(res_num)

	coord = []

	n_atom = len(x_coord)
	nres = extract_pdb.count_unique(res_num)
	for i in range(n_atom):
		coord.append([x_coord[i], y_coord[i], z_coord[i]])

	donors = []
	acceptors = []
	for index in range(n_atom):
		if chain[index] == 'OW':
			donors.append(index)
			acceptors.append(index)
	ndon = len(donors)
	nacc = len(acceptors)
	

#
	# print(f'Frame: {nf} \nX: {x_coord[:10]}\nY: {y_coord[:10]}\nZ: {z_coord[:10]}\nRes: {res[:10]}\nResNum: {res_num[:10]}\nChain: {chain[:10]}\nnres = {nres}\nnatoms = {n_atom}')
	#Out=main_calculation.inner_surface(x_coord,y_coord,z_coord,res,res_num,chain,ftime,5)
	cmd2 = f"rm {frame_file_name}"
	os.system(cmd2)
	
    	### The analysis part starts now

	print("Working on frame " + str(nf))
	count = 0
	for i in range (0,nres):
		rname = res[i]
		hbnet.add_node(i,name=rname)
	#
	for index1 in donors:
		rid1 = res_num[index1]
		for index2 in acceptors:
			rid2 = res_num[index2]
			if rid1 != rid2:
				DAdist_sq = extract_pdb.distance_sq(coord,box_len,nf,index1,index2)
				if DAdist_sq <= 0.1225:
					angle_h1 = extract_pdb.angle(coord,box_len,nf,index1+1,index1,index2)
					if angle_h1 <= 30:
						count = count + 1
						wt = 1.0
						hbnet.add_edge(rid1,rid2,weight=wt)
					angle_h2 = extract_pdb.angle(coord,box_len,nf,index1+2,index1,index2)
					if angle_h2 <= 30:
						count = count + 1
						wt = 1.0 
						hbnet.add_edge(rid1,rid2,weight=wt)
					#Angle loops
				#heavy atom distance loop
		#end of acceptor loop
	#end of donor loop
	#in_degree = number of h-bonds accepted by the molecule
	#out_degree = number of h-bonds donated by the molecule
	print("Network created for frame %d"%(nf))
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
	file1.write("%d\t\t%d\t\t%d\t\t%d\n"%(nf,count,hbnet.number_of_nodes(), hbnet.number_of_edges()))
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
	file6.write("%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\n"%(nf,waccsum,wdonsum,whbsum,pcount,mecount,mecount2,nrcount1))
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

