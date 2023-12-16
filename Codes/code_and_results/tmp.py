import time
start_time = time.time()
#USE PYTHON 3.* and Networkx 2.3
import math
import numpy as np
import mdtraj as md
from mdtraj.formats import XTCTrajectoryFile
import networkx as nx
import matplotlib.pyplot as plt
import concurrent.futures ## For implementing parallelization in centrality
from scipy.spatial import cKDTree ## Using kdtree data structure for optimisation
from concurrent.futures import ProcessPoolExecutor
import copy

print("Version NetworkX = " + nx.__version__)

#input filenames
xtcname = 'water_traj/bulk_tip3p_all_part1.xtc'
groname = 'water_traj/bulk_tip3p_all.gro'

#mode of execution
approximate = False
threading = False
multiprocessing = True
community_analysis = True
optimize_loops = True
test = False
stop = False

shortest_path_mode = "unweighted-all-pair-shortest-paths"
# shortest_path_mode = "all-pair-dijkstras"

print("filenames read")
#output filenames
head = "shaun_outputs/"
if threading:
    head = head + "with_threading_"
elif multiprocessing:
    head = head + "with_multiprocessing_"
else:
    head = head + "without_threading_"
    
if approximate:
    head = head + "with_approximation/"
else:
    head = head + "without_approximation/"
    
extnname='.dat'
file1name=head + 'Network_Stats_'+extnname
file2name=head + 'molhb_histogram_'+extnname
file3name=head + 'hbweight_histogram_'+extnname
file5name=head + 'Pathlength_histogram_'+extnname
file6name=head + 'PerFrame_HistCount_'+extnname
file7name=head + 'Centrality_histogram_'+extnname
file9name = head + "Average Bond Strength.csv"
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
if test: 
    file8 = open("test_shaun.txt", "w")
file9 = open(file9name, "w")
file9.write(f'Key, Avg. Strength\n')
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
kd_tree_threshold = 0.4

old_hbond_list = []
latest_hbond_list = []
average_hbond_life_mapping = {}


def update_mapping(frame_no, num_nodes, max_frames = 2):
    print(f'Updating Mapping: {frame_no+1}/{max_frames}')
    for key in latest_hbond_list:
        
        if key not in average_hbond_life_mapping.keys():
            average_hbond_life_mapping[key] = {'Start': None, 'End': None, 'Time_List': []}
        
        if (key not in old_hbond_list) or (average_hbond_life_mapping[key]['Start'] is None):
            average_hbond_life_mapping[key]['Start'] = frame_no
        
        average_hbond_life_mapping[key]['End'] = frame_no
        
            
    for key in old_hbond_list:
        
        if key not in average_hbond_life_mapping.keys():
            average_hbond_life_mapping[key] = {'Start': None, 'End': None, 'Time_List': []}
        
        if (key not in latest_hbond_list):
            average_hbond_life_mapping[key]['End'] = frame_no
            average_hbond_life_mapping[key]['Time_List'].append(
                average_hbond_life_mapping[key]['End'] - average_hbond_life_mapping[key]['Start']
            )

    if ((frame_no + 1) == max_frames):
        for key in latest_hbond_list:
            average_hbond_life_mapping[key]['End'] = frame_no + 1
            average_hbond_life_mapping[key]['Time_List'].append(
                average_hbond_life_mapping[key]['End'] - average_hbond_life_mapping[key]['Start']
            )

                
            
#

# Defining some centrality functions (with optimisation)
def calculate_closeness_centrality(graph):
    return nx.closeness_centrality(graph, distance="weight")

def calculate_betweenness_centrality(graph):
    return nx.betweenness_centrality(graph, weight="weight")

def calculate_degree_centrality(graph):
    return nx.degree_centrality(graph)

def calculate_eigenvector_centrality(graph):
    return nx.eigenvector_centrality(graph, weight="weight", max_iter=1000)

def calculate_katz_centrality(graph):
    return nx.katz_centrality(graph, weight="weight")

def calculate_pagerank_centrality(graph):
    return nx.pagerank(graph, weight="weight")
        
#
print("No. of frames: %d"%(n_frame))

num_iter = min(10, n_frame)
for frame in range(0,num_iter):
    if test: 
        print(coord[frame].shape)
        print(acceptors.shape)
        print(acceptors[0])
    # if stop:
    #     break
    if optimize_loops: 
        kdtree = cKDTree(coord[frame]) 
    count =0
    for i in range (0,nres) :
        rname = topol.residue(i).name
        hbnet.add_node(i,name=rname)
        if(rname == "HOH"):
            watnet.add_node(i)
	#
    ## Applying an optimisation that converts the code from O(n2) to O(nlogn)
    for donor_index in donors:
        index1 = donor_index
        rid1 = topol.atom(index1).residue.index
        if optimize_loops:
           potential_acceptors_indices = [acc_index for acc_index in acceptors if distance_sq(frame, index1, acc_index) <= kd_tree_threshold]
        else:
            potential_acceptors_indices = acceptors
        for acceptor_index in potential_acceptors_indices:
            index2 = acceptor_index
            rid2 = topol.atom(index2).residue.index
			#if index1 != index2:
            if rid1 != rid2:
                DAdist_sq = distance_sq(frame,index1,index2)
                if DAdist_sq <= 0.1225:
                    if test:
                        file8.write("Donor Index: %d, Acceptor Index: %d, DistSquared: %d\n" %( donor_index, acceptor_index, DAdist_sq))
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
       
    # Store copies as deep copies
    old_hbond_list = copy.deepcopy(latest_hbond_list)

    # Update latest_hbond_list (assuming hbnet.edges is a list)
    latest_hbond_list = copy.deepcopy(hbnet.edges)
    print(f'After first update, OldSize: {len(old_hbond_list)}, NewSize: {len(latest_hbond_list)}')
    # print(topol.n_atoms)
    update_mapping(frame_no=frame, num_nodes=topol.n_atoms, max_frames=num_iter)
    
    print("after path updations")
    print("--- %s seconds ---" % (time.time() - start_time))   
    
    if stop:
        continue
                            
    wwcount = 0
    # weight distribution
    for n, nbrsdict in hbnet.adj.items():
        for nbr, keydict in nbrsdict.items():
            for key, eattr in keydict.items():
                wt = eattr['weight']
                bin = int(wt/wbin_width)
                weight_hist[bin] = weight_hist[bin] + 1
                wwcount = wwcount + 1
                
    # path length properties using Floyd-Warshall
    if shortest_path_mode == "floyd-warshall":
        path = dict(nx.floyd_warshall(hbnet))   ## Applying Floyd Warshall on the graph
    elif shortest_path_mode == "all-pair-dijkstras":
        path = dict(nx.all_pairs_dijkstra_path_length(hbnet))
    elif shortest_path_mode == "unweighted-all-pair-shortest-paths":
        path = dict(nx.all_pairs_shortest_path_length(hbnet))

    print("after Shortest path %s path for 1st network", shortest_path_mode)
    print("--- %s seconds ---" % (time.time() - start_time))
                        
    #
    waccsum = 0
    wdonsum = 0
    whbsum = 0
    pcount = 0
    mecount=0
    mecount2=0
    file1.write("%d\t\t%d\t\t%d\t\t%d\n"%(frame,count,hbnet.number_of_nodes(), hbnet.number_of_edges()))
	
    
    waccsum = sum(hbnet.in_degree(nd) for nd in hbnet.nodes)
    wdonsum = sum(hbnet.out_degree(nd) for nd in hbnet.nodes)
    whbsum = sum(hbnet.degree(nd) for nd in hbnet.nodes)

    mecount = sum(1 for nd in hbnet.nodes for nd2 in hbnet.nodes if nd != nd2 and hbnet.number_of_edges(nd, nd2) > 1)
    mecount2 = sum(1 for nd in hbnet.nodes for nd2 in hbnet.nodes if nd != nd2 and hbnet.number_of_edges(nd, nd2) > 2)

    for value, hist_index in zip([hbnet.in_degree(nd) for nd in hbnet.nodes], [0, 1, 2]):
        if value < nhb:
            nhb_hist[value][hist_index] += 1
        else:
            print("Number of Hbonds exceeds %d" % nhb)

		#
    #hbnet node loop
    print("after loop network for 1st network")
    print("--- %s seconds ---" % (time.time() - start_time))   
    
    if shortest_path_mode == "floyd-warshall":
        for ss in path:
            for tt, nconn in path[ss].items():
                if math.isfinite(nconn) and nconn > 0:
                    pcount += 1
                    if nconn < mpl:
                        path_hist[int(nconn)] += 1  # Convert to int to get the number of edges
                    else:
                        print("Total path length %d exceeds maximum path length of %d" % (nconn, mpl))
                        # print frame, ss, tt, nconn, path[ss][tt]
    elif shortest_path_mode == "all-pair-dijkstras":
        for ss in path.keys():
            for tt, path_length in path[ss].items():
                if path_length > 0:  # Exclude self-loops
                    path_length = int(path_length)
                    pcount += 1
                    if path_length < mpl:
                        path_hist[path_length] += 1
                    else:
                        print("Total path length %d exceeds maximum path length of %d" % (path_length, mpl))
                        # print frame, ss, tt, path_length
    elif shortest_path_mode == "unweighted-all-pair-shortest-paths":
        for ss in path:
            for tt, path_length in path[ss].items():
                if path_length > 0:  # Exclude self-loops
                    pcount += 1
                    if path_length < mpl:
                        path_hist[int(path_length)] += 1  # Convert to int to use as an index
                    else:
                        print("Total path length %d exceeds maximum path length of %d" % (path_length, mpl))
                        # print frame, ss, tt, path_length
                
    print("after path histogram for 1st network")
    print("--- %s seconds ---" % (time.time() - start_time))     
                            
    #Centrality Analysis
    #clcen - closeness centrality gives the product of ratio of nodes that can reach a given node (inward paths) and the inverse of the average distance to the given node. closeness distance is based on the incoming distance 
    #btcen - between centrality

    # Parallelize centrality calculations
    if threading:
        with concurrent.futures.ThreadPoolExecutor() as executor:
            clcen_future = executor.submit(calculate_closeness_centrality, hbnet)
            btcen_future = executor.submit(calculate_betweenness_centrality, hbnet)

            clcen1 = clcen_future.result()
            btcen1 = btcen_future.result()
        ## Has it been parallelized (NO - Doesnt work well on CPU bound tasks)
    elif multiprocessing:
        with ProcessPoolExecutor() as executor:
            clcen_future = executor.submit(calculate_closeness_centrality, hbnet)
            btcen_future = executor.submit(calculate_betweenness_centrality, hbnet)

            clcen1 = clcen_future.result()
            btcen1 = btcen_future.result()
    else:
        clcen1 = calculate_closeness_centrality(hbnet)
        btcen1 = calculate_betweenness_centrality(hbnet)
    # Avoiding redundant storage and repeated calculations
    hnodes = hbnet.nodes()
    nrcount1 = 0

    for i, nd in enumerate(hnodes):
        in_degree_nd = hbnet.in_degree(nd)  # Store in_degree result
        if in_degree_nd == 0:
            nrcount1 += 1

        clcen_nd = clcen1[nd]
        btbin_nd = btcen1[nd]

        clbin = int(clcen_nd / clbin_width)
        btbin = int(btbin_nd / btbin_width)

        if clbin < mcen:
            cen_hist[clbin][0] += 1
        else:
            print(f"Water {nd} Closeness {clcen_nd} binid {clbin} is greater than {mcen}")

        if btbin < mcen:
            cen_hist[btbin][1] += 1
        else:
            print(f"Water {nd} Betweenness {btbin_nd} binid {btbin} is greater than {mcen}")


	#plotting
	#if(frame==0):
	#	nx.draw(hbnet,nx.circular_layout(hbnet),nodecolor='r', edge_color='b')
	#	plt.show()
	#
    print("after centrality histogram for networks")
    print("--- %s seconds ---" % (time.time() - start_time))                              
    #print count,len(nodes),minwt,maxwt,x

    if community_analysis:

        # Community Detection
        communities = nx.algorithms.community.greedy_modularity_communities(hbnet, cutoff=1)

        # Create a mapping of node to community index
        community_mapping = {node: i for i, community in enumerate(communities) for node in community}

        # Add community information to nodes
        nx.set_node_attributes(hbnet, community_mapping, "community")

        # Filter nodes based on degree centrality (optional)
        # degree_threshold = 10
        # filtered_nodes = [node for node, degree in hbnet.degree() if degree >= degree_threshold]
        # filtered_hbnet = hbnet.subgraph(filtered_nodes)

        filtered_hbnet = hbnet
        # Draw the graph with nodes colored by community
        plt.figure(figsize=(10, 8))
        pos = nx.spring_layout(filtered_hbnet)  # You can use a different layout algorithm if needed
        node_colors = [filtered_hbnet.nodes[node]["community"] for node in filtered_hbnet]
        nx.draw(filtered_hbnet, pos, node_color=node_colors, cmap=plt.cm.get_cmap("viridis"), with_labels=True, font_size=8)

        # Save the plot as an image file (adjust the filename and format as needed)
        
        plt.savefig(f"{head}_network_plot_frame{frame}.png")
        # plt.show()
        # Clearing networks
        print("after motif and community analysis for networks")
        print("--- %s seconds ---" % (time.time() - start_time))

    hbnet.clear()

#frame loop
#Normalization could use 0,1,2,6,7,8 - #water nodes * nframes; 3,4,5 - #sugar nodes*nframes

# file9.write(f'{average_hbond_life_mapping.keys()}')
    
for key in average_hbond_life_mapping.keys():
    len_ = len(average_hbond_life_mapping[key]['Time_List'])
    if(len_ == 0):
        continue
    # file9.write(f'Nonzero Key: {key}')
    sum_ = sum(average_hbond_life_mapping[key]['Time_List'])
    average = sum_ / len_
    
    file9.write(f'{key}, {average}\n')
    
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
if test: 
    file8.close()
file9.close()