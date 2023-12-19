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
from concurrent.futures import ThreadPoolExecutor
import copy
from matplotlib import colormaps

print("Version NetworkX = " + nx.__version__)

#input filenames
xtcname = 'water_traj/bulk_tip3p_all_part1.xtc'
groname = 'water_traj/bulk_tip3p_all.gro'

# Execution Parameters
cmd_print = False
optimize_loops = False
subgraph_threading = False

approximate = False
approximation_factor = 0.7
'''
Approximation factor <= 1
1 -> perfect accuracy
'''

threading = False
'''
Useful for File Operations and I/O based operations
'''

multiprocessing = False
'''
Useful for CPU intensive operations
'''


weighted = None
'''
Choices: None, 
         "weight"
'''

default_shortest_path_mode = "unweighted-all-pair-shortest-path"
'''
Choices: "unweighted-all-pair-shortest-path", [preferred if unweighted graph]
         "floyd-warshall", [preferred if weighted dense graph]
         "all-pair-dijkstras" [preferred if weighted sparse graph]
'''

print("filenames read")
#output filenames
head = "shaun_outputs/"

status = ""
if subgraph_threading:
    status += "subgraphtheading_"
if threading:
    status += "threading_"
if multiprocessing:
    status += "multiprocessing_"
if optimize_loops:
    status += "kd_optim_"
status += f"{default_shortest_path_mode}"
    
extnname='.dat'
file1name=head + 'Network_Stats_'+extnname
file2name=head + 'molhb_histogram_'+extnname
file3name=head + 'hbweight_histogram_'+extnname
file5name=head + 'Pathlength_histogram_'+extnname
file6name=head + 'PerFrame_HistCount_'+extnname
file7name=head + 'Centrality_histogram_'+extnname
file8name = head + 'Average_Bond_Strength.csv'
file9name = head + 'Community_Analysis.txt'
file10name = head + f'Time_Analysis_{status}.csv'
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

## Mathematical Operations

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
##

## File Handling 

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
file8 = open(file8name, "w")
file8.write(f'Key, Avg. Strength\n')
file9 = open(file9name, "w")
file10 = open(file10name, "w")
file10.write(f"Graph Construction, Shortest Path, Centrality, Community\n")

##

## Graph Creation
hbnet = nx.MultiDiGraph()
# watnet = nx.MultiDiGraph() 

## Defining the graph based constraints
wmin = 0.4 # 0.14/0.35 as 3.5 A is the h-bond cut-off distance
wmax = 0.6 # 0.14/0.23 as 2.3 A is the absolute hard sphere radius because gr of wat-wat oxygens is absolutely 0.0 upto 2.3 A
cosmin = np.cos(30*np.pi/180)
cosnorm = 1.0 - cosmin

## Defining histogram specific constraints
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

## KD - Tree Loop Optimisation distance threshold
kd_tree_threshold = 0.4

## H-bond lifetime analysis global params
old_hbond_list = []
latest_hbond_list = []
average_hbond_life_mapping = {}

## Useful Functions 

def matrix_to_dict(mat):
    dict_ = {}
    rows_ = mat.shape[0]  
    cols_ = mat.shape[1] 
    for i in range(rows_):
        entry_ = {}
        for j in range(cols_):
            entry_.update({j: mat[i, j]})  
        dict_.update({i: entry_})

    return dict_  

def calculate_shortest_path(G, mode = default_shortest_path_mode):
    if mode == "floyd-warshall":
        path = matrix_to_dict(nx.floyd_warshall_numpy(G))  
    elif mode == "all-pair-dijkstras":
        path = dict(nx.all_pairs_dijkstra_path_length(G))
    elif mode == "unweighted-all-pair-shortest-paths":
        path = dict(nx.all_pairs_shortest_path_length(G))
    
    return path

def calculate_centrality(G, mode = "closeness", weight_mode = weighted, approximate_status = approximate, factor = approximation_factor):
    num_nodes = G.number_of_nodes()
    if not approximate_status:
        factor = 1
    k_factor = int(factor * num_nodes)
    if mode == "closeness":
        centrality_value = nx.closeness_centrality(G, distance=weight_mode)
    elif mode == "betweeness":
        centrality_value = nx.betweenness_centrality(G, weight=weight_mode, k = k_factor)
    elif mode == "degree":
        centrality_value = nx.degree_centrality(G)
    elif mode == "eigen":
        centrality_value = nx.eigenvector_centrality_numpy(G, weight=weight_mode, max_iter=1000)
    elif mode == "katz":
        centrality_value = nx.katz_centrality(G, weight=weight_mode)
    elif mode == "pagerank":
        centrality_value = nx.pagerank(G, weight=weight_mode)
    return centrality_value
    
def calculate_centrality_for_subgraphs(hbnet, centrality_type):
    connected_components = list(nx.weakly_connected_components(hbnet))
    subgraphs = [hbnet.subgraph(component).copy() for component in connected_components]

    # Use parallel processing to calculate centrality for each connected component
    with ProcessPoolExecutor() as executor:
        try: 
            centrality_list = list(executor.map(calculate_centrality, subgraphs, [centrality_type]*len(subgraphs)))
        except:
            centrality_list = [calculate_centrality(hbnet, centrality_type)]

    # Combine results into a single dictionary
    centrality_combined = {}
    for node_centrality_detail in centrality_list:
        centrality_combined.update(node_centrality_detail)

    centrality_combined = dict(sorted(centrality_combined.items()))
    return centrality_combined

def update_mapping(frame_no, max_frames=2):
    print(f'Updating Mapping: {frame_no+1}/{max_frames}')
    
    for key in set(latest_hbond_list) - set(old_hbond_list):
        average_hbond_life_mapping.setdefault(key, {'Start': frame_no, 'End': frame_no, 'Time_List': []})
        
    for key in set(old_hbond_list) - set(latest_hbond_list):
        average_hbond_life_mapping[key]['End'] = frame_no
        average_hbond_life_mapping[key]['Time_List'].append(
            average_hbond_life_mapping[key]['End'] - average_hbond_life_mapping[key]['Start']
        )

    if frame_no + 1 == max_frames:
        for key in set(latest_hbond_list):
            average_hbond_life_mapping[key]['End'] = frame_no + 1
            average_hbond_life_mapping[key]['Time_List'].append(
                average_hbond_life_mapping[key]['End'] - average_hbond_life_mapping[key]['Start']
            )


        
#
print("No. of frames: %d"%(n_frame))

for i in range (0,nres) :
    rname = topol.residue(i).name
    hbnet.add_node(i,name=rname)
    if(rname == "HOH"):
        ## Can add it as a node in water network
        pass

## For testing purposes    
num_iter = min(5, n_frame)

for frame in range(0,num_iter):
    local_start_time = time.time()
    if optimize_loops: 
        kdtree = cKDTree(coord[frame]) 
    count =0

    ## Applying an optimisation that converts the code from O(n2) to O(nlogn)
    for donor_index in donors:
        index1 = donor_index
        rid1 = topol.atom(index1).residue.index
        
        ## Optimization
        if optimize_loops:
           potential_acceptors_indices = [acc_index for acc_index in acceptors if distance_sq(frame, index1, acc_index) <= kd_tree_threshold]
        else:
            potential_acceptors_indices = acceptors
            
        for acceptor_index in potential_acceptors_indices:
            index2 = acceptor_index
            rid2 = topol.atom(index2).residue.index
            
            if rid1 != rid2:
                DAdist_sq = distance_sq(frame,index1,index2)
                if DAdist_sq <= 0.1225:
                    angle_h1 = angle(frame,index1+1,index1,index2)
                    if angle_h1 <= 30:
                        count = count + 1
                        wt = 1.0
                        hbnet.add_edge(rid1,rid2,weight=wt)
                    angle_h2 = angle(frame,index1+2,index1,index2)
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
    if cmd_print:
        print("Network created for frame %d"%(frame))
        print("--- %s seconds ---" % (time.time() - start_time))
    
    graph_creation_time = time.time() - local_start_time   
                            
    wwcount = 0
    # weight distribution
    for n, nbrsdict in hbnet.adj.items():
        for nbr, keydict in nbrsdict.items():
            for key, eattr in keydict.items():
                wt = eattr['weight']
                bin = int(wt/wbin_width)
                weight_hist[bin] = weight_hist[bin] + 1
                wwcount = wwcount + 1
                
        ## H-bond Time Analysis
    local_start_time = time.time()
    old_hbond_list = copy.deepcopy(latest_hbond_list)
    latest_hbond_list = copy.deepcopy(hbnet.edges)
    
    update_mapping(frame_no=frame, max_frames=num_iter)
    
    if cmd_print:
        print("After the Hydrogen Bond Time Analysis")
        print("--- %s seconds ---" % (time.time() - start_time))   
    
    h_bond_time_analysis = time.time() - local_start_time
    
    ## Shortest Path Calculations
    local_start_time = time.time()
    
    shortest_path_mode = "unweighted-all-pair-shortest-paths"
    path = calculate_shortest_path(hbnet, mode = shortest_path_mode)

    if cmd_print:
        print("after Shortest path %s method for 1st network" % shortest_path_mode)
        print("--- %s seconds ---" % (time.time() - start_time))

    shortest_path_time = time.time() - local_start_time             
    
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
    if cmd_print:
        print("after loop network for 1st network")
        print("--- %s seconds ---" % (time.time() - start_time))   
    
    for ss in path.keys():
            for tt, path_length in path[ss].items():
                if math.isfinite(path_length) and path_length > 0:
                    path_length = int(path_length)
                    pcount += 1
                    if path_length < mpl:
                        path_hist[int(path_length)] += 1  # Convert to int to get the number of edges
                    else:
                        print("Total path length %d exceeds maximum path length of %d" % (path_length, mpl))
        
    
    if cmd_print:
        print("after path histogram for 1st network")
        print("--- %s seconds ---" % (time.time() - start_time))     
                      
    ## Centrality Analysis
    #clcen - closeness centrality gives the product of ratio of nodes that can reach a given node (inward paths) and the inverse of the average distance to the given node. closeness distance is based on the incoming distance 
    #btcen - between centrality

    # Parallelize centrality calculations
    
    local_start_time = time.time()
    centrality_types = ["closeness", "betweeness", "eigen"]
    
    if subgraph_threading:
        clcen1 = calculate_centrality_for_subgraphs(hbnet, "closeness")
        btcen1 = calculate_centrality_for_subgraphs(hbnet, "betweeness")
        egcen1 = calculate_centrality_for_subgraphs(hbnet, "eigen")
        
    elif multiprocessing:
        with ProcessPoolExecutor() as executor:
            results = list(executor.map(calculate_centrality, [hbnet]*3, centrality_types))
        clcen1, btcen1, egcen1 = results
    
    elif threading:
        with ThreadPoolExecutor() as executor:
            results = list(executor.map(calculate_centrality, [hbnet]*3, centrality_types))
        clcen1, btcen1, egcen1 = results
        
    else:
        clcen1 = calculate_centrality(hbnet, "closeness")
        btcen1 = calculate_centrality(hbnet, "betweeness")
        egcen1 = calculate_centrality(hbnet, "eigen")
        
    centrality_analysis_time = time.time() - local_start_time
    
    #
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
    if cmd_print:
        print("after centrality histogram for networks")
        print("--- %s seconds ---" % (time.time() - start_time))                              
    #print count,len(nodes),minwt,maxwt,x

    
    ## Eigen Value Analysis
    # Add eigenvector centrality as node attribute
    nx.set_node_attributes(hbnet, egcen1, 'eigenvector_centrality')

    # Visualize the graph with node colors based on eigenvector centrality
    node_colors = [egcen1[node] for node in hbnet.nodes()]
    pos = nx.spring_layout(hbnet)
    
    nx.draw(hbnet, pos, node_color=node_colors, cmap=colormaps["plasma"], with_labels=True, font_size=4, node_size = 200)
    plt.title("Graph with Eigenvector Centrality")
    plt.savefig(f"{head}_network_plot_frame_eigen_{frame}.png")


    ## Community Analysis
    local_start_time = time.time()
    communities = nx.algorithms.community.greedy_modularity_communities(hbnet, cutoff=1)
    community_analysis_time = time.time() - local_start_time
    # Create a mapping of node to community index
    community_mapping = {node: i for i, community in enumerate(communities) for node in community}

    ## Plotting based details
    
    nx.set_node_attributes(hbnet, community_mapping, "community")

    # Filter nodes based on degree centrality (optional)
    average_degree = 0
    max_degree = 0
    for node, degree in hbnet.degree():
        average_degree += degree
        max_degree = max(degree, max_degree)
    average_degree /= hbnet.number_of_nodes()
    degree_threshold = min(average_degree * 1.1, max_degree*0.8)
    print(f'Average Degree = {average_degree}, Max Degree = {max_degree}')
    filtered_nodes = [node for node, degree in hbnet.degree() if degree >= degree_threshold]
    filtered_hbnet = hbnet.subgraph(filtered_nodes)
    # filtered_hbnet = hbnet
    
    plt.figure(figsize=(10, 8))
    pos = nx.spring_layout(filtered_hbnet)  # You can use a different layout algorithm if needed
    node_colors = [filtered_hbnet.nodes[node]["community"] for node in filtered_hbnet]
    nx.draw(filtered_hbnet, pos, node_color=node_colors, cmap=colormaps["viridis"], with_labels=True, font_size=4, node_size = 200)

    plt.savefig(f"{head}_network_plot_frame_comm1_{frame}.png")
        
    # Print the detected communities
    for i, community in enumerate(communities):
        file9.write(f"Community {i + 1}: {list(community)}")
    
    if cmd_print:
        print("after community analysis for networks")
        print("--- %s seconds ---" % (time.time() - start_time))
    file10.write(f'{graph_creation_time}, {shortest_path_time}, {centrality_analysis_time}, {community_analysis_time}\n')
    hbnet.clear_edges()

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
    
    file8.write(f'{key}, {average}\n')
    
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
file5.close()
file6.close()
file7.close()
file8.close()
file9.close()
file10.close()