#USE PYTHON 3.* 
import os,sys,string
import math,time
import numpy as np
import extract_pdb
import matplotlib.pyplot as plt
def min_rad(x,y,z,ftime):
	'''
	min_z = []
	min_z = min(z)
	z_min = float(min_z[0])
	print ("Minimum Z value:", min_z)

	max_z = []
	max_z = max(z)
	z_max = float(max_z[0])
	'''
	z_min = np.min(z)
	z_max = np.max(z)
	print(z_min, z_max)
	z_value = [2]
	z_Max2=[]
	all_t_r=[]
	for j in z_value:
		b = np.arange(np.floor(z_min),np.ceil( z_max), j)
		rad_append = []
		z_rad = []
		z_max2=[]
		r=[]

		for dz in range(len(b)-1):
			x_st = []
			y_st = []
			z_st = []
			Z_mid=b[dz]
			for k in range(len(z)) :
				if(z[k]>=b[dz] and z[k]<=b[dz+1]):
					x_st.append(x[k])
					y_st.append(y[k])
					z_st.append(z[k])
			r=[0.0]*len(x_st)
			for i in range(len(x_st)):
				r[i]=((x_st[i]**2)+(y_st[i]**2))**0.5
			ind_min=np.argmin(r)
		
			z_max2.append([min(r),z[ind_min],ftime])
		z_Max2.append(z_max2)
	return z_Max2

'''
#
def distance_sq(n1,n2):
    coor1 = coord[n1]
    coor2 = coord[n2]
    dx = coor2[0] - coor1[0]
    dy = coor2[1] - coor1[1]
    dz = coor2[2] - coor1[2]
    dx = dx - (boxl * round(dx/boxl))
    dy = dy - (boxl * round(dy/boxl))
    dz = dz - (boxl * round(dz/boxl))
    new_dist_sq = (dx)**2 + (dy)**2 + (dz)**2
    return new_dist_sq
#
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
#
def dot(A,B):
    dot = A[0]*B[0] + A[1]*B[1] + A[2]*B[2]
    return dot

def angle(n1,n2,n3):
    coor1 = coord[n1]
    coor2 = coord[n2]
    coor3 = get_pbcacc(n2,n3)
    vec1 = vecsub(coor1,coor2)
    vec2 = vecsub(coor3,coor2)
    x = (dot(vec1,vec2))/(mag(vec1)*mag(vec2))
    angle_rad = (np.arccos(x))  
    angle = np.rad2deg(angle_rad)      # angle in degrees
    return angle
'''

start_time = time.time()
#
#input filenames
#
xtcname='md_noPBC.xtc'
groname='md_0_1.gro'
#indfile='system.ndx'
deltat = 100 #changed from 1.0
nframes=11 #changed from 2
frame_start=0.0
coord=[]
atomname=[]
resid=[]
R_a=[]
#########################

for nf in range(0,nframes):
	ftime=str(frame_start+nf*deltat)
	cmd="echo 1 | gmx trjconv -f "+xtcname+" -s "+groname+" -o frame.pdb -b "+ftime+" -e "+ftime
	os.system(cmd)
	ffile=open("frame.pdb",'r')
	lcount = 0
	natoms = 0
	coord=[]
    
	x_main=[]
	y_main=[]
	z_main=[]      
########################################################################################################################################################################            
	C=extract_pdb.coord_extractor("frame.pdb")
	xx1=C[0,:]
	yy1=C[1,:]
	zz1=C[2,:]
	rr= min_rad(xx1,yy1,zz1,float(ftime))
	ffile.close()
	R_a.append(rr)
	cmd2="rm frame.pdb"
	os.system(cmd2)
    #do the analysis you want to do with the coordinates 
    #
    #
    #
for i in range(11):
	R_a1=R_a[:][:][i][:]
	R_a2=R_a1[0][:]
	R_a3=np.array(R_a2)
	R_a4=np.transpose(R_a3)
	r_val=R_a4[0,:]
	z_val=R_a4[1,:]
	plt.plot(z_val,r_val,'o')
	plt.show()
#frame lofro

#for i in range(10):
	
		
