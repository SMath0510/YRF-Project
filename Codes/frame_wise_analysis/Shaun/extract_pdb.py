import math
import numpy as np
import re


def extractor(subunit):
    	# Initialize an empty list to store coordinates
	features = []

	# Open the file and read line by line
	with open(subunit, 'r') as f:
		for line in f:
			# Check if the line starts with "ATOM"
			if line.startswith('ATOM'):
				# Split the line using both tabs and spaces as delimiters
				char_list = re.split(r'\t| ', line.strip())

				# Remove empty strings and spaces from the list
				char_list = [item for item in char_list if item]

				# Append it to the list of coordinates
				features.append(char_list)

	# Convert the list of coordinates to a numpy array
	feature_array = np.array(features).T

	return feature_array
	

def str_to_float(X):
    	return [float(x) for x in X]

def str_to_int(X):
    	return [int(x) for x in X]

def count_unique(input_list):
	empty_list = []
	for ele in input_list:
		if(ele not in empty_list):
			empty_list.append(ele) 
	return len(empty_list)

def anint(value):
	return round(value)

def distance_sq(coord, box_len,frame,n1,n2):
	coor1 = coord[n1]
	coor2 = coord[n2]
	boxl = box_len[frame]
	dx = coor2[0] - coor1[0]
	dy = coor2[1] - coor1[1]
	dz = coor2[2] - coor1[2]
	dx = dx - (boxl * anint(dx/(boxl)))
	dy = dy - (boxl * anint(dy/(boxl)))
	dz = dz - (boxl * anint(dz/(boxl)))
	new_dist_sq = (dx)**2 + (dy)**2 + (dz)**2
	return new_dist_sq

def get_pbcacc(coord,box_len,frame,n1,n2):
	coor1 = coord[n1]
	coor2 = coord[n2]
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

def angle(coord,box_len,frame,n1,n2,n3):
	coor1 = coord[n1]
	coor2 = coord[n2]
	coor3 = get_pbcacc(coord,box_len,frame,n2,n3)
	vec1 = vecsub(coor1,coor2)
	vec2 = vecsub(coor3,coor2)
	x = (dot(vec1,vec2))/(mag(vec1)*mag(vec2))
	angle_rad = (np.arccos(x))  
	angle = np.rad2deg(angle_rad)      # angle in degrees
	return angle
