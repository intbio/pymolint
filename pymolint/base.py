"""
This is a module which provides basic analysis routines
for interactions in molecular systems, it is independent from VMD or any software.
exlucde_bonded might not work as expected
"""

from scipy.spatial import KDTree
from scipy.spatial import cKDTree
import numpy as np
from collections import OrderedDict

__author__="Alexey Shaytan"


def find_contacts(A_coordinates,A_indices,B_coordinates,B_indices,d_threshold=3.9,exclude_bonded=False,half_matrix=True):
	""" Function finds contacts between atoms in groups A and B.
	If exclude_bonded is set to True, we assume that set A = set B,
	and we need to exclude bonded neighbors up to three bonds away (inclusive: 1-4 interactions are also excluded).
	Self-interaction has also to be excluded.
    and only half-matrix is output by default (change parameter to override)
	the bonded threshold is set to colser than 1.7 A 
	Atoms must be heavy!
	Parameters
	----------
	A_coordinates - list of tuples with (x,y,z) of group A.
	A_indices - list of indices of atoms corresponding to the list of coordinates. This facilitates using and avoids unneeded conversion by user.
	B_coordinates - same as for A.
	B_indices - same as for A.
	d_threshold - distance threshold for contacts
	returns - list of tuples (A_ind, B_ind, distance)
	"""
	bond_threshold=1.7

	A_tree=cKDTree(A_coordinates)
	B_tree=cKDTree(B_coordinates)

	# print A_indices
	# print B_indices
	#We can calculate neighbors but it is not needed since it is automatically done by distance matrix calculation
	#neighbors=A_tree.query_ball_tree(B_tree,d_threshold)
	#print neighbors #it is a list of lists

	dok_dist=A_tree.sparse_distance_matrix(B_tree,d_threshold)
	# Now we have a sparse distance matrix in a dictionary of keys format (see scipy.sparse.dok_matrix)

	# In case we need to exclude bonded let's do it
	if(exclude_bonded):
		bonded=A_tree.query_pairs(bond_threshold)
		bonded_upto_1_4=bonded_1_4(bonded)
		# print bonded_upto_1_4
		#Bonded is a set
		#Now we need also to construct 1-3, 1-4 bonded sets

	#Now let's iterate through the matrix and form an output list of tuples

	contact_list=list()
	# print dok_dist
	for key in iter(dok_dist.keys()):
		if(exclude_bonded):
			if((key not in bonded_upto_1_4) and ((key[1],key[0]) not in bonded_upto_1_4)):
				if(dok_dist.get(key)>0.1):
					# print dok_dist.get(key) # self-exclusion
					if(half_matrix):
						if(key[0]<=key[1]): #only half of matrix is needed
							contact_list.append((A_indices[key[0]],B_indices[key[1]],dok_dist.get(key)))
					else:
						contact_list.append((A_indices[key[0]],B_indices[key[1]],dok_dist.get(key)))

			# else:
				# print "Excluding bonded"
		else:
			contact_list.append((A_indices[key[0]],B_indices[key[1]],dok_dist.get(key)))
	


	return contact_list


def bonded_1_4(bonding_set):
	"""This function builds a bonding set upto 1-4 neighbours (inclusive)
	Parameters
	----------
	bonding_set - is a set, returned by KDTree.query_pairs, of bonded pairs (atoms connected by bonds)
	Returns: 1-4 bonded set, a set including pairs of atoms connected by up to 3 bonds (inclusive).
	"""
	bonding_1_3_set=set()
	bonding_1_4_set=set()

#the keys are ordered i<j
# we exploit that the cycle visits every pair two times
# i<k is indeed visited twice
	for i in bonding_set:
		for k in bonding_set:
			if(i[0]==k[0]):
				if(i[1]<k[1]): bonding_1_3_set.update([(i[1],k[1])])

			if(i[0]==k[1]):
				if(i[1]<k[0]): bonding_1_3_set.update([(i[1],k[0])])

			if(i[1]==k[0]): 
				if(i[0]<k[1]): bonding_1_3_set.update([(i[0],k[1])])

			if(i[1]==k[1]):
				if(i[0]<k[0]): bonding_1_3_set.update([(i[0],k[0])])

	for i in bonding_1_3_set:
		for k in bonding_set:
			if(i[0]==k[0]):
				if(i[1]<k[1]): bonding_1_4_set.update([(i[1],k[1])])

			if(i[0]==k[1]):
				if(i[1]<k[0]): bonding_1_4_set.update([(i[1],k[0])])

			if(i[1]==k[0]): 
				if(i[0]<k[1]): bonding_1_4_set.update([(i[0],k[1])])

			if(i[1]==k[1]):
				if(i[0]<k[0]): bonding_1_4_set.update([(i[0],k[0])])

	combined_set=bonding_set | bonding_1_3_set | bonding_1_4_set
	return combined_set


def find_hbonds(A_heavy_coord,A_heavy_indices,A_hydr_coord,B_heavy_coord,B_heavy_indices,B_hydr_coord,D_A_thresh=3.5,D_H_A_ang_thresh=30,same=False):
	""" Function finds Hydrogen bonds between atoms in groups A and B.
	Parameters
	----------
	A_heavy_coord - list of tuples with (x,y,z) of heavy atoms (donors and aceptors) in group A.
	A_heavy_indices - list of indices of heavy atoms in group A atoms corresponding to the list of coordinates. This facilitates function calls and avoids unneeded conversion by user.
	A_hydr_coord - list of tuple with (x,y,z) of hydrogen atoms in group A. The bonding will be determined atomatically.
	B - everything the same as for A.
	D_A_thresh - distance threshold for H-bonds, the distance between donor and aceptor.
	D_H_A_ang_thresh - angle threshold for deviation of D-H-A angle from 180 degrees. E.g. 30 would mean that the angle should be more than 150 degrees.
	
	same - means A and B are the same
	Output
	----------
	Returns - list of tuples (A_ind, B_ind, A_type, distance, angle)
	A_type is either D (donor) or A (acceptor) and determines the bond direction. 
	"""

	D_H_bond_threshold=1.3

	A_tree=cKDTree(A_heavy_coord)
	B_tree=cKDTree(B_heavy_coord)

	if(len(A_hydr_coord)>0):
		A_hydr_tree=cKDTree(A_hydr_coord)
	if(len(B_hydr_coord)>0):
		B_hydr_tree=cKDTree(B_hydr_coord)

	#First filter is based on distance
	# Let's get the pairs with distances < D-A threshold, and save their distances to a sparse matrix
	dok_dist=A_tree.sparse_distance_matrix(B_tree,D_A_thresh)
	# Now we have a sparse distance matrix in a dictionary of keys format (see scipy.sparse.dok_matrix)

	#We need not to implement the angle criterion
	# Mathematically: we have a pair (a,b):
	# we need to determine if their exists a hydrogen bonded to atom 
	# a or b, such as the D-H-A angle criterion is satisfied.

	# We need to know which heavy atoms are bonded to which hydrogens
	# TODO: ideally we would need to make a new KDtree only for atoms that are closer than D_A_thresh
	# for now we will go a brute force way.
	if(len(A_hydr_coord)>0):
		A_hydr_list=A_tree.query_ball_tree(A_hydr_tree,D_H_bond_threshold)
	#will return a list of lists with "internal" indices of hydrogens in A_hydr_tree
	if(len(B_hydr_coord)>0):
		B_hydr_list=B_tree.query_ball_tree(B_hydr_tree,D_H_bond_threshold)

	#now comes the search for hydrogen and D-H-A angle satisfaction

	hbond_list=list()

	for key in iter(dok_dist.keys()):

		if(same):
			if(key[0]>key[1]): #only half of matrix is needed
				continue

		#suppose atom a with idex key[0] is D, and b with index key[1] is A
		#let's loop through hydrogens of a
		if(len(A_hydr_coord)>0):
			for h_ind in A_hydr_list[key[0]]:
				ang=get_D_H_A(A_heavy_coord[key[0]],A_hydr_coord[h_ind],B_heavy_coord[key[1]])
				if(ang>(180-D_H_A_ang_thresh)):
					hbond_list.append((A_heavy_indices[key[0]],B_heavy_indices[key[1]],'D',dok_dist.get(key),ang))

		#suppose atom a with idex key[0] is A, and b with index key[1] is D
		#let's loop through hydrogens of b
		if(len(B_hydr_coord)>0):
			for h_ind in B_hydr_list[key[1]]:
				ang=get_D_H_A(B_heavy_coord[key[1]],B_hydr_coord[h_ind],A_heavy_coord[key[0]])
				if(ang>(180-D_H_A_ang_thresh)):
					hbond_list.append((A_heavy_indices[key[0]],B_heavy_indices[key[1]],'A',dok_dist.get(key),ang))


	return hbond_list

def get_D_H_A(D_coord,H_coord,A_coord):
	""" Function calculates the D-H-A angle
	Parameters
	----------
	Input are just three tuples with coordinates.
	Output is angle in degrees
	"""
	D=np.array(D_coord)
	H=np.array(H_coord)
	A=np.array(A_coord)
	HD=D-H
	HA=A-H
	return (np.arccos(np.dot(HD,HA)/(np.linalg.norm(HD)*np.linalg.norm(HA)))/3.14159*180)



def find_watermed_bonds2(A_heavy_coord,A_heavy_indices,A_hydr_coord,B_heavy_coord,B_heavy_indices,B_hydr_coord,W_heavy_coord,W_heavy_indices,W_hydr_coord,D_A_thresh=3.5,D_H_A_ang_thresh=30):
	""" Function finds water mediated bonds between atoms in groups A and B.
	Variant 2, subject to testing, and speed testing.
	Parameters
	----------
	A_heavy_coord - list of tuples with (x,y,z) of heavy atoms (donors and aceptors) in group A.
	A_heavy_indices - list of indices of heavy atoms in group A atoms corresponding to the list of coordinates. This facilitates function calls and avoids unneeded conversion by user.
	A_hydr_coord - list of tuple with (x,y,z) of hydrogen atoms in group A. The bonding will be determined atomatically.
	B - everything the same as for A.
	W - everything for water
	D_A_thresh - distance threshold for H-bonds, the distance between donor and aceptor.
	D_H_A_ang_thresh - angle threshold for deviation of D-H-A angle from 180 degrees. E.g. 30 would mean that the angle should be more than 150 degrees.
	
	Output
	----------
	Returns - list of tuples (A_ind, B_ind, W_ind)
	"""

	#Let's rely on find_hbonds function
	A_W_list=find_hbonds(A_heavy_coord,A_heavy_indices,A_hydr_coord,W_heavy_coord,W_heavy_indices,W_hydr_coord,D_A_thresh,D_H_A_ang_thresh)
	#Now we have a list of A W atoms, assuming their are a lot of waters
	# let's now construct a list of W atoms, that are bonded and then do hbond search between them
	# and B
	a_list,water_list,c,d,e=zip(*A_W_list)
	#Now let's get rid of duplicates through this hack
	WR_heavy_indices=list(OrderedDict.fromkeys(water_list))
	#Now we will need to construct a reduced water coordinates set corresponding to this list
	WR_heavy_coord=list()
	W_ind_dict=dict((j,i) for i,j in enumerate(W_heavy_indices))
	for i in WR_heavy_indices:
		WR_heavy_coord.append(W_heavy_coord[W_ind_dict[i]])

	#Now we need to run find hbonds once more
	B_W_list=find_hbonds(B_heavy_coord,B_heavy_indices,B_hydr_coord,WR_heavy_coord,WR_heavy_indices,W_hydr_coord,D_A_thresh,D_H_A_ang_thresh)

	#Now let's construct a final list
	water_med_list=list()
	b_list,water2_list,c,d,e=zip(*B_W_list)
	#we will need to get all ai given wi

	for bi,wi2 in zip(b_list,water2_list):
		for ai,wi in zip(a_list,water_list):
			if(wi==wi2):
				water_med_list.append((ai,bi,wi))

	return water_med_list



def find_watermed_bonds(A_heavy_coord,A_heavy_indices,A_hydr_coord,B_heavy_coord,B_heavy_indices,B_hydr_coord,W_heavy_coord,W_heavy_indices,W_hydr_coord,D_A_thresh=3.5,D_H_A_ang_thresh=30,same=False):
	""" Function finds water mediated bonds between atoms in groups A and B.
	A brute force version.
	Water mediated bond is a pair of hydrogen bonds connecting two atoms in group A and B with
	the same water molecule.
	Note: in case atom A has 2 HB with water molecule W, and B has 3 HB with the same molecule,
	the total number of water mediated H-bonds reported is 2*3=6
	So normally, you would like to process the output of this function to get only unique waters.
	Parameters
	----------
	A_heavy_coord - list of tuples with (x,y,z) of heavy atoms (donors and aceptors) in group A.
	A_heavy_indices - list of indices of heavy atoms in group A atoms corresponding to the list of coordinates. This facilitates function calls and avoids unneeded conversion by user.
	A_hydr_coord - list of tuple with (x,y,z) of hydrogen atoms in group A. The bonding will be determined atomatically.
	B - everything the same as for A.
	W - everything for water
	D_A_thresh - distance threshold for H-bonds, the distance between donor and aceptor.
	D_H_A_ang_thresh - angle threshold for deviation of D-H-A angle from 180 degrees. E.g. 30 would mean that the angle should be more than 150 degrees.
	
	Output
	----------
	Returns - list of tuples (A_ind, B_ind, W_ind)
	"""

	#Let's rely on find_hbonds function
	A_W_list=find_hbonds(A_heavy_coord,A_heavy_indices,A_hydr_coord,W_heavy_coord,W_heavy_indices,W_hydr_coord,D_A_thresh,D_H_A_ang_thresh)
	B_W_list=find_hbonds(B_heavy_coord,B_heavy_indices,B_hydr_coord,W_heavy_coord,W_heavy_indices,W_hydr_coord,D_A_thresh,D_H_A_ang_thresh)
	
	water_med_list=list()

	if(len(A_W_list)*len(B_W_list)==0):
		return water_med_list
	#Now brute force
	a_list,water_list,c,d,e=zip(*A_W_list)
	b_list,water2_list,c,d,e=zip(*B_W_list)

	a_bonds=zip(a_list,water_list)
	b_bonds=zip(b_list,water2_list)

	if(same):
		for i in range(len(a_bonds)):
			for j in range(i,len(b_bonds)):
				if(a_bonds[i][1]==b_bonds[j][1]):
					water_med_list.append((a_bonds[i][0],b_bonds[j][0],a_bonds[i][1]))
	else:
		for ai,wi in a_bonds:
			for bi,wi2 in b_bonds:
				if(wi==wi2):
					water_med_list.append((ai,bi,wi))

	return water_med_list

def find_ionmed_bonds(A_coordinates,A_indices,B_coordinates,B_indices,ION_coordinates,ION_indices,d_threshold=3.9,same=False):
	""" Function finds ion mediated bonds between atoms in groups A and B.
	A brute force version
	Parameters
	----------
	A_coordinates - list of tuples with (x,y,z) of atoms in group A.
	A_indices - list of indices of atoms in group A atoms corresponding to the list of coordinates. This facilitates function calls and avoids unneeded conversion by user.
	B - everything the same as for A.
	ION - everything for ions
	d_threshold - distance threshold for contacts
	Output
	----------
	Returns - list of tuples (A_ind, B_ind, ION_ind)
	"""

	#Let's rely on find_contacts function
	A_I_list=find_contacts(A_coordinates,A_indices,ION_coordinates,ION_indices,d_threshold,False)
	B_I_list=find_contacts(B_coordinates,B_indices,ION_coordinates,ION_indices,d_threshold,False)
	
	ion_med_list=list()
	
	if(len(A_I_list)*len(B_I_list)==0):
		return ion_med_list
	#Now bute force
	a_list,ion_list,c=zip(*A_I_list)
	b_list,ion2_list,c=zip(*B_I_list)

	a_bonds=zip(a_list,ion_list)
	b_bonds=zip(b_list,ion2_list)

	if(same):
		for i in range(len(a_bonds)):
			for j in range(i,len(b_bonds)):
				if(a_bonds[i][1]==b_bonds[j][1]):
					ion_med_list.append((a_bonds[i][0],b_bonds[j][0],a_bonds[i][1]))
	else:
		for ai,ii in zip(a_list,ion_list):
			for bi,ii2 in zip(b_list,ion2_list):
				if(ii==ii2):
					ion_med_list.append((ai,bi,ii))

	return ion_med_list





if __name__ == '__main__':
	
	#Understanding KDTree
	import numpy as np                                                                                                                                                                                             
	import scipy.spatial as ss                                                                                                                                                                                     
	from itertools import combinations
	from itertools import permutations
	import time
	d_threshold = 1.5
	x=[1,2,3,4,5,6,7,8,9,10]
	y=[1,2,3,4,5,6,7,8,9,10]
	points = zip(x, y)
	x2=[1,3,4,5,6,7,8,8,10,16]
	y2=[1,2,3,4,5,6,7,8,9,16.1]
	points2 = zip(x, y)
	start = time.clock()
	A_tree=KDTree(points)
	B_tree=KDTree(points2)
	neighbors=A_tree.query_ball_tree(B_tree,d_threshold)
	# print neighbors #it is a list of lists
	distances=A_tree.sparse_distance_matrix(B_tree,d_threshold)
	# print distances
	A_tree=cKDTree(points)
	B_tree=cKDTree(points2)
	# print points[0]
	cneighbors=A_tree.query_ball_tree(B_tree,d_threshold)
	# print neighbors #it is a list of lists
	cdistances=A_tree.sparse_distance_matrix(B_tree,d_threshold)
	# print distances
	# print neighbors
	# print cneighbors
	#print distances
	# if(neighbors==cneighbors):
		# print "Yes!!!"

	# di=distances.items()
	# cdi=cdistances.items()
	# if(cdi==di):
		# print "Yes!!!"
	# for i in distances.iterkeys():
		# print i
	# print distances
	# print cdistances
	D=(-2,0,0)
	H=(0,0,0)
	A=(-1,0.1,0)
	print(get_D_H_A(D,H,A))

	elapsed = (time.clock() - start)
	print("Time elapsed: %f s" % elapsed)
