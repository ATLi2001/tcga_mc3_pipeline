#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 18:35:47 2018

@author: Tianyan
"""

from numpy	  import *
from numpy.random import *
from tssb		import *
from util		import *
from util2 import *

from ete2 import *

from subprocess import call
import sys
import os

import csv #need to write out to csv

from r_medicc_ce import *

#node counters
ctr=0
ctr2=0

#number of ssms
num_ssms = 0

#fout is now the file name
def print_top_trees(tree_archive,fout,outloc,k):
    global ctr;
    global ctr2
    tree_reader = TreeReader(tree_archive)
    
    #files for depth, numbranches, and ce
    depthFile = open(os.path.join(outloc,os.path.basename(fout[:-4])+"_depth.txt"), 'w+')
    numbranchesFile = open(os.path.join(outloc,os.path.basename(fout[:-4])+"_numbranches.txt"), 'w+')
    ceFile = open(os.path.join(outloc,os.path.basename(fout[:-4])+"_CE.txt"), 'w+')
    
    #files for max muts, num nodes, and total muts
    maxmutsFile = open(os.path.join(outloc,os.path.basename(fout[:-4])+"_maxmuts.txt"), 'w+')
    numnodesFile = open(os.path.join(outloc,os.path.basename(fout[:-4])+"_numnodes.txt"), 'w+')
    totalmutsFile = open(os.path.join(outloc,os.path.basename(fout[:-4])+"_totalmuts.txt"), 'w+')

    #file for trunk proportion
    trunkpropFile = open(os.path.join(outloc,os.path.basename(fout[:-4])+"_trunkprop.txt"), 'w+')
    cnvtrunkpropFile = open(os.path.join(outloc,os.path.basename(fout[:-4])+"_trunkprop_cnv.txt"), 'w+')
    
    for idx, (tidx, llh, tree) in enumerate(tree_reader.load_trees_and_metadata(k)):
        ctr=0
        ctr2=0
        remove_empty_nodes(tree.root, None)
        
#        print(idx)
        
        #number of ssms in the tree
        global num_ssms
        num_ssms = get_num_ssms(tree.root, None, tree)
#        print(num_ssms)
        
        
        
        #current before matrix
        currBtMat = get_before_Matrix(tree)
        
        #current equal matrix
        currEqMat = get_equal_Matrix(tree)
        
        #distance matrix for the MEDICC analysis
        dist = get_distance_Matrix(tree)
        
        #first time is just the unchanged matrix
        if idx == 0:
            before_result = add_blank_Matrix(currBtMat, 0)
            equal_result = add_blank_Matrix(currEqMat, 0)
            gdict = get_gene_dict(tree)
        
                        
        #rest of the times add the current tree to the existing
        else:
            before_result = add_Matrix(currBtMat, before_result)
            equal_result = add_Matrix(currEqMat, equal_result)


        #write the depth, numbranches, and ce
        depthFile.write(str(get_depth(tree)) + '\n')
        numbranchesFile.write(str(get_num_branches(tree)) + '\n')
        ceFile.write(str(medicc_ce(dist, idx))[4:])
        #write the maxmuts, numnodes, and totalmuts
        maxmutsFile.write(str(get_max_muts(tree)) + '\n')
        numnodesFile.write(str(get_num_subclones(tree.root, None, tree)) + '\n')
        totalmutsFile.write(str(get_num_ssms(tree.root, None, tree)) + '\n')
        #write the trunkprop
        trunkpropFile.write(str(get_trunk_prop(tree)) + '\n')
        cnvtrunkpropFile.write(str(get_cnv_trunk_prop(tree)) + '\n')
        
    #after matrix is just before matrix transposed
    after_result = transpose_Matrix(before_result)
        
    #write to csv
    before_result_csv = add_matrix_row_col_names(before_result, gdict)
    after_result_csv = add_matrix_row_col_names(after_result, gdict)
    equal_result_csv = add_matrix_row_col_names(equal_result, gdict)
    
    #write after matrix results to csv
    #use [:-4] to get rid of .txt extension in fout
    #use basename to prevent fout from containing the input directory
    with open(os.path.join(outloc,os.path.basename(fout[:-4])+"_after.csv"), 'w') as afterFile:
        writer = csv.writer(afterFile)
        writer.writerows(after_result_csv)
        
    #write before matrix results to csv
    with open(os.path.join(outloc,os.path.basename(fout[:-4])+"_before.csv"), 'w') as beforeFile:
        writer = csv.writer(beforeFile)
        writer.writerows(before_result_csv)
    
    #write equal matrix results to csv
    with open(os.path.join(outloc,os.path.basename(fout[:-4])+"_equal.csv"), 'w') as equalFile:
        writer = csv.writer(equalFile)
        writer.writerows(equal_result_csv)
        
    
    #when combining for ordering matrix, before will be represented by negatives
    before_result = scalar_matrix_division(before_result, -1)
    #convert to ints
    for i in range(len(before_result)):
        for j in range(len(before_result[0])):
            before_result[i][j] = int(before_result[i][j])
    
    #combine for the overall ordering matrix
    ordering_matrix = add_Matrix(before_result, after_result)
    ordering_matrix = scalar_matrix_division(ordering_matrix, k)
    ordering_matrix_csv = add_matrix_row_col_names(ordering_matrix, gdict)
    
    #wrtie ordering matrix to csv
    with open(os.path.join(outloc,os.path.basename(fout[:-4])+"_ordering_matrix")+".csv", 'w') as outFile:
        writer = csv.writer(outFile)
        writer.writerows(ordering_matrix_csv)
    
    #close all files
    depthFile.close()
    numbranchesFile.close()
    ceFile.close()
    maxmutsFile.close()
    numnodesFile.close()
    totalmutsFile.close()
    trunkpropFile.close()
    cnvtrunkpropFile.close()
    tree_reader.close()
    
    

def get_before_Matrix(tssb):
	
    t = Tree();t.name='0'
    
    #total number of genes in the tree
    total_genes = get_num_genes(tssb.root, None, t)
    
#    print(total_genes)
    
    #matrix to hold the before pairs
    before_count_Matrix = [[0 for x in range(total_genes)] for y in range(total_genes)]
    
    #fill the matrix
    fill_before_Matrix(tssb.root, None, t, before_count_Matrix)
        
#    print(before_count_Matrix)
    
    return before_count_Matrix

def get_equal_Matrix(tssb):
    
    t = Tree();t.name='0'
    
    #total number of genes in the tree
    total_genes = get_num_genes(tssb.root, None, t)
    
    #matrix to hold the equal pairs
    equal_count_Matrix = [[0 for x in range(total_genes)] for y in range(total_genes)]
    
    #fill the matrix
    fill_equal_Matrix(tssb.root, None, t, equal_count_Matrix)
    
#    print(equal_count_Matrix)
    
    return equal_count_Matrix

#before than matrix
#note that this function does not return anything, it fills in btMat
#note that this matrix wouldn't include non-descendant counts
#i.e. if the structure were
#   
#     /---1---3
#   0
#     \---2---4
# the genes in 4 would not be included in the descendants of 1 
def fill_before_Matrix(node, parent,tree, btMat):
    global ctr;
    node_name  = ctr ; ctr+=1;
    
    #genes in the node
    genes = node['node'].get_data()
    
    #descendant genes
    desc = get_gene_descendants(node, parent, tree)
    
    global num_ssms
    
    #loop through the descendants and change the names to numbers
    for ids in range(len(desc)):
        desc[ids] = convert_gene_name_to_int(desc[ids], num_ssms)
    
    #change names to numbers
    gene_ids = []
    for g in range(len(genes)):
        gene_ids.append(convert_gene_name_to_int(genes[g].id, num_ssms))
    
    #need to have genes to work with
    if(len(genes))>0:
        #loop through the genes 
        for i in range(len(genes)):
            #loop through the descendants
            for d in range(len(desc)):
                #compare each of this level's genes to all lower level
                btMat[gene_ids[i]][desc[d]] = 1
            
    #recursion
    for child in node['children']:
		name_string = str(ctr)#+'('+str(len(child['node'].get_data()))+')'
		fill_before_Matrix(child, node_name,tree.add_child(name=name_string), btMat)
    

#equal level matrix
#i.e. if the structure were
#   
#     /---1---3
#   0
#     \---2---4
# the genes in 1 one would only be equal to the others within it (not 2)
def fill_equal_Matrix(node, parent, tree, eqMat):
    global ctr;
    node_name  = ctr ; ctr+=1;
    
    #genes in the node
    genes = node['node'].get_data()
    
    #number of ssms in the tree
    global num_ssms
    
    #change names to numbers
    gene_ids = []
    for g in range(len(genes)):
        gene_ids.append(convert_gene_name_to_int(genes[g].id, num_ssms))
        
    #always end up going through the node internally
    #need genes to work with
    if(len(genes))>0:
        #loop through genes twice
        for i in range(len(genes)):
            for j in range(len(genes)):
                #each pair within is equal level
                #make sure that both i,j and j,i are set to 1
                eqMat[gene_ids[i]][gene_ids[j]] = 1
                eqMat[gene_ids[j]][gene_ids[i]] = 1    

    #recursion
    for child in node['children']:
        name_string = str(ctr)
        fill_equal_Matrix(child, node,tree.add_child(name=name_string), eqMat)
    


############## additional statistics to get #################################
#DO NOT include subclone 0, which has no mutations, in any statistics


#get the depth of a tree
def get_depth(tree):
    return get_depth_helper(tree.root, None, tree)
def get_depth_helper(node, parent, tree):
    depth = 0
    i=0
    
    genes = node['node'].get_data()
    
    for child in node['children']:
        #only start doing depth if the node is not empty
        if(len(genes)>0):
            #make sure only add one depth for all the children
            if(i==0):
                depth+=1    
                i+=1
        depth += get_depth_helper(child, node, tree)
        
    return depth

#get the number of branches in a tree
def get_num_branches(tree):
    return get_num_branches_helper(tree.root, None, tree)
def get_num_branches_helper(node, parent, tree):
    num_branches = 0
    genes = node['node'].get_data()
    for child in node['children']:
        #only start counting when the current node is not empty
        if(len(genes)>0):
            num_branches += 1
        num_branches += get_num_branches_helper(child, node, tree)
        
    return num_branches

def get_distance_Matrix(tssb):
    t = Tree();t.name='0'
    #number of subclones
    num_subclones = get_num_subclones(tssb.root, None, t)
    #blank distance matrix
    dist_mat = [[0 for x in range(num_subclones)] for y in range(num_subclones)]
    
    #fill the matrix
    fill_distance_Matrix(tssb.root, None, t, dist_mat)
        
    return dist_mat

def fill_distance_Matrix(node, parent, tree, distMat):
    node_gene_dict = get_node_genes(node, parent, tree, [])
    
    for i in range(len(distMat)):
        for j in range(len(distMat[0])):
            #got rid of subclone 0; thus need to use i+1 and j+1
            distMat[i][j]=get_uncommon(node_gene_dict.get(i+1), node_gene_dict.get(j+1))
    
    

#gets all genes that have been mutated in a node
#includes all mutations of the parent
def get_node_genes(node, parent, tree, parent_gene_set):
    global ctr2;
    node_name  = ctr2 ; ctr2+=1;
    
    global num_ssms
    
    #output is a dictionary with key = nodename and values = genes
    out = {}
    
    #all data
    genes = node['node'].get_data()
    
    #set containing the genes in this subclone
    gset = []
    
    #need genes to work with
    if len(genes) > 0:
        #for each gene
        for g in range(len(genes)): 
            #values are the genes
            gnum = convert_gene_name_to_int(genes[g].id, num_ssms)
            gset.append(gnum)
    
    #value is node's mutated genes + parent's
    total_gene_set = gset + parent_gene_set
    #only include if this node has genes
    if(len(genes)>0):
        out[node_name] = total_gene_set
    
    #recursively go through the tree
    for child in node['children']:
        name_string = str(ctr2)
        #use update to continuously edit the original dictionary
        out.update(get_node_genes(child, node_name, tree.add_child(name=name_string),total_gene_set))
    
    #after all recursions, should be all genes in the dictionary
    return out

############## more statistics to get part 2 #################################

def get_max_muts(tree):
    return max(get_max_muts_helper(tree.root, None, tree))
def get_max_muts_helper(node, parent, tree):
    num_muts = []
    
    #all data
    genes = node['node'].get_data()
    
    #number of genes that are ssms
    to_add = len(genes)
    
    #if the node has cnv, don't count it
    for g in range(len(genes)):
        if(isinstance(genes[g].id, basestring)):
            if(genes[g].id[0]=='c'):
                to_add -= 1
    
    
        
    global num_ssms
        
    #number of genes in this node
    num_muts.append(to_add)
    
    #recursively go through the tree
    for child in node['children']:
		num_muts.extend(get_max_muts_helper(child, node, tree))

    #after all recursions, should be total number
    return num_muts

#proportion of mutations in the trunk of tree
def get_trunk_prop(tree):
    return get_trunk_prop_helper(tree.root, None, tree)
def get_trunk_prop_helper(node, parent, tree):
    genes = node['node'].get_data()
    
    #number of genes that are ssms
    to_add = len(genes)
    
    #if the node has cnv, don't count it
    for g in range(len(genes)):
        if(isinstance(genes[g].id, basestring)):
            if(genes[g].id[0]=='c'):
                to_add -= 1
        
    
    if(len(genes)>0):
        #use float to prevent int division
        return float(to_add)/get_num_ssms(tree.root, None, tree)
    else:
        for child in node['children']:
            return get_trunk_prop_helper(child, node, tree)

#proportion of cnvs in the trunk of the tree
def get_cnv_trunk_prop(tree):
    return get_cnv_trunk_prop_helper(tree.root,None,tree)
def get_cnv_trunk_prop_helper(node, parent, tree):
    genes = node['node'].get_data()
    
    num_cnvs = 0
    for g in range(len(genes)):
        if(isinstance(genes[g].id, basestring)):
            if(genes[g].id[0]=='c'):
                num_cnvs += 1
    
    if(len(genes) > 0):
        #use float to prevent int division
        if(get_num_cnvs(tree.root, None, tree)==0):
            return "NA"
        else:
            return float(num_cnvs)/get_num_cnvs(tree.root, None, tree)
    else:
        for child in node['children']:
            return get_cnv_trunk_prop_helper(child, node, tree)

#number of cnvs in the tree
def get_num_cnvs(node, parent, tree):
    #all data
    genes = node['node'].get_data()
    
    #number of genes in this node
    num_cnvs = 0
    for g in range(len(genes)):
        if(isinstance(genes[g].id, basestring)):
            if(genes[g].id[0]=='c'):
                num_cnvs += 1
    
    #recursively go through the tree
    for child in node['children']:
		num_cnvs += get_num_cnvs(child, node, tree)
    
    #after all recursions, should be total number
    return num_cnvs

############## miscellaneous helper functions ################################

#get the total number of genes (includes cnvs) in the tree
def get_num_genes(node, parent, tree):
    #all data
    genes = node['node'].get_data()
    
    #number of genes in this node
    num_genes = 0
    num_genes += len(genes)
    
    #recursively go through the tree
    for child in node['children']:
		num_genes += get_num_genes(child, node, tree)
    
    #after all recursions, should be total number
    return num_genes

#get the total number of ssms in the tree
def get_num_ssms(node,parent,tree):
    #all data
    genes = node['node'].get_data()
    
    #number of ssms in this node
    num_genes = 0
    to_add = len(genes)
    
    #if the node has cnv, don't count it 
    for g in range(len(genes)):
        if(isinstance(genes[g].id, basestring)):
            if(genes[g].id[0]=='c'):
                to_add -= 1
    
    num_genes += to_add
    
    #recursively go through the tree
    for child in node['children']:
		num_genes += get_num_ssms(child, node, tree)
    
    #after all recursions, should be total number
    return num_genes

#has one parameter, calls the helper
def get_gene_dict(tssb):
    t = Tree();t.name='0'
    return get_gene_dict_helper(tssb.root, None, t)

#create a dictionary to map from the gene id to the name
def get_gene_dict_helper(node, parent, tree):
    global num_ssms
    
    #gene dictionary
    gene_dict = {}
    
    #all data
    genes = node['node'].get_data()
    
    #need genes to work with
    if len(genes) > 0:
        #for each gene
        for g in range(len(genes)):
            #key is the gene id number
            gnum = convert_gene_name_to_int(genes[g].id, num_ssms)
            #value is the name of the gene
            gene_dict[gnum] = genes[g].name
    
    #recursively go through the tree
    for child in node['children']:
        #use update to continuously edit the original dictionary
		gene_dict.update(get_gene_dict_helper(child, node, tree))
    
    #after all recursions, should be all genes in the dictionary
    return gene_dict

#get a list of the genes below a node
#direct descendants only
def get_gene_descendants(node, parent, tree):
    #output list
    out = []
    
    #recursively go through the tree
    for child in node['children']:
        #genes of the child
        genes = child['node'].get_data()
        
        #has genes
        if len(genes)>0:
            #for each gene, add it to out
            for g in range(len(genes)):
                out.append(genes[g].id)
                
        #use extend to prevent lists in lists    
        out.extend(get_gene_descendants(child, node, tree))
    
    return remove_duplicates(out)


#gets the number of subclones in a tree
def get_num_subclones(node, parent, tree):
    num_subclones = 0
    
    genes = node['node'].get_data()
    
    #don't include emtpy nodes in count
    if(len(genes)>0):
        num_subclones+=1
    
    #recursively go through the tree
    for child in node['children']:
		num_subclones += get_num_subclones(child, node, tree)
    
    #after all recursions, should be total number
    return num_subclones

#gets the number of uncommon elements between 2 lists
def get_uncommon(X, Y):
    total = len(X) + len(Y)
    #length of the intersection
    common = len(set(X).intersection(Y))
    
    #principle of inclusion exclusion; don't want to count common at all
    return total - common*2

#convert string a to int
#if a isn't already int, it is of the form s10
#n is the number of ssms in the tree
def convert_gene_name_to_int(a, n):
    if type(a) == int:
        return a
    #numbering for cnvs will start after last ssm number
    #i.e. if ssms go up from s0 to s11, c0 will become 12
    elif a[0]=="c":
        return int(n) + int(a[1:])
    else:
        return int(a[1:])
        

#note that X, Y are of same size
def add_Matrix(X, Y):
    #resultant matrix
    resultant = [[0 for x in range(len(X))] for y in range(len(X[0]))]
    
    # iterate through rows
    for i in range(len(X)):
        # iterate through columns
        for j in range(len(X[0])):
            resultant[i][j] = X[i][j] + Y[i][j]
    
    return resultant

#take in an int blank and add it to a matrix
def add_blank_Matrix(X, blank):
    #create a matrix of same size as X
    blank_Matrix = [[blank for x in range(len(X))] for y in range(len(X[0]))]
    #add the two matrices
    return add_Matrix(X, blank_Matrix)

#transpose a matrix (switch the rows and columns)
def transpose_Matrix(X):
    trans = [[X[j][i] for j in range(len(X))] for i in range(len(X[0]))]
    return trans

#adds row and col names to the output matrix using the gene dictionary
def add_matrix_row_col_names(X, d):
    gnames = d.values()
    
    new_row_len = len(X) + 1
    new_col_len = len(X[0]) + 1
    
    newMat = [[0 for x in range(new_row_len)] for y in range(new_col_len)]
    
    for i in range(len(newMat)):
        for j in range(len(newMat[0])):
            if i == 0 and j == 0:
                newMat[i][j] = ''
            elif i == 0 and j != 0:
                newMat[i][j] = gnames[j-1]
            elif j == 0 and i != 0:
                newMat[i][j] = gnames[i-1]
            else:
                newMat[i][j] = X[i-1][j-1]
    
    return newMat

#divide matrix X by s
def scalar_matrix_division(X,s):
    Y = X
    for i in range(len(X)):
        for j in range(len(X[0])):
            if(not isinstance(X[i][j], basestring)):
                Y[i][j] = X[i][j] / (s * 1.0)
    
    return Y


def remove_duplicates(l):
    return list(set(l))
