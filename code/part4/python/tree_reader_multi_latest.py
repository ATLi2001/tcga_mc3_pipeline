#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 16:06:52 2018

@author: Tianyan
"""

#!/usr/bin/env python2

#new printo
import printo2_multi_latest
#old printo
import printo

import argparse


#call all functions from printo2 so it does the work
#this program takes in the zip and an output directory
def start_new_run(trees, ssm, out_dir, k):
#    state = {}
#    
#    state['trees'] = trees
#    state['out_dir'] = out_dir
#    
    printo2_multi_latest.print_top_trees(trees, ssm, out_dir, k)
    printo.print_top_trees(trees,ssm,out_dir, k)

def parse_args():
    parser = argparse.ArgumentParser(
        description='Get tree information from a zip file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
   
    parser.add_argument('-t', '--trees', dest='trees',
        help='Zip file for trees')
    parser.add_argument('-s', '--ssm', dest='ssm',
        help='SSM file so output files can have patient id')
    parser.add_argument('-O', '--output-dir', dest='output_dir',
        help='Path to directory for output files')
    parser.add_argument('-k', '--num-trees', dest='k',type=int,
        help='Number of trees to be read')
    
    args = parser.parse_args()
    return args

def run():
    
    args = parse_args()
    
    start_new_run(
                args.trees,
                args.ssm,
                args.output_dir,
                args.k
            )
              

if __name__ == "__main__":
     run()
