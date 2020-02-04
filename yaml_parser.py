# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 18:19:51 2019

@author: kenelm
"""

from Bio.SeqIO.FastaIO import SimpleFastaParser
import argparse
import itertools
import re
import os

def check_input_file(parser, file_name):
    full_file_name = os.path.abspath(file_name)
    if not os.path.exists(full_file_name):
        parser.error("The file {} does not exist!".format(file_name))
    else:
        return full_file_name

def make_arguments_parser():
    parser = argparse.ArgumentParser(
        description='Filter Partis *.yaml output.')
        
    parser.add_argument('-i', dest = 'input_file', 
                        help='Input FASTA', 
                        required=True,
                        type=lambda x: check_input_file(parser, x))

    parser.add_argument('-y', dest = 'yaml_file', 
                        help='Input *.yaml', 
                        required=True,
                        type=lambda x: check_input_file(parser, x))

    parser.add_argument('-o', dest = 'output_file', 
                        help='Output file name (without extension)', 
                        required=True,
                        type=str,
                        default='0')
    
    return parser.parse_args()

if __name__ == "__main__":
    
    args = make_arguments_parser()
    
    with open(args.yaml_file,'r') as yfile:
        yfile.readline()
        yfile.readline()
        clone_description = yfile.readline()
    
    clone_description = clone_description.split("    ")    
    tags = {current_row.split(":")[0] for current_row in clone_description if re.search('[a-zA-Z]', current_row.split(":")[0])}
#    tags = list( itertools.chain.from_iterable(clone_description) )
    
    print("Processing is started.")
    print("Number of tags is {}.".format(len(tags)))
    print(tags)
    
    counter = 0
    with open(args.input_file,'r') as handle_input, open(args.output_file,'w') as handle_output:
        for title, seq in SimpleFastaParser(handle_input):                
            if any([current_tag in title for current_tag in tags]):
                counter += 1
                handle_output.write(">%s\n%s\n" % (title, seq))
    
    print("Processing is finished.")
    print("Number of tags is {}.".format(counter))    
