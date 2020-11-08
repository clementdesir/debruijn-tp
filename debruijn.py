#!/usr/local/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""



import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
#from collections import Counter

__author__ = "Clément Désir"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Clément Désir"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Clément Désir"
__email__ = "clementdesir@gmail.com"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()

def read_fastq(fasta_file):
    cpt = 0
    with open(fasta_file, "r") as file:
        for line in file:
            if cpt % 2 != 0:
                sequence = line.strip('\n')
                yield(sequence) 
            cpt += 1

def cut_kmer(seq, k):
    kmer = ""
    for i in range(len(seq)+1-k):
        kmer = seq[i:i+k]
        yield(kmer)

def build_kmer_dict(fasta_file, kmer):
    kmer_dict = dict()
    kmer_set = []
    for line in read_fastq(fasta_file):
        kmer_set.append(cut_kmer(line,kmer))
    cpt = 0
    for x in kmer_set:
        if x in kmer_dict:
            kmer_dict[x] += 1
        else:
            kmer_dict[x] = 1
    return kmer_dict

def build_graph(kmer_dict):
    graph = nx.DiGraph()
    for k in kmer_dict.keys():
        prefixe = k[:-1]
        suffixe = k[1:]
        graph.add_edge(prefixe, suffixe, weight = kmer_dict[k])
    return graph

def get_starting_nodes(graph):
    entry_node = []
    for i in graph.nodes():
        if i not in graph.predecessors(i):
            entry_node.append(i)
    return entry_node 

def get_sink_nodes(graph):
    sink_node = []
    for i in graph.nodes():
        if i not in graph.successors(i):
            sink_node.append(i)
    return sink_node 

def get_contigs(graph, start_nodes, end_nodes):
    contigs = []
    for start_node in start_nodes:
        for end_node in end_nodes:
            path = nx.shortest_path(graph, start_node, end_node)
            contigs.append(path, len(path))
    return contigs

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def save_contigs(contig_list, output_file):
    with open(output_file, "w") as final_file:
        for contig in contig_list:
            final_file.write(">contig_Numero{} len = {}".format(contig_list.index(),contig))
            final_file.write(fill(contig[0]))

def std(list):
    return statistics.stdev(list)

def path_average_weight(graph, path):
    av_w = 0
    for node in path:
        av_w += node["weight"]
    return av_w/(len(path)-1)

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    for path in path_list:
        if delete_entry_node == True and delete_sink_node == True:
            graph.remove_nodes_from(path)
        elif delete_entry_node == False and delete_sink_node == True:
            graph.remove_nodes_from(path[:len(path)-1])
        elif delete_entry_node == False and delete_sink_node == False:
            graph.remove_nodes_from(path[1:])    
        else:
            graph.remove_nodes_from(path[1:len(path)-1])
    return graph


def select_best_path(graph, path_list, path_length, weight_average_list, delete_entry_node=False, delete_sink_node=False):   
    best_weight_av_indexes = []     
    for weight in weight_average_list:
        if weight == max(weight_average_list):
            best_weight_av_indexes.append(weight.index())
    best_length_indexes = []
    for length in path_length:
        if length.index() in best_weight_av_indexes:
            best_length_indexes.append(length.index())
    best_path_indexes = []
    for path in path_list:
        if path.index() == max(best_length_indexes).index():
            best_path_indexes.append(path.index())
    best_path_index = random.choice(best_path_indexes)
    graph.remove_paths(graph,path_list[:best_path_index] + path_list[best_path_index + 1:],delete_entry_node,delete_sink_node)
    return graph

def solve_bubble(graph, ancestor_node, descendant_node):
    path_list = nx.all_simple_paths(graph, ancestor_node, descendant_node)
    path_lengths = []
    avg_weight_list = []
    for path in path_list:
        path_lengths.append(len(path))
        avg_weight_list.append(path_average_weight(path))
    graph.select_best_path(graph, path_list, path_lengths,avg_weight_list)
    return graph

#def simplify_bubble(graph):



        





        





#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    for line in read_fastq(args.fastq_file):
        print(line)
    print("\n")
    print("kmer size: ", args.kmer_size)
    print(build_kmer_dict(args.fastq_file, args.kmer_size))
    #print("number of kmers : ", len(build_kmer_dict(args.fastq_file, args.kmer_size)))
    #graph = build_graph(build_kmer_dict(args.fastq_file, args.kmer_size))
    

    #starting_nodes = get_starting_nodes(graph)
    #exit_nodes = get_sink_nodes(graph)
    #print(starting_nodes)
    #print(exit_nodes)

if __name__ == '__main__':
    main()
