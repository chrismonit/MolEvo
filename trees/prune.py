#!/usr/bin/python

# Author: Christopher Monit

import argparse
from Bio import Phylo


def readSequenceNames(sequenceNamesFilePath):
    """ Read the names of the taxa from the given file path """
    f = open(sequenceNamesFilePath)
    sequenceNames = [ line.rstrip() for line in f.readlines() ]
    f.close()
    return sequenceNames 




def pruneTreeSequenceNames(tree, sequenceNames):
    """ Prune tree of every terminal sequence whose name is given in sequenceNames """
    
    for tip in tree.get_terminals():
        if sequenceNames.__contains__(tip.name):
            tree.prune(tip)
    return tree


def pruneTreeSequencePrefixes(tree, sequencePrefixes):
    """ Prune tree of every terminal sequence whose name is prefixed by a label given in seqeuencePrefixes"""
    
    for prefix in sequencePrefixes: # can this be made more efficient?? are two for loops necessary?  
        for tip in tree.get_terminals():
            if tip.name.startswith(prefix):
                tree.prune(tip) 
    return tree




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Prune one of more Newick trees of specified taxa.')
    parser.add_argument( 'treePath' )
    parser.add_argument( 'newTreePath' )
    parser.add_argument( '-p', nargs='+', help='prefixes for sequence names to be pruned (e.g. \'hu.\' for human)' )
    parser.add_argument( '-f', help="file containing names of taxa to be pruned from tree (in a single column with no header)")
    parser.add_argument( '-t', nargs='+', help='names of individual taxa to be pruned' )
    args = parser.parse_args() 

    
    if args.p != None and  args.f != None and args.t != None:
        print "ERROR no sequences provided. Specify individual sequence names with -f or -t, or prefixes for groups of sequences with -p"
        exit()

    if args.p != None and ( args.f != None or args.t != None ):
        print "\n\nERROR cannot prune sequences specified both individually (-f or -t) AND as prefixes (-p)\n\n"
        exit()

   
    if args.p != None: # we are specifiying sequences using prefixes 
        
        trees = list(Phylo.parse(args.treePath, 'newick')) # NB may only be one tree
        for tree in trees:
            tree = pruneTreeSequencePrefixes(tree, args.p)
        
        Phylo.write( trees, args.newTreePath, 'newick' )
        

    else: # we are using full sequence names (from file and/or command line)

        sequenceNames = []
        if args.f != None: sequenceNames += readSequenceNames(args.f)
        if args.t != None: sequenceNames += args.t

        trees = list(Phylo.parse(args.treePath, 'newick')) # NB may only be one tree
        for tree in trees:
            tree = pruneTreeSequenceNames(tree, sequenceNames)
        
        Phylo.write( trees, args.newTreePath, 'newick' )



