#!/usr/bin/python

from Bio import Phylo
import argparse
import os

description = """
Labels an internal node of a Newick tree, specified by its descendent taxa.
This is especially useful for PAML models (by Ziheng Yang, UCL), where branches of interest are labelled with "#1" etc.

Taxa whose most recent common ancestor you want to label can be specified by name prefix (-p), a file listing them (-f) or listing on the command line (-t). -f and -t can be used in conjunction, but neither can be used with -p.

For example:

$ cat tree.tre 
((a.taxon:0.01,a.taxon:0.2):0.2,(b.taxon:0.2,c.taxon:0.2):0.2);
$ label_clades.py tree.tre out.tre -p a -l \#1
$ cat out.tre 
((a.taxon:0.01000,a.taxon:0.20000)#1:0.20000,(b.taxon:0.20000,c.taxon:0.20000):0.20000):1.00000;

Note that "\#1" is necessary, so the shell does not mistake "#" for a comment. 
Alternatively:

$ cat taxa.txt
b.taxon
c.taxon
$ label_clades.py tree.tre out.tre -f taxa.txt -l \#1
$ cat out.tre 
((a.taxon:0.01000,a.taxon:0.20000):0.20000,(b.taxon:0.20000,c.taxon:0.20000)#1:0.20000):1.00000;

Or,

$ label_clades.py tree.tre out.tre -t b.taxon c.taxon -l \#1
$ cat out.tre 
((a.taxon:0.01000,a.taxon:0.20000):0.20000,(b.taxon:0.20000,c.taxon:0.20000)#1:0.20000):1.00000;

The input tree file can contain multiple trees, which will all be labelled in the output. If it already exists, output is appended to the speficied output file.

Please get in touch if you encounter bugs
Please acknowledge if you find this code useful

Christopher Monit
November 2016
"""

parser = argparse.ArgumentParser(description=description)
parser.add_argument( 'treePath' )
parser.add_argument( 'newTreePath' )
parser.add_argument( '-p', help='prefix denoting taxa whose common ancestor node is to be labelled' )
parser.add_argument( '-l', help='label to add to clade of interest' )
parser.add_argument( '-f', help="file containing names of taxa whose common ancestor node is to be labelled (in a single column with no header)")
parser.add_argument( '-t', nargs='+', help='names of individual taxa whose common ancestor node is to be labelled' )
args = parser.parse_args() 

NEWLINE = "\n"

def readSequenceNames(sequenceNamesFilePath):
        """ Read the names of the taxa from the given file path """
        f = open(sequenceNamesFilePath)
        sequenceNames = [ line.rstrip() for line in f.readlines() ]
        f.close()
        return sequenceNames 


def getNamesByPrefix(tree, prefix):
        """ Returns list of taxa whose names start with the given prefix """
        
        sequenceNames = []
        for taxon in tree.get_terminals():
	      if taxon.name.startswith(prefix):
		    sequenceNames.append(taxon.name)
        return sequenceNames

def getCladeLongNames(tree, terminalNames): 
    """ The Biopython.Phylo method Tree.find_clades has the limitation that it can't identify terminal clades with long taxon names.
        Some of the names in my HIV datasets are too long for this, so I've had to develop my own work around.
    """
    outgroupTerminals = [] # will be list of clade objects
    for clade in tree.get_terminals():
        if terminalNames.__contains__(clade.name): # NB assumes all names in terminalNames are unique (a reasonable assumption)
            outgroupTerminals.append(clade)    
        if len(outgroupTerminals) == len(terminalNames):
            break
    else:
        return None #couldn't find terminal taxa
    outgroupClade = tree.common_ancestor(outgroupTerminals)
    return outgroupClade


def labelClade(tree, taxa, labelString, confidenceValue=262019063.00): # confidence value here is a randomly generated large number, which is unlikely to appear in the taxon names as a string. Serves as a place holder only
	""" Takes tree object, labels the ancestor node with confidence value, then converts to string and replaces confidence for string label """	
        clade = getCladeLongNames(tree, taxa)
	clade.confidence = confidenceValue
        treeString = tree.format('newick')
	return treeString.replace(str(confidenceValue)+"0", labelString).rstrip() # need an extra "0" because str() method only gives representation to x.0 decimal place, but biopython works with x.00. Also the biopython method adds an undesirable newline at the end


def getNamesFromFile(filePath):
    f = open(filePath)
    lines = [ l.rstrip() for l in f.readlines() ]
    f.close()
    return lines

def Main():
        if args.p != None and args.f != None and args.t != None:
	      print "ERROR no taxon identifiers provided. See options flags by running with -h"
	      exit()

        if args.p != None and ( args.f != None or args.t != None ):
	      print "\n\nERROR cannot identify taxa specified both individually (-f or -t) AND as prefixes (-p)\n\n"
	      exit()

        # get the total set of sequence names to consider
        sequenceNames = []

        if args.f != None:
            sequenceNames += getNamesFromFile(args.f)
        if args.t != None:
            sequenceNames += args.t

        # read in the trees
        trees = list(Phylo.parse(args.treePath, 'newick')) # NB may only be one tree
        
        # do the replacements and write strings to file
        
        f = open(args.newTreePath, 'a+')
        for tree in trees:
	      if args.p != None: # the taxa possessing the prefixes may vary between trees, so can change the set of names for each tree
                  sequenceNames = []
                  sequenceNames = getNamesByPrefix(tree, args.p)
              else:
                  pass # using the sequenceNames set as defined above, with args.f and/or args.t options
              
              labelledTreeString = labelClade(tree, sequenceNames, args.l)
              f.write(labelledTreeString+os.linesep)
        f.close()

if __name__ == '__main__':
        Main()
