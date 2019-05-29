#!/usr/bin/python

description = """ 
Author: Christopher Monit

Reroot one or more newick trees.
Tree is rooted on the branch leading to the most recent common ancestor of a set of taxa, whose names are provided in a separate text file.

NB root is placed exactly half way along the given branch;
i.e. rooting on a branch of length t, the branches leading from root to outgroup clade and ingroup clade will each be of length t/2.0

"""


from Bio import Phylo
import Bioplus.TreeExtras as TreeExtras

# rooting a tree
def rootByOutgroupNames(tree, outgroupNames):
    """ Root a tree on the branch leading to the common ancestor of the outgroup taxa
    NB clade labels may be disturbed by this approach - check result carefully if your tree has clade labels/support values""" 
    # get outgroup clade
    outgroupTerminals = TreeExtras.findCladesLongNames(tree, outgroupNames)
    outgroupClade = tree.common_ancestor(outgroupTerminals)
    bl = outgroupClade.branch_length / 2.0

    tree.root_with_outgroup(outgroupClade, outgroup_branch_length=bl)
    return tree # this this necessary?


def readOutgroupNames(outgroupNamesFilePath):
    """ Read the names of the outgroup taxa from the given file path """
    outgroupFile = open(outgroupNamesFilePath)
    outgroupNames = [ line.rstrip() for line in outgroupFile.readlines() ]
    outgroupFile.close()
    return outgroupNames


# for rooting a single tree from a single file
def rootTreeNamesFiles(inputTreeFilePath, outputTreeFilePath, outgroupNamesFilePath):
    """ Given a file containing a single newick tree and a set of taxon names which comprise the outgroup, root the tree on the branch leading to the outgroup common ancestor and write to the given output file name """
    tree = Phylo.read(inputTreeFilePath, 'newick')
    
    outgroupNames = readOutgroupNames(outgroupNamesFilePath) 

    treeNew = rootByOutgroupNames(tree, outgroupNames)
    Phylo.write(treeNew, outputTreeFilePath, 'newick') 
    
# for rooting multiple trees - can accept a file containing several
def multipleTreeFiles(inputTreeFilePath, outputTreeFilePath, outgroupNamesFilePath):
    """ Given a file containing a set of newick trees and a set of taxon names which comprise the outgroup, root each tree on the branch leading to the outgroup common ancestor and write the trees to the given output file name """
    trees = list(Phylo.parse(inputTreeFilePath, 'newick'))
    outgroupNames = readOutgroupNames(outgroupNamesFilePath) 
    for tree in trees:
        tree = rootByOutgroupNames(tree, outgroupNames)
    
    Phylo.write(trees, outputTreeFilePath, 'newick') 


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument( "-m", action='store_true', help="Multiple trees. Input file contains more than one tree which will need rooting with the same outgroup" )
    parser.add_argument( "tree_file", help="file containing one or more Newick format trees" )
    parser.add_argument( "outgroup_taxa", help="file containing names of outgroup taxa. The tree will be rooted on the branch leading to the node which is their most recent common ancestor" )
    parser.add_argument( "new_tree_file", help="name for new file to which the resulting tree will be written (in Newick format)" )
    args = parser.parse_args()

    if not args.m: # not dealing with multiple trees
        rootTreeNamesFiles( args.tree_file, args.new_tree_file, args.outgroup_taxa )
    else:
        multipleTreeFiles( args.tree_file, args.new_tree_file, args.outgroup_taxa )


