#!/usr/bin/python
 
description = """ 
Author: Christopher Monit

While BioPython's PAML package provides a parsing tool, it does not (currently) provide parsing of the Bayes Emprical Bayes results

This code is for conventiently obtaining BEB sites, including in 'machine friendly' formats. At present it works with the main results file only, so can present data for sites whose BEB prob>0.5

[Improvements could (should) be made to access exact BEB prob values, and interpreset comprehensive output in rst file]
"""

import argparse
import re

def getSection(allLines, startPattern, endPattern): # this is modified (improved) version of something in get_codeml_results.py
    
    """ Given all of the lines which comprise a codeml results file, return those lines which are beween start and end patterns """

    section = []

    for i in range(len(allLines)):
        if allLines[i].__contains__(startPattern):
            
            for j in range(i, len(allLines)):
                
                if allLines[j].__contains__(endPattern): # we've found the end of the section we want, so let's return what we have now
                    return section
                else:
                    section.append(allLines[j]) # we have not reached the end pattern, so this is a line in between (which we should include)
            else:
                raise Exception('getSection: could not find endPattern')
    else:
        raise Exception('getSection: could not find startPattern')


def getSignificant(bebSection, sigPattern="*"): # need "* " as pattern because beb header includes "*"
    sigSites = []
    for line in bebSection:
        if bool(re.search(r'\d\*', line)): # if the string contains a digit, followed by an astrisk
            split = line.rsplit()
            site = split[0]
            sigSites.append(site)
            #beb = split[2].remove("*") # the prob value, stored as string
    return sigSites

def printSection(section):
    for line in section:
        if line != "":
            print line

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument( "codeml_output" )
    parser.add_argument( "-n", action="store_true", help="Print the number of sites which have BEB>0.95 (i.e. are denoted by * or ** in file)" )
    parser.add_argument( "-p", action="store_true", help="Print the BEB section of results file, rather than just the sites" )
    parser.add_argument( "-l", action="store_true", help="Print the sites which have BEB>0.95 on one single line, delimited by comma" )
    parser.add_argument( "-start", default="Bayes Empirical Bayes", help="Start pattern for section (default=\"Bayes Empirical Bayes\")" )
    parser.add_argument( "-end", default="The grid",help="End pattern for section (default=\"The grid\")" )
    args = parser.parse_args()
    
    f = open(args.codeml_output)
    lines = [ l.rstrip() for l in f.readlines() ]
    f.close()

    bebSection = getSection(lines, args.start, args.end)
    
    significant = getSignificant(bebSection)

    if args.n:
        print len( significant )

    if args.p:
        printSection(bebSection)
    
    delimiter = ","

    if args.l:
        print delimiter.join(significant)

    if not args.n and not args.p and not args.l: # default
        for site in significant:
            print site


