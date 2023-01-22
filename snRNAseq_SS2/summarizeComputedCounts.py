# collect counts
import sys

dataDir = sys.argv[1]

import os,glob

os.chdir(dataDir)

counts = {}
geneNames = []
header = []

# this script is adjusted for otuput from feature counts



for f in glob.glob("*.counts"):
    
    countsFile = open(f)
    
    i = 0   
    for line in countsFile:
        
        if line.startswith("#"):
            continue

        if line.startswith("Geneid"):
            sName = f.split(".")[0]
            #if sName in sampleNames:
            header.append( sName )
            continue
            #else:
            #    print "Sample name %s is not found!" % sName
            #    break
                 
        items = line.split()     
        geneName = items[0]
        c = items[-1]
        if geneName in counts:
            counts[geneName].append(c)
        else:
            counts[geneName] = [ c ]
            geneNames.append(geneName)
        i += 1
    print "Processed %d genes from %s" % (i,f) 
    
outfile = open("counts_full.txt", "w")

outfile.write( "\t".join(header) + "\n") 

for g in geneNames:
    c = counts[g]
    outfile.write( g + "\t" + "\t".join(c) + "\n")


