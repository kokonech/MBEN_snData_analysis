import sys
import argparse
import HTSeq
import numpy as np


def idsContainGiven(givenId, transcriptIds):
    for tId in transcriptIds:
        if givenId.find(tId) != -1:
            return True

    return False


ANN_FILE="/omics/odcf/analysis/OE0290_projects/Ependymoma/annotations/star_index_gencode19/gencode.v19.annotation.clean.gtf"

if __name__ == "__main__":
    
    descriptionText = "The script performs convertion of gene IDs to names in counts file."

    parser = argparse.ArgumentParser(description = descriptionText,formatter_class=argparse.RawDescriptionHelpFormatter)
    
         
    parser.add_argument("-i", action="store", required="true", dest="inFile",
        help="Input file.")
        
    parser.add_argument("-g", action="store",  dest="gtfFile",
        default=ANN_FILE, help="Input file with list of genes in GTF format. Default: H.Sapiens gencode v19 ")

    parser.add_argument("-o", action="store", required="true", dest="outFile",
        help="Output file.")

    parser.add_argument("--keep-coords", action="store_true", default=False, dest="keepCoords",
        help="Output file.")

   
    args = parser.parse_args()
    
    print args
    
    inFileName = args.inFile
    gtfFileName = args.gtfFile
    outFileName = args.outFile

    # parse GTF file

    gtf_file = HTSeq.GFF_Reader( gtfFileName )

    geneIdMap = {}

    geneStarts = {}
    geneEnds = {}
    geneChrs = {}

    for feature in gtf_file:
        if feature.type == 'exon':
            geneId = feature.attr[ "gene_id" ]
            geneName = feature.attr[ "gene_name" ]
            geneIdMap[geneId] = geneName

            if args.keepCoords:
                
                geneChrs[geneId] = feature.iv.chrom

                if geneId in geneStarts:
                    startPos = feature.iv.start
                    prevStartPos = geneStarts[geneId]
                    if startPos < prevStartPos:
                        geneStarts[geneId] = startPos    
                else:
                    geneStarts[geneId] = feature.iv.start      

                if geneId in geneEnds:
                    endPos = feature.iv.end
                    prevEndPos = geneEnds[geneId]
                    if endPos > prevEndPos:
                        geneEnds[geneId] = endPos    
                else:
                    geneEnds[geneId] = feature.iv.end



    
    outFile = open(outFileName, "w")
    
    print "Writing result ...."

    renamed = 0
    total = 0

    genenames = set()

    for line in open(inFileName):
    
        geneId = line.split()[0]
        total += 1
        if total == 1:
            # write header
            outFile.write(line)
            continue   

        if geneId in geneIdMap:
            geneName = geneIdMap[geneId]
            if geneName in genenames:
                print "WARNING! Gene name is already present:", geneName
                continue

            if args.keepCoords:
                adjustedGeneName = "%s\t%s\t%s\t%s" % (geneName,geneChrs[geneId],geneStarts[geneId], geneEnds[geneId])
                line = line.replace(geneId, adjustedGeneName)
            else:
                line = line.replace(geneId, geneName)
    
            renamed += 1
            genenames.add(geneName)

            outFile.write(line)
        

    print "Extracted %d genes out of %d" % (renamed,total)
