#!/usr/bin/env python

from Bio import SeqIO
from Bio.Blast import NCBIXML
import subprocess
import multiprocessing
import argparse
import itertools
#import networkx
import sys
import os
import md5

# functions

def getMd5(fileTuple):
    """
    
    """

    text = fileTuple[0] + '_' + fileTuple[1]
    md5sum = md5.new()
    md5sum.update(text)

    return md5sum.hexdigest()

def getBestHits(blastout, cS=0.7, cL=0.7, I=0.7):

    besthitsList = []

    for blastRecord in blastout:

        for alignment in blastRecord.alignments:

            identity = float(alignment.hsps[0].identities) / alignment.hsps[0].align_length

            if blastRecord.query_length > alignment.length:

                smallerSequenceCoverage = float(alignment.length)/alignment.hsps[0].align_length
                largerSequenceCoverage = float(blastRecord.query_length)/alignment.hsps[0].align_length

            else:

                smallerSequenceCoverage = float(blastRecord.query_length)/alignment.hsps[0].align_length
                largerSequenceCoverage = float(alignment.length)/alignment.hsps[0].align_length

            if identity >= I and smallerSequenceCoverage >= cS and largerSequenceCoverage >= cL:

                besthitsList.append((blastRecord.query, alignment.hit_def, identity))

                break

        else:

            besthitsList.append((blastRecord.query, None, None))
        
    return besthitsList

def getANI(bestHitsA, bestHitsB, cS=0.7, cL=0.7, I=0.7):

    selectedBestHits = []

    for bestHitA in bestHitsA:

        for bestHitB in bestHitsB:

            if bestHitA[0] == bestHitB[1] and bestHitA[1] == bestHitB[0]:

                selectedBestHits.append((bestHitA[0], bestHitB[0], float(bestHitA[2] + bestHitB[2])/2))
    
    AverageNucleotideIdentity = sum([record[2] for record in selectedBestHits])/len(selectedBestHits)

    return AverageNucleotideIdentity

# parse arguments

argumentParser = argparse.ArgumentParser(description="getANI: a simple tool to calculate average nucleotide identity (ANI)")
argumentParser.add_argument('-i', '--input_dir', required=True, help="a directory with two or more files in .gbk format")
argumentParser.add_argument('-o', '--output_dir', default=os.getcwd())
argumentParser.add_argument('-t', '--threads', default=multiprocessing.cpu_count(), type=int)
argumentParser.add_argument('-f', '--force', help="use to overwrite the output directory if it exists", action='store_true')
argumentParser.add_argument('-cS', '--coverage_smaller', default=0.7, type=float)
argumentParser.add_argument('-cL', '--coverage_larger',default=0.7, type=float)
argumentParser.add_argument('-I', '--identity_treshold', default=0.7, type=float)
arguments = argumentParser.parse_args()

# check if output directory exists

sys.stderr.write("Checking directories ... ")

if os.path.isdir(arguments.output_dir) and not arguments.force:

    sys.stderr.write("Error: output directory already exists! Use '-f/--force' if you want to overwrite it.\n")
    exit()

sys.stderr.write('ok!\n')

# parse inputs

sys.stderr.write("Parsing inputs ... ")

if os.path.isdir(arguments.input_dir):

    splitext = lambda inputFile:  os.path.splitext(inputFile)[1] == '.gbk'
    inputFiles = [inputFile for inputFile in os.listdir(arguments.input_dir) if splitext(inputFile)]
    inputFiles = [os.path.abspath("%s/%s"%(arguments.input_dir, inputFile)) for inputFile in inputFiles]

else:

    sys.stderr.write("Error: input directory not found!\n")
    exit()

sys.stderr.write('ok!\n')

# separate CDSs

sys.stderr.write("Extracting CDSs ... ")

cdsDir = "%s/cds_data"%(arguments.output_dir)

if not os.path.isdir(cdsDir):

    os.mkdir(cdsDir)

tempBlastDir = os.path.abspath("%s/tempblast_data"%(arguments.output_dir))

if not os.path.isdir(tempBlastDir):

    os.mkdir(tempBlastDir)

inputFilesCDSs = []
for inputFile in inputFiles:

    # parse input data

    inputFileHandler = open(inputFile)
    inputFileParser = SeqIO.parse(inputFileHandler, 'genbank')
    inputFileCDSsText = ""

    # iterate features

    for inputFileRecord in inputFileParser:

        locus_tag = ""

        for recordFeature in inputFileRecord.features:

            if 'locus_tag' in recordFeature.qualifiers.keys():

                locus_tag = recordFeature.qualifiers['locus_tag'][0]
            
            if recordFeature.type == 'CDS':

                inputFileCDSsText += '>{0}:{1}:{2}\n{3}\n'.format(
                    os.path.basename(inputFile),
                    inputFileRecord.id,
                    locus_tag,
                    str(recordFeature.extract(inputFileRecord.seq))
                )
    

    inputFileCdsPath = "%s/%s"%(cdsDir, os.path.basename(inputFile))
    inputFileCdsHandler = open(inputFileCdsPath, 'w')
    inputFileCdsHandler.write(inputFileCDSsText)
    inputFileCdsHandler.close()

    inputFilesCDSs.append(inputFileCdsPath)

sys.stderr.write('ok!\n')

# run blast

sys.stderr.write('Running BLASTn ... ')

bestHits = {}

for gA, genomeA in enumerate(inputFilesCDSs):

    makeblastdbCommandline = ("makeblastdb -in {0} "
                              "-dbtype nucl "
                              "-out {1}/{2} ").format(genomeA, 
                                                  tempBlastDir, 
                                                  os.path.basename(genomeA))

    subprocess.call(makeblastdbCommandline, shell=True,
                    stderr=open(os.devnull,'w'), stdout=open(os.devnull,'w')
                    )

    for gB, genomeB in enumerate(inputFilesCDSs):

        if genomeA == genomeB:

            continue

        blastnCommandline = ("blastn -query {0} "
                              "-db {1}/{2} "
                              "-out {1}/{3}_{4}.xml "
                              "-outfmt 5 "
                              "-num_threads {5}").format(genomeB, 
                                                         tempBlastDir,
                                                         os.path.basename(genomeA),
                                                         os.path.basename(genomeB),
                                                         os.path.basename(genomeA),
                                                         arguments.threads)
    
        subprocess.call(blastnCommandline,
                    shell=True#, 
                    #stderr=open(os.devnull,'w'), stdout=open(os.devnull,'w')
                    )
        blastoutAPath = "%s/%s_%s.xml"%(tempBlastDir, os.path.basename(genomeB), os.path.basename(genomeA))

        blastoutAParser = NCBIXML.parse(open(blastoutAPath))
        bestHits[(os.path.basename(genomeA), os.path.basename(genomeB))] = getBestHits(blastoutAParser)

        sys.stderr.write('\r\r\r%s %s                         '%(gA, gB))
    
    os.system("rm -r %s/*"%(tempBlastDir))

sys.stderr.write('ok!\n')

filePrefixes = [os.path.basename(file) for file in inputFilesCDSs]

fileCombinations = [combination for combination in itertools.combinations(filePrefixes, 2)]

ANIresults = []

outputTable = ""

for c, fileCombination in enumerate(fileCombinations):

    genomeA = fileCombination[0]
    genomeB = fileCombination[1]

    ANI = getANI(bestHits[(genomeA, genomeB)], bestHits[(genomeB, genomeA)])

    ANIresults.append((genomeA, genomeB, ANI))
    outputTable += "%s\t%s\t%s\n"%(genomeA, genomeB, (ANI * 100))

outputTableHandler = open("%s/result.tab"%(arguments.output_dir), 'w')
outputTableHandler.write(outputTable)
outputTableHandler.close()
