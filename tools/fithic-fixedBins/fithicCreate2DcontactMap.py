#!/bin/python
######################################
# Generate contact maps contactCount.gz
#
# Author: Fabian Buske (22/05/2015)
######################################

import os, sys, re
import traceback
from optparse import OptionParser
import fileinput
import datetime
from readData import *
from quicksect import IntervalTree
import gzip
from scipy.sparse import lil_matrix
import numpy

# manage option and arguments processing
def main():
    global options
    global args
    usage = '''usage: %prog [options] <contactCounts.gz>

    reads a fithic contactCounts.gz file and produces a full 2D contact matrix per chromosome
    '''
    parser = OptionParser(usage)
    parser.add_option("-q", "--quiet", action="store_false", dest="verbose", default=True,
                    help="don't print status messages to stdout")
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False,
                    help="print status messages to stdout")
    parser.add_option("-V", "--veryverbose", action="store_true", dest="vverbose", default=False,
                    help="print lots of status messages to stdout")
    parser.add_option("-O", "--onlycis", action="store_true", dest="onlycis", default=False,
                    help="only consider intra chromosomal contacts (cis)")
    parser.add_option("-r", "--resolution", type=int, dest="resolution", default=1000000,
                    help="size of a fragment in bp [default 1000000]")
    parser.add_option("-c", "--chromsizes", type="string", dest="chromSizes", default="",
                    help="tab separated file containing chromosome sizes")
    parser.add_option("-C", "--chrompattern", type="string", dest="chromPattern", default="",
                    help="pattern of chromosomes to filter for [default all]")
    parser.add_option("-o", "--outputDir", type="string", dest="outputDir", default="",
                    help="output directory [default: %default]")
    parser.add_option("-n", "--outputFilename", type="string", dest="outputFilename", default="",
                    help="output filename [default: extracted from first input file")
    parser.add_option("-t", "--tmpDir", type="string", dest="tmpDir", default="/tmp",
                    help="directory for temp files [default: %default]")
    parser.add_option("-s", "--sep", type="string", dest="separator", default=" ",
                    help="delimiter to use when reading the input [default: %default]")
    parser.add_option("--matrixFormat", type="string", dest="matrixFormat", default="tadbit",
                    help="either tadbit or domainfinder")
    parser.add_option("--inputIsFragmentPairs", action="store_true", dest="inputIsFragmentPairs", default=False,
                    help="input is a gzipped fragment pair file")
    parser.add_option("--inputIsReadPairs", type="string", dest="inputIsReadPairs", default="",
                    help="gzipped files with mapped read pair information, requires 4 column identifier corresponding to chrA,posA,chrB,posB,chrPrefix (separated buy comma), e.g. 2,3,6,7,chr")



    (options, args) = parser.parse_args()
    if (len(args) < 1):
        parser.print_help()
        parser.error("[ERROR] Incorrect number of arguments, need a dataset")

    if (options.resolution < 1):
        parser.error("[ERROR] resolution must be a positive integer, was :"+str(options.resolution))
        sys.exit(1)
    elif (options.chromSizes == "" or not os.path.isfile(options.chromSizes)):
        parser.error("[ERROR] chromSizes not given or not existing, was :"+str(options.chromSizes))
        sys.exit(1)

    if (options.outputDir != ""):
        options.outputDir += os.sep

    if (options.verbose):
        print >> sys.stdout, "resolution:            %s" % (options.resolution)
        print >> sys.stdout, "chromSizes:            %s" % (options.chromSizes)
        print >> sys.stdout, "outputDir:             %s" % (options.outputDir)
        print >> sys.stdout, "tmpDir:                %s" % (options.tmpDir)

    process()


def output(fragmentsMap , fragmentList, fragmentPairs, fragmentCount, fragmentsChrom):
    '''
    outputs the 2D contact matrix
    '''

    if (options.verbose):
        print >> sys.stdout, "- %s START   : output data " % (timeStamp())

    fragmentIds = fragmentsMap.keys()
    fragmentIds.sort()

    # populate sparse matrix
    A = lil_matrix((fragmentCount, fragmentCount), dtype='i')
    for fragmentIds, contactCounts in fragmentPairs.iteritems():
        A[fragmentIds[0],fragmentIds[1]] = contactCounts
        A[fragmentIds[1],fragmentIds[0]] = contactCounts
    # convert to coordinate format
    B = A.tocoo()

    for chr in fragmentsChrom.keys():

        C = B.tocsc()[:,fragmentsChrom[chr][0]:fragmentsChrom[chr][1]].tocsr()[fragmentsChrom[chr][0]:fragmentsChrom[chr][1],:]

        fragmentRange=fragmentsChrom[chr][1]-fragmentsChrom[chr][0]
        header=['d']+[ "%s%d" % i for i in zip(['r']*fragmentRange,range(fragmentRange))]

        if ( options.outputFilename != "" ):
            outfile3 = options.outputDir+options.outputFilename+"."+chr+".matrix"
        else:
            outfile3 = options.outputDir+os.path.basename(args[0])+"."+chr+".matrix"

        if (options.verbose):
            print >> sys.stdout, "- save 2Dmatrix for chromosome %s to %s " % (chr, outfile3)

        f_handle=open(outfile3,'w')
        if (options.matrixFormat == "domainfinder"):
            for i in xrange(fragmentRange):
                binStart = fragmentsMap[fragmentsChrom[chr][0]+i][1] - options.resolution/2
                binEnd = binStart + options.resolution

                f_handle.write(chr+"\t"+str(binStart)+"\t"+str(binEnd)+"\t")
                numpy.savetxt(f_handle, C[i].toarray(),fmt='%i', delimiter='\t')
        else:
            f_handle.write('\t'.join(header)+"\n")

            for i in xrange(fragmentRange):
                f_handle.write(header[i+1]+"\t")
                numpy.savetxt(f_handle, C[i].toarray(),fmt='%i', delimiter='\t')

        f_handle.close()

    if (options.verbose):
        print >> sys.stdout, "- %s FINISHED: output data" % (timeStamp())


def process():
    global options
    global args

    [ fragmentsMap, lookup_structure, fragmentCount, fragmentsChrom ] = createIntervalTreesFragmentResolution(options)

    [ fragmentList, fragmentPairs ] = countReadsPerFragment(lookup_structure, options, args)

    output(fragmentsMap, fragmentList, fragmentPairs, fragmentCount, fragmentsChrom)


######################################
# main
######################################
if __name__ == "__main__":
    main()
