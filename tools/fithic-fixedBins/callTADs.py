#!/bin/python
######################################
# Calls topological domains 
# With TADbit
#
# Author: Fabian Buske (4/12/2014)
######################################

# deactivate interactive plotting
import matplotlib
matplotlib.use('Agg')

import os, sys, re
import traceback
from optparse import OptionParser
import fileinput
import datetime
from pytadbit import Chromosome

######################################
# Timestamp
######################################
def timeStamp():
    return datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S').format()

def process():

    if ( options.outputFilename != "" ):
        outfilefileprefix=options.outputDir+options.outputFilename
    else:
        outfilefileprefix=options.outputDir+os.path.basename(args[0])
    
    for matrixFile in xrange(len(args)):
        sample=os.path.splitext(os.path.basename(args[matrixFile]))[0].split(".matrix")[0]
        chr = sample.rsplit(".",1)[-1]
        sample = sample.rsplit(".",1)[0]
        chrom = Chromosome(name=chr, centromere_search=True, species=options.species, assembly=options.assembly)
        chrom.set_max_tad_size(5000000)
        chrom.add_experiment(sample, exp_type='Hi-C', identifier=sample,
                        hic_data=args[matrixFile], resolution=options.resolution)
        
        exp = chrom.experiments[sample]
        exp.normalize_hic(silent=True)
        chrom.find_tad(sample, n_cpus=options.threads, normalized=True, verbose=False)
        exp.write_tad_borders(outfilefileprefix+"."+chr+".border")

        chrom.tad_density_plot(sample,savefig=outfilefileprefix+".density."+chr+".pdf")
        chrom.visualize(exp.name, paint_tads=True, savefig=outfilefileprefix+"chr."+chr+".pdf")
        chrom.save_chromosome(outfilefileprefix+"chr."+chr+".tdb", force=True)

# manage option and arguments processing
def main():
    global options
    global args
    usage = '''usage: %prog [options] [contactMatrix]+
    
    calls TADbit functions to infer topological domains from contact matrices
    '''
    parser = OptionParser(usage)
    parser.add_option("-q", "--quiet", action="store_false", dest="verbose", default=True,
                    help="don't print status messages to stdout")
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False,
                    help="print status messages to stdout")
    parser.add_option("-V", "--veryverbose", action="store_true", dest="vverbose", default=False,
                    help="print lots of status messages to stdout")
    parser.add_option("-r", "--resolution", type=int, dest="resolution", default=1000000, 
                    help="size of a fragment in bp if no genomeFragmentFile is given [default %default]")
    parser.add_option("-s", "--species", type="string", dest="species", default="", 
                    help="species, e.g. homosapiens, [default %default]")
    parser.add_option("-a", "--assembly", type="string", dest="assembly", default="", 
                    help="genome assembly, e.g. hg19 [default %default]")
    parser.add_option("-o", "--outputDir", type="string", dest="outputDir", default="", 
                    help="output directory [default: %default]")
    parser.add_option("-t", "--threads", type=int, dest="threads", default=1, 
                    help="number of CPUs to use [default: %default]")
    parser.add_option("-n", "--outputFilename", type="string", dest="outputFilename", default="", 
                    help="prefix for the output files [default: extracted from first input file")

    (options, args) = parser.parse_args()
    if (len(args) < 1):
        parser.print_help()
        parser.error("[ERROR] Incorrect number of arguments, need a dataset")

    else:
        if (options.resolution < 1):
            parser.error("[ERROR] resolution must be a positive integer, was :"+str(options.resolution))
            sys.exit(1)
        
    if (options.outputDir != ""): 
        options.outputDir += os.sep
    
    if (options.verbose):
        print >> sys.stdout, "resolution:            %s" % (options.resolution)
        print >> sys.stdout, "outputDir:             %s" % (options.outputDir)
        print >> sys.stdout, "tmpDir:                %s" % (options.tmpDir)
    
    process()


######################################
# main
######################################
if __name__ == "__main__":
    main()

