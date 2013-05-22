import os, sys
import logging
from optparse import OptionParser
import matplotlib
matplotlib.use('PDF')

from hiclib import mapping
from mirnylib import h5dict, genome
from hiclib import fragmentHiC

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from mirnylib import plotting
from hiclib import binnedData


# manage option and arguments processing
def main():
	global options
	global args
	usage = '''usage: %prog [options] reads.[fastq|sra|bam]+

takes fastq or sra files and runs the hiclib pipeline on it
Note, read pairs in fastq format (possible gzipped) or bam need to be stated next to each other, i.e. fastq_r1 fastq_r2
	'''
	parser = OptionParser(usage)
	parser.add_option("-q", "--quiet", action="store_false", dest="verbose", default=True,
					help="don't print status messages to stdout")
	parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
					help="print status messages to stdout")
	parser.add_option("-e", "--restrictionEnzyme", type="string", dest="enzyme", default="", 
					help="Name of the restriction enzyme, e.g. BglII")
	parser.add_option("-n", "--experimentName", type="string", dest="experiment", default="", 
					help="Name of the experiment")
	parser.add_option("-b", "--bowtie", type="string", dest="bowtie", default="", 
					help="location of bowtie [default: %default]")
	parser.add_option("-r", "--referenceGenome", type="string", dest="genome", default="", 
					help="genome in fasta format [default: %default]")
	parser.add_option("-g", "--gapFile", type="string", dest="gapFile", default="",
					help="location of the gapfile [default: %default]")
	parser.add_option("-i", "--index", type="string", dest="index", default="", 
					help="location of genome index including the basename")
	parser.add_option("-l", "--readLength", type="int", dest="readLength", default=100, 
					help="length of the reads [default: %default]")
	parser.add_option("-f", "--inputFormat", type="string", dest="inputFormat", default="fastq", 
					help="format of the input file, either fastq, sra or bam [default: %default]")
	parser.add_option("-o", "--outputDir", type="string", dest="outputDir", default="", 
					help="output directory [default: %default]")
	parser.add_option("-c", "--cpus", type="int", dest="cpus", default=1, 
					help="number of cpus to use [default: %default]")
	parser.add_option("-t", "--tmpDir", type="string", dest="tmpDir", default="/tmp", 
					help="directory for temp files [default: %default]")
	parser.add_option("-s", "--sra-reader", type="string", dest="sra", default="fastq-dump", 
					help="location of sra reader fastq-dump in case input is SRA [default: %default]")
	
	(options, args) = parser.parse_args()
	if (len(args) < 1):
		parser.print_help()
		parser.error("[ERROR] Incorrect number of arguments, need at least one read file")

	if (options.inputFormat != 'fastq' and options.inputFormat != 'sra' and options.inputFormat != 'bam'):
		print >> sys.stderr, "[ERROR] Input format not supported: %s" % (options.inputFormat)
		sys.exit(1)	

	if ((options.inputFormat == 'fastq' or options.inputFormat == 'bam') and len(args) % 2 != 0):
		print >> sys.stderr, "[ERROR] Both reads are required for files in fastq"
		sys.exit(1)	

	if (options.genome == ""):
		print >> sys.stderr, "[ERROR] Please specify the location of the reference genome in fasta format"
		sys.exit(1)

	if (options.inputFormat != 'bam' and options.index == ""):
		print >> sys.stderr, "[ERROR] Please specify the location of the bowtie2 index for the reference genome"
		sys.exit(1)
		
	if (options.enzyme == ""):
		print >> sys.stderr, "[ERROR] Please specify the restriction enzyme (supported enzymes: http://www.biopython.org/DIST/docs/api/Bio.Restriction-module.html)"
		sys.exit(1)

	if (options.experiment == ""):
		print >> sys.stderr, "[ERROR] Please provide a name for the experiment, e.g. [Cellline]_[Enzymename]_[Replica]"
		sys.exit(1)
	
	if (options.outputDir != ""): 
		options.outputDir += os.sep
	

	if (options.verbose):
		print >> sys.stdout, "restrictionEnzyme:  %s" % (options.enzyme)
		print >> sys.stdout, "experimentName:     %s" % (options.experiment)
		print >> sys.stdout, "bowtie:             %s" % (options.bowtie)
		print >> sys.stdout, "referenceGenome:    %s" % (options.genome)
		print >> sys.stdout, "index:              %s" % (options.index)
		print >> sys.stdout, "readLength:         %d" % (options.readLength)
		print >> sys.stdout, "outputDir:          %s" % (options.outputDir)
		print >> sys.stdout, "tmpDir:             %s" % (options.tmpDir)
		print >> sys.stdout, "cpus:               %s" % (options.cpus)
		print >> sys.stdout, "inputFormat:        %s" % (options.inputFormat)
		print >> sys.stdout, "sra-reader:         %s" % (options.sra)

	process()

def correctedScalingPlot(resolution, filename, experiment, genome, mouse=False, **kwargs):
    "Paper figure to compare scaling before/after correction"
    
    global pp
    if (options.verbose):
        print >> sys.stdout, "correctedScalingPlot: res: %d file1: %s exp1:%s gen:%s" % (resolution, filename, experiment, genome)

    #plt.figure(figsize=(4, 4))
    Tanay = binnedDataAnalysis(resolution, genome)
    Tanay.simpleLoad(filename, experiment)
    Tanay.removePoorRegions()
    Tanay.removeDiagonal()
    Tanay.plotScaling(experiment, label="Raw data", color="#A7A241")
    Tanay.iterativeCorrectWithSS()
    Tanay.plotScaling(experiment, label="Corrected", color="#344370")
    ax = plt.gca()
    mirnylib.plotting.removeAxes()
    fs = 6
    plt.xlabel("Genomic distance (MB)", fontsize=6)
    plt.ylabel("Contact probability", fontsize=6)
    for xlabel_i in ax.get_xticklabels():
        xlabel_i.set_fontsize(fs)
    for xlabel_i in ax.get_yticklabels():
        xlabel_i.set_fontsize(fs)
    legend = plt.legend(loc=0, prop={"size": 6})
    legend.draw_frame(False)
    plt.xscale("log")
    plt.yscale("log")
    plt.show()
    pp.savefig()

def doArmPlot(resolution, filename, experiment, genome, mouse=False, **kwargs):
    "Plot an single interarm map - paper figure"

    global pp    
    if (options.verbose):
        print >> sys.stdout, "doArmPlot: res: %d file: %s exp:%s gen:%s" % (resolution, filename, experiment, genome)

    Tanay = binnedDataAnalysis(resolution, genome)
    Tanay.simpleLoad(filename, experiment)
    if mouse == True:
        Tanay.fakeTranslocations([(0, 0, None, 12, 52000000, None),
                                  (4, 45000000, None, 12, 0, 30000000),
                                  (9, 0, 50000000, 12, 0, 35000000)])
        Tanay.removeChromosome(19)
    else:
        Tanay.removeChromosome(22)
    Tanay.removeDiagonal(1)
    Tanay.removePoorRegions()
    Tanay.truncTrans()
    Tanay.fakeCis()
    #mat_img(Tanay.dataDict["GM-all"])
    #plt.figure(figsize = (3.6,3.6))
    Tanay.averageTransMap(experiment, **kwargs)

    #plotting.removeBorder()
    cb = plt.colorbar(orientation="vertical")
    #cb.set_ticks([-0.05,0.05,0.15])
    for xlabel_i in cb.ax.get_xticklabels():
        xlabel_i.set_fontsize(6)

def mapFile(fastq, read):
	global options
	global args

	fileName, fileExtension = os.path.splitext(fastq)
	bamOutput = options.outputDir+fileName.split(os.sep)[-1]+'_R'+str(read)+'.bam'
	
	if (fileExtension == '.sra'):
		if (options.verbose):
			print >> sys.stdout, "Map short read archive %s utilizing %s" % (fastq, options.sra)

		mapping.iterative_mapping(
		    bowtie_path=options.bowtie,
		    bowtie_index_path=options.index,
		    fastq_path=fastq,
		    out_sam_path=bamOutput,
		    min_seq_len=25,
		    len_step=5,
		    seq_start=options.readLength*(read-1),
		    seq_end=options.readLength*(read),
		    nthreads=options.cpus,
		    temp_dir=options.tmpDir, 
		    bowtie_flags='--very-sensitive',
		    bash_reader=options.sra+' -Z')
	
	else:
		if (options.verbose):
			print >> sys.stdout, "Map fastq %s" % (fastq)
		
		mapping.iterative_mapping(
		    bowtie_path=options.bowtie,
		    bowtie_index_path=options.index,
		    fastq_path=fastq,
		    out_sam_path=bamOutput,
		    min_seq_len=25,
		    len_step=5,
		    nthreads=options.cpus,
		    temp_dir=options.tmpDir, 
		    bowtie_flags='--very-sensitive')
		    
	return bamOutput


def mapFiles():

	bams = []
	if (options.inputFormat == 'fastq'):
	
		if (options.verbose):
			print >> sys.stdout, "**  Process fastq files"

		for i in range(0, len(args),2):
			
			if (options.verbose):
				print >> sys.stdout, "**  Map first input file"
			bams+=[mapFile(args[i], 1)]

			if (options.verbose):
				print >> sys.stdout, "**  Map second input file"
		
			bams+=[mapFile(args[i+1], 2)]
	else:
		if (options.verbose):
			print >> sys.stdout, "**  Process sra files"

		for i in range(0, len(args)):
			
			if (options.verbose):
				print >> sys.stdout, "**  Map first input file"
			bams+=[mapFile(args[i], 1)]

			if (options.verbose):
				print >> sys.stdout, "**  Map second input file"
		
			bams+=[mapFile(args[i], 2)]
	
	return bams


def collectMappedReads(bam_read1, bam_read2, mapped_reads, genome_db):
	global options
	global args
	
	mapping.parse_sam(
	    sam_basename1=bam_read1,
	    sam_basename2=bam_read2,
	    out_dict=mapped_reads,
	    genome_db=genome_db, 
	    enzyme_name=options.enzyme)

def filterFragments(genome_db):
	'''
	Filter the data at the level of individual restriction fragments

	The following reads are remove from the dataset:

	- the reads that start within the 5 bp range from the restriction site
	- the identical read pairs, with both ends starting at exactly the same positions
	- the reads coming from extremely large and extremely small restriction fragments (length > 10^5 bp or length < 100 bp)
	- the reads coming from the top 0.5% most frequently detected restriction fragments

	The rationale behind each of the filters is discussed in the hiclib publication. The API documentation contains the description of the filters.
	'''
	
	fragments = fragmentHiC.HiCdataset(
	    filename=options.outputDir+options.experiment+'-fragment_dataset.hdf5',
	    genome=genome_db,
	    maximumMoleculeLength=500,
	    mode='w')
	
	# Load the parsed reads into the HiCdataset. The dangling-end filter is applied
	# at this stage, with maximumMoleculeLength specified at the initiation of the 
	# object.
	fragments.parseInputData(dictLike=options.outputDir+options.experiment+'-mapped_reads.hdf5')
	
	fragments.filterRsiteStart(offset=5)
	fragments.filterDuplicates()
	
	fragments.filterLarge()
	fragments.filterExtreme(cutH=0.005, cutL=0)
	
	fragments.saveHeatmap(options.outputDir+options.experiment+'-1M.hdf5', resolution=1000000)

	fragments.saveHeatmap(options.outputDir+options.experiment+'-200k.hdf5', resolution=200000)
	
	return fragments

def iterativeFiltering(genome_db, fragments, filesuffix):
	'''
	Filter the data at the binned level and perform the iterative correction.
	'''
	
	# Read resolution from the dataset.
	raw_heatmap = h5dict.h5dict(options.outputDir+options.experiment+filesuffix, mode='r') 
	resolution = int(raw_heatmap['resolution'])
	
	# Create a binnedData object, load the data.
	BD = binnedData.binnedData(resolution, genome_db)
	BD.simpleLoad(options.outputDir+options.experiment+filesuffix, options.experiment)

	# Remove the contacts between loci located within the same bin.
	BD.removeDiagonal()
	
	# Remove bins with less than half of a bin sequenced.
	BD.removeBySequencedCount(0.5)
	
	# Remove 1% of regions with low coverage.
	BD.removePoorRegions(cutoff=1)
	
	# Truncate top 0.05% of interchromosomal counts (possibly, PCR blowouts).
	BD.truncTrans(high=0.0005)
	
	# Perform iterative correction.
	BD.iterativeCorrectWithoutSS()

	# Save the iteratively corrected heatmap.
	BD.export(options.experiment, options.outputDir+options.experiment+'-IC'+filesuffix)

	plotting.plot_matrix(np.log(BD.dataDict[options.experiment]))
	pp.savefig()

def process():
	global options
	global args
	global pp
	
	if (options.verbose):
		print >> sys.stdout, "*** START processing"

	fig = plt.gcf()
	pp = PdfPages(options.outputDir+options.experiment+'.pdf')
	
	logging.basicConfig(level=logging.DEBUG)
	
	if (options.verbose):
		print >> sys.stdout, "**  Create directories"

	if not os.path.exists(options.tmpDir):
		os.mkdir(options.tmpDir)

	if not os.path.exists(options.outputDir):
		os.mkdir(options.outputDir)
	
	if (options.verbose):
		print >> sys.stdout, "**  Create data objects"

	mapped_reads = h5dict.h5dict(options.outputDir+options.experiment+'-mapped_reads.hdf5')
	genome_db    = genome.Genome(options.genome, gapFile=options.gapFile, readChrms=['#', 'X', 'Y'])

	bams = []
	if (options.inputFormat != 'bam'):
		bams = mapFiles()
	else:
		bams = args[0:]

	if (options.verbose):
		print >> sys.stdout, "**  Collect mapped reads"
		
	collectMappedReads(bams[0], bams[1], mapped_reads, genome_db)
	
	if (options.verbose):
		print >> sys.stdout, "**  Filter fragments"
	
	fragments = filterFragments(genome_db)
	
	if (options.verbose):
		print >> sys.stdout, "**  Iterative filtering of fragments"

	iterativeFiltering(genome_db, fragments, '-1M.hdf5')
	
	iterativeFiltering(genome_db, fragments, '-200k.hdf5')

	# plotting
	correctedScalingPlot(200000, options.outputDir+options.experiment+'-200k.hdf5', options.experiment, genome_db)

	doArmPlot(1000000, options.outputDir+options.experiment+'-1M.hdf5', options.experiment, genome_db)

	if (options.verbose):
		print >> sys.stdout, "*** FINISHED processing"
	
	pp.close()
	
	
######################################
# main
######################################
if __name__ == "__main__":
	main()
