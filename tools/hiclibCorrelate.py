import os
import sys
from optparse import OptionParser

import matplotlib
matplotlib.use('PDF')
from matplotlib.backends.backend_pdf import PdfPages

from mirnylib.systemutils import setExceptionHook
sys.path.append(os.path.split(os.getcwd())[0])
from hiclib.binnedData import binnedData, binnedDataAnalysis,\
    experimentalBinnedData
from mirnylib import h5dict, genome

import mirnylib.systemutils
mirnylib.systemutils
from hiclib.fragmentHiC import HiCdataset
from mirnylib.numutils import EIG, coarsegrain, project, arrayInArray
import numpy
import mirnylib.plotting
import scipy.stats
import scipy.ndimage
from hiclib import fragmentHiC
cr = scipy.stats.spearmanr
import cPickle
from mirnylib.plotting import mat_img, removeAxes, removeBorder, niceShow

import matplotlib.pyplot as plt

# manage option and arguments processing
def main():
	global options
	global args
	usage = '''usage: %prog [options] dataset1 dataset2 [datasetX]*

takes multiple hiclib output folder and compares the experiments in a pairwise manner
	'''
	parser = OptionParser(usage)
	parser.add_option("-q", "--quiet", action="store_false", dest="verbose", default=True,
					help="don't print status messages to stdout")
	parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
					help="print status messages to stdout")
	parser.add_option("-e", "--restrictionEnzyme", type="string", dest="enzyme", default="", 
					help="Name (or comma separated list) of the restriction enzyme, e.g. BglII")
	parser.add_option("-n", "--experimentNames", type="string", dest="experiment", default="", 
					help="Comma separated list of the experiment names")
	parser.add_option("-r", "--referenceGenome", type="string", dest="genome", default="", 
					help="genome in fasta format [default: %default]")
	parser.add_option("-g", "--gapFile", type="string", dest="gapFile", default="gap.txt",
					help="location of the gapfile [default: %default]")
	parser.add_option("-o", "--outputDir", type="string", dest="outputDir", default="", 
					help="output directory [default: %default]")
	parser.add_option("-t", "--tmpDir", type="string", dest="tmpDir", default="/tmp", 
					help="directory for temp files [default: %default]")

	(options, args) = parser.parse_args()
	if (len(args) < 2):
		parser.print_help()
		parser.error("[ERROR] Incorrect number of arguments, need at least two datasets")
		
	if (options.enzyme == ""):
		print >> sys.stderr, "[ERROR] Please specify the restriction enzyme(s). Supported enzymes: http://www.biopython.org/DIST/docs/api/Bio.Restriction-module.html)"
		sys.exit(1)

	if (len(options.enzyme.split(',')) !=  1 and len(options.enzyme.split(',')) != len(args)):
		print >> sys.stderr, "[ERROR] Number of restriction enzymes must be either 1 (same enzymes for all data) or match the number of datasets"
		sys.exit(1)

	if (len(options.experiment.split(',')) != len(args)):
		print >> sys.stderr, "[ERROR] Please provide the (base-)name for each dataset, e.g. [Cellline]_[Enzymename]_[Replica]"
		sys.exit(1)
	
	if (options.genome == ""):
		print >> sys.stderr, "[ERROR] Please specify the location of the reference genome in fasta format"
		sys.exit(1)

	if (options.outputDir != ""): 
		options.outputDir += os.sep
	

	if (options.verbose):
		print >> sys.stdout, "restrictionEnzyme:  %s" % (options.enzyme)
		print >> sys.stdout, "experimentName:     %s" % (options.experiment)
		print >> sys.stdout, "outputDir:          %s" % (options.outputDir)
		print >> sys.stdout, "tmpDir:             %s" % (options.tmpDir)

	process()


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

def calculateTanayCorrelation(resolution, filename1, filename2, experiment1, experiment2, genome, mouse=False, **kwargs):
    "Calculates correlation between datasets, smoothed in a Tanay way"

	global pp
    if (options.verbose):
        print >> sys.stdout, "calculateTanayCorrelation: res: %d file1: %s file2: %s exp1:%s exp2:%s gen:%s" % (resolution, filename1, filename2, experiment1, experiment2, genome)

    BD = binnedData(resolution, genome)
    BD.simpleLoad(filename1, experiment1)
    BD.simpleLoad(filename2, experiment2)

    def tanaySmooth(matrix):
        matrix = numpy.array(matrix, dtype="double")
        a = numpy.arange(-9, 10)
        mat = 1 / (1. + numpy.abs(a[:, None]) + numpy.abs(a[None, :]))
        return scipy.ndimage.filters.convolve(input=matrix,
                                              weights=mat,
                                              mode="constant")

    def propagateSmooth(data):
        mask1 = numpy.sum(data, axis=0) > 0
        mask = mask1[:, None] * mask1[None, :]
        ret = numpy.zeros_like(data, dtype=float)
        for i in xrange(BD.genome.chrmCount):
            for j in xrange(BD.genome.chrmCount):
                beg1 = BD.chromosomeStarts[i]
                beg2 = BD.chromosomeStarts[j]
                end1 = BD.chromosomeEnds[i]
                end2 = BD.chromosomeEnds[j]
                mymask = mask[beg1:end1, beg2:end2]
                d = data[beg1:end1, beg2:end2]
                toret = tanaySmooth(d) / tanaySmooth(mymask)
                toret[mymask == 0] = 0
                ret[beg1:end1, beg2:end2] = toret
        return ret
        

    BD.removePoorRegions(cutoff=2)

    BD.removeCis()

    BD.iterativeCorrectWithoutSS()
    data1 = BD.dataDict[experiment1]
    data2 = BD.dataDict[experiment2]

    mask = (numpy.sum(data1, axis=0) > 0) * (numpy.sum(data2, axis=0) > 0)
    validMask = mask[:, None] * mask[None, :]
    transmask = BD.chromosomeIndex[:, None] != BD.chromosomeIndex[None, :]
    cormask = transmask * validMask

    d1 = propagateSmooth(data1)
    d2 = propagateSmooth(data2)
    print scipy.stats.spearmanr(d1[cormask], d2[cormask])

def correctedScalingPlot(resolution, filename1, experiment1, genome, mouse=False, **kwargs):
    "Paper figure to compare scaling before/after correction"
    
   	global pp
    if (options.verbose):
        print >> sys.stdout, "correctedScalingPlot: res: %d file1: %s exp1:%s gen:%s" % (resolution, filename1, experiment1, genome)

    #plt.figure(figsize=(4, 4))
    Tanay = binnedDataAnalysis(resolution, genome)
    Tanay.simpleLoad(filename1, experiment1)
    Tanay.removePoorRegions()
    Tanay.removeDiagonal()
    Tanay.plotScaling(experiment1, label="Raw data", color="#A7A241")
    Tanay.iterativeCorrectWithSS()
    Tanay.plotScaling(experiment1, label="Corrected", color="#344370")
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
    
def compareInterarmMaps(resolution, filename1, filename2, experiment1, experiment2, genome, mouse=False, **kwargs):
    "plots witn 8 inetrarm maps - paper supplement figure"
	global pp
	
    if (options.verbose):
        print >> sys.stdout, "compareInterarmMaps: res: %d file1: %s file2: %s exp1:%s exp2:%s gen:%s" % (resolution, filename1, filename2, experiment1, experiment2, genome)

    Tanay = binnedDataAnalysis(resolution, genome)

    Tanay.simpleLoad(filename1, experiment1)
    Tanay.simpleLoad(filename2, experiment2)
    Tanay.removeDiagonal()
    Tanay.removePoorRegions(cutoff=2)
    #Tanay.removeStandalone(3)
    fs = 10
    vmin = None
    vmax = None

    plt.subplot(421)
    plt.title(experiment1+", raw", fontsize=fs)
    Tanay.averageTransMap(experiment1, vmin=vmin, vmax=vmax)
    plt.colorbar()
    plt.subplot(422)
    plt.title(experiment2+", raw", fontsize=fs)
    Tanay.averageTransMap(experiment2, vmin=vmin, vmax=vmax)
    plt.colorbar()

    Tanay.iterativeCorrectWithSS()
    vmin = None
    vmax = None

    plt.subplot(425)

    plt.title(experiment1+", with SS reads", fontsize=fs)
    Tanay.averageTransMap(experiment1, vmin=vmin, vmax=vmax)
    plt.colorbar()
    plt.subplot(426)

    plt.title(experiment2+", with SS reads", fontsize=fs)
    Tanay.averageTransMap(experiment2, vmin=vmin, vmax=vmax)
    plt.colorbar()

    Tanay.iterativeCorrectWithoutSS()
    vmin = None
    vmax = None
    plt.subplot(423)
    plt.title(experiment1+", no SS reads", fontsize=fs)
    Tanay.averageTransMap(experiment2, vmin=vmin, vmax=vmax)
    plt.colorbar()

    plt.subplot(424)
    plt.title(experiment2+", no ss reads", fontsize=fs)
    Tanay.averageTransMap(experiment2, vmin=vmin, vmax=vmax)
    plt.colorbar()
    Tanay.fakeCis()

    vmin = None
    vmax = None
    plt.subplot(427)

    plt.title(experiment1+", trans only", fontsize=fs)
    Tanay.averageTransMap(experiment1, vmin=vmin, vmax=vmax)
    plt.colorbar()
    plt.subplot(428)

    plt.title(experiment1+", trans only", fontsize=fs)
    Tanay.averageTransMap(experiment2, vmin=vmin, vmax=vmax)
    plt.colorbar()

    plt.show()
   	pp.savefig()
    
def compareCorrelationOfEigenvectors(resolution, filename1, filename2, experiment1, experiment2, genome, mouse=False, **kwargs):
	"""Plot correlation figure with eigenvector correlation between datasets
	paper figure """
	global pp
	if (options.verbose):
		print >> sys.stdout, "compareCorrelationOfEigenvectors: res: %d file1: %s file2: %s exp1:%s exp2:%s gen:%s" % (resolution, filename1, filename2, experiment1, experiment2, genome)

	Tanay = binnedDataAnalysis(resolution, genome)
	Tanay.simpleLoad(filename1, experiment1)
	Tanay.simpleLoad(filename2, experiment2)
	
	Tanay.removeDiagonal()
	Tanay.removePoorRegions()
	Tanay.removeZeros()
	Tanay.truncTrans()
	Tanay.fakeCis()
	M = 10
	Tanay.doEig(numPCs=M)
	
	
	E1 = Tanay.EigDict[experiment1]
	E2 = Tanay.EigDict[experiment2]
	
	data = numpy.zeros((M, M))
	
	for i in xrange(M):
	    for j in xrange(M):
	        data[i][j] = abs(numpy.corrcoef(E2[i], E1[j])[0, 1])
	
	plt.figure(figsize=(7.5, 2.5))
	plt.gcf().subplots_adjust(0.2, 0.2, 0.85, 0.85)
	plt.subplot(111)
	plt.xlabel(experiment1)
	plt.ylabel(experiment2)
	#plt.title("Abs. correlation between eigenvectors")
	plt.imshow(data, interpolation="nearest", vmin=0, vmax=1)
	plt.colorbar()
	plt.show()
	pp.savefig()
    
def plotDiagonalCorrelation(resolution, filename1, filename2, experiment1, experiment2, genome, mouse=False, **kwargs):
    "Correlation of diagonal bins - paper figure"
	global pp

    if (options.verbose):
        print >> sys.stdout, "plotDiagonalCorrelation: res: %d file1: %s file2: %s exp1:%s exp2:%s gen:%s" % (resolution, filename1, filename2, experiment1, experiment2, genome)

    S = 50
    x = numpy.arange(2, S)
    Tanay = binnedData(resolution, genome)
    Tanay.simpleLoad(filename1, experiment1)
    Tanay.simpleLoad(filename2, experiment2)
    Tanay.removeDiagonal(1)
    Tanay.removePoorRegions()
    Tanay.removeZeros()

    pairs = [(experiment1, experiment2)]
    
    cors = [[] for _ in pairs]
    for i in x:
        for j, pair in enumerate(pairs):
            cors[j].append(cr(
                           numpy.diagonal(Tanay.dataDict[pair[0]], i),
                           numpy.diagonal(Tanay.dataDict[pair[1]], i)
                           )[0])

    Tanay.iterativeCorrectWithoutSS(M=1)
    cors2 = [[] for _ in pairs]
    for i in x:
        for j, pair in enumerate(pairs):
            cors2[j].append(cr(
                            numpy.diagonal(Tanay.dataDict[pair[0]], i),
                            numpy.diagonal(Tanay.dataDict[pair[1]], i)
                            )[0])
    Tanay.iterativeCorrectWithoutSS(M=20)
    cors3 = [[] for _ in pairs]
    for i in x:
        for j, pair in enumerate(pairs):
            cors3[j].append(cr(
                            numpy.diagonal(Tanay.dataDict[pair[0]], i),
                            numpy.diagonal(Tanay.dataDict[pair[1]], i)
                            )[0])

    matplotlib.rcParams['font.sans-serif'] = 'Arial'

    #plt.figure(figsize = (2.3,1.8))
    print cors
    print cors2
    print cors3
    plt.figure(figsize=(10, 3))
    ax = plt.gca()
    for j, pair in enumerate(pairs):
        plt.subplot(1, len(pairs), j)
        fs = 8
        for xlabel_i in ax.get_xticklabels():
            xlabel_i.set_fontsize(fs)
        for xlabel_i in ax.get_yticklabels():
            xlabel_i.set_fontsize(fs)
        plt.title("%s vs %s" % pair)
        plt.plot(x / 5., cors3[j], color="#E5A826", label="Iterative")
        plt.plot(x / 5., cors2[j], color="#28459A", label="Single")
        plt.plot(x / 5., cors[j], color="#E55726", label="Raw")
        plt.xlabel("Genomic Separation, MB", fontsize=8)
        plt.ylabel("Spearman correlation", fontsize=8)
        plt.legend()

        legend = plt.legend(prop={"size": 6}, loc=9, handlelength=2)
        legend.draw_frame(False)
        plt.ylim((0, 1))
        removeAxes(shift=0)

    plt.show()
	pp.savefig()
    
def process():
	global options
	global args
	global pp 
	
	fig = plt.gcf()
	pp = PdfPages(options.outputDir+'HiC-correlate.pdf')
	
	experiments = options.experiment.split(',')
	enzymes = options.enzyme.split(',')
	if (len(enzymes) == 1):
		enzymes = enzymes * len(args) 

	# check dataset exist
	for i in xrange(len(args)):
		if (not os.path.isfile(args[i] + os.sep + enzymes[i] +'_' + experiments[i] + '-1M.hdf5')):
			print '[ERROR] Could not find: '+ args[i] + os.sep + enzymes[i] +'_' + experiments[i] + '-1M.hdf5'
			sys.exit(1)
			
		if (not os.path.isfile(args[i] + os.sep + enzymes[i] +'_' + experiments[i] + '-200k.hdf5')):
			print '[ERROR] Could not find: '+ args[i] + os.sep + enzymes[i] +'_' + experiments[i] + '-1M.hdf5'
			sys.exit(1)
			
		if (not os.path.isfile(args[i] + os.sep + enzymes[i] +'_' + experiments[i] + '-IC-1M.hdf5')):
			print '[ERROR] Could not find: '+ args[i] + os.sep + enzymes[i] +'_' + experiments[i] + '-IC-1M.hdf5'
			sys.exit(1)

	genome_db    = genome.Genome(options.genome, gapFile=options.gapFile, readChrms=['#', 'X', 'Y'])
	
	for i in xrange(len(args)):
		print " Process file "+str(i)+":"+	args[i] + os.sep + enzymes[i] +'_' + experiments[i]
		
		correctedScalingPlot(200000, args[i] + os.sep + enzymes[i] +'_' + experiments[i] + '-200k.hdf5', experiments[i], genome_db)

		doArmPlot(1000000, args[i] + os.sep + enzymes[i] +'_' + experiments[i]+ '-1M.hdf5', experiments[i], genome_db)

		for j in xrange(i+1, len(args)):
		
			compareCorrelationOfEigenvectors(1000000, args[i] + os.sep + enzymes[i] +'_' +experiments[i] + '-1M.hdf5', args[j] + os.sep + enzymes[i] +'_' +experiments[j] + '-1M.hdf5', experiments[i], experiments[j], genome_db)

			calculateTanayCorrelation(1000000, args[i] + os.sep + enzymes[i] +'_' + experiments[i] + '-1M.hdf5', args[j] + os.sep + enzymes[i] +'_' +experiments[j] + '-1M.hdf5', experiments[i], experiments[j], genome_db)

			plotDiagonalCorrelation(200000, args[i] + os.sep + enzymes[i] +'_' + experiments[i] + '-200k.hdf5', args[j] + os.sep + enzymes[i] +'_' +experiments[j] + '-200k.hdf5', experiments[i], experiments[j], genome_db)
						
			compareInterarmMaps(1000000, args[i] + os.sep + enzymes[i] +'_' + experiments[i] + '-1M.hdf5', args[j] + os.sep + enzymes[i] +'_' +experiments[j] + '-1M.hdf5', experiments[i], experiments[j], genome_db)
			
	if (options.verbose):
		print >> sys.stdout, "print plots into pdf:%s" % (options.outputDir+'HiC-correlate.pdf')
	
	pp.close()
	
######################################
# main
######################################
if __name__ == "__main__":
	main()
