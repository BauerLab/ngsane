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
dataset1 should point to the "[ENZYME]-[EXPERIMENT]-fragment-dataset.hdf5" file
	'''
	parser = OptionParser(usage)
	parser.add_option("-q", "--quiet", action="store_false", dest="verbose", default=True,
					help="don't print status messages to stdout")
	parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
					help="print status messages to stdout")
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
	
	if (options.genome == ""):
		print >> sys.stderr, "[ERROR] Please specify the location of the reference genome in fasta format"
		sys.exit(1)

	if (options.outputDir != ""): 
		options.outputDir += os.sep
	
	if (options.verbose):
		print >> sys.stdout, "outputDir:          %s" % (options.outputDir)
		print >> sys.stdout, "tmpDir:             %s" % (options.tmpDir)

	process()

def calculateTanayCorrelation(resolution, filename1, filename2, experiment1, experiment2, genome, outfile, mouse=False, **kwargs):
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
    (scorr, pvalue) = scipy.stats.spearmanr(d1[cormask], d2[cormask])
    outfile.write("Spearman corrleation	%s	%s %.4f	%.4f" % (filename1, filename2, scorr, pvalue))

   
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

    plt.figure(figsize=(12, 16))
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
	
	plt.figure(figsize=(8, 8))
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

    print "Eigenvectors"
    print cors
    print cors2
    print cors3
    plt.figure(figsize=(10, 4))
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
	
	# check dataset exist
	for i in xrange(len(args)):
		if (not os.path.isfile(args[i].replace('-fragment_dataset.hdf5','-1M.hdf5'))):
			print '[ERROR] Could not find: '+ args[i].replace('-fragment_dataset.hdf5','-1M.hdf5')
			sys.exit(1)
			
		if (not os.path.isfile(args[i].replace('-fragment_dataset.hdf5','-200k.hdf5'))):
			print '[ERROR] Could not find: '+ args[i].replace('-fragment_dataset.hdf5','-200k.hdf5')
			sys.exit(1)
			
		if (not os.path.isfile(args[i].replace('-fragment_dataset.hdf5','-IC-1M.hdf5'))):
			print '[ERROR] Could not find: '+ args[i].replace('-fragment_dataset.hdf5','-IC-1M.hdf5')
			sys.exit(1)

	genome_db = genome.Genome(options.genome, gapFile=options.gapFile, readChrms=['#', 'X', 'Y'])
	
	outfile = open(options.outputDir+'HiC-correlate.txt',"w")
	for i in xrange(len(args)):
		print " Process file "+str(i)+":"+ args[i]
		enzyme_i = os.path.basename(args[i]).split("_")[0]
		experiment_i = "_".join(os.path.basename(args[i]).strip("-fragment_dataset.hdf5").split("_")[1:])
		for j in xrange(i+1, len(args)):
			enzyme_j = os.path.basename(args[j]).split("_")[0]
	                experiment_j = "_".join(os.path.basename(args[j]).strip("-fragment_dataset.hdf5").split("_")[1:])

			compareCorrelationOfEigenvectors(1000000, args[i].replace('-fragment_dataset.hdf5','-1M.hdf5'), args[j].replace('-fragment_dataset.hdf5','-1M.hdf5'), experiment_i, experiment_j, genome_db)

			calculateTanayCorrelation(1000000, args[i].replace('-fragment_dataset.hdf5','-1M.hdf5'), args[j].replace('-fragment_dataset.hdf5','-1M.hdf5'), experiment_i, experiment_j, genome_db, outfile)

			plotDiagonalCorrelation(200000, args[i].replace('-fragment_dataset.hdf5','-200k.hdf5'), args[j].replace('-fragment_dataset.hdf5','-200k.hdf5'), experiment_i, experiment_j, genome_db)
			
			compareInterarmMaps(1000000, args[i].replace('-fragment_dataset.hdf5','-1M.hdf5'), args[j].replace('-fragment_dataset.hdf5','-1M.hdf5'), experiment_i, experiment_j, genome_db)
			
	if (options.verbose):
		print >> sys.stdout, "print plots into pdf:%s" % (options.outputDir+'HiC-correlate.pdf')
	outfile.close()
	pp.close()
	
######################################
# main
######################################
if __name__ == "__main__":
	main()
