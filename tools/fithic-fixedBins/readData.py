# author Fabian Buske 2015
import os, sys, re
import traceback
from optparse import OptionParser
import fileinput
import datetime
from quicksect import IntervalTree
import gzip
import numpy as np
import collections
from functools import partial
from rosetta.parallel.parallel_easy import imap_easy
try:
   import cPickle as pickle
except:
   import pickle
######################################
# Read
######################################

class Read():
    def __init__(self, read):

        if (read==""):
            self.seq=""
            self.qname="dummy"
            self.is_unmapped=True
            self.tid=None
            self.is_read1=None
            self.is_reverse=None
        else:
            self.is_duplicate=read.is_duplicate
            self.is_unmapped=read.is_unmapped
            self.tid=read.tid
            self.qname=read.qname
            self.is_read1=True
            if (read.is_read2):
                self.is_read1=False
            self.is_reverse=read.is_reverse
            self.alen=read.alen
            self.pos=read.pos


    def check(self):
        if(self.is_reverse):
            self.revcomp()

    def isPair(self,number):
        if(self.is_read1 and number==0):
            return True
        if(not(self.is_read1) and number==1):
            return True
        print("Reads not paired up correctly: paired/single ended? not namesorted?")
        print(str(self)+" "+str(number))

    def __str__(self):
        if(self.is_unmapped):
            return "%s read1 %s %s" % (self.qname,self.is_read1, self.seq)
        else:
            return "%s read1 %s %s %i %i" % (self.qname,self.is_read1, self.tid, self.pos,self.alen)

    def revcomp(self):
        basecomplement = {'N':'N','A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        rc=""
        for i in reversed(self.seq):
            rc+=basecomplement[i]
        self.seq=rc


    def getInfo(self, string):
        for i in self.tags:
            if( i[0]==string):
                return i[1]

######################################
# Interval
######################################
class Interval():
    def __init__(self, chrom, start, end):
        self.chrom=chrom
        self.start=start
        self.end=end

######################################
# Timestamp
######################################
def timeStamp():
    return datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S').format()


def createIntervalTreesFragmentResolution(options):
    '''
    creates lookup tables
    returns
        fragmentsMap[fragmentId] = [tuple(chrom, start, end)]
        chromBinMap - lookup table for fixed bin sizes
    '''
    if (options.verbose):
        print >> sys.stdout, "- %s START   : populate lookup table with given resolution for chromosomes matching pattern %s" % (timeStamp(), options.chromPattern)

    fragmentsCount = 0
    fragmentsMap = {}
    fragmentsChrom = {}
    chromBinMap = {}

    for line in fileinput.input([options.chromSizes]):
        chrom=line.split("\t")[0]
        # check if chromosome needs to be filtered out or not
        if (options.chromPattern != "" and not re.match("^"+options.chromPattern+"$", chrom)):
            # skip this one
            if (options.vverbose):
                print "skipping pattern %s" % (line)
            continue

        fragmentsStart=fragmentsCount
        chromlen=int(line.split("\t")[1])

        for start in range(0, chromlen, options.resolution):
            end=min(int(start+options.resolution), chromlen)
            chromBinMap[tuple([chrom,int(start/options.resolution)])] = fragmentsCount
            fragmentsMap[fragmentsCount] = tuple([chrom, start, end])
            fragmentsCount += 1
            if (options.vverbose):
                print >> sys.stdout, "-- chromBinMap.add %s:%d (%d)" % (chrom, start, int(start/options.resolution) )

        fragmentsEnd=fragmentsCount
        fragmentsChrom[chrom] = tuple([fragmentsStart, fragmentsEnd])
    if (options.verbose):
        print >> sys.stdout, "- %s FINISHED: lookup table populated" % (timeStamp())

    return [ fragmentsMap, chromBinMap, fragmentsCount, fragmentsChrom ]

def createIntervalTreesFragmentFile(options):
    '''
        creates one interval tree for quick lookups
        returns
            fragmentsMap[fragmentId] = [tuple(chrom, start, end)]
            intersect_tree - intersect Tree for interval matching

    '''

    if (options.verbose):
        print >> sys.stdout, "- %s START   : populate intervaltree from fragmented genome" % (timeStamp())

    intersect_tree = IntervalTree()
    fragmentsCount = 0
    fragmentsMap = {}
    fragmentsChrom = {} # lookp table for fragment ranges of a chromosome
    fragmentsStart = 0

    start = 0
    end = 0
    counter = 0
    chrom = ""

    for line in fileinput.input([options.fragmentFile]):
        line = line.strip()
        if (len(line)==0 or line.startswith("Genome") or line.startswith("Chromosome")):
            continue

        cols = line.split("\t")
        try:
            # check if chromosome changed from last
            if (cols[0] != chrom):
                # do we have to finish the last chromosome?
                if (end > 0):
                    interval = Interval(chrom, start, end)
                    intersect_tree.insert(interval, fragmentsCount)
                    fragmentsMap[fragmentsCount] = tuple([chrom, start, end])
                    fragmentsCount += 1

                    fragmentsChrom[chrom] = tuple([fragmentsStart, fragmentsCount])
                    fragmentsStart = fragmentsCount

                    if (options.vverbose):
                        print >> sys.stdout,  "-- intervaltree.add %s:%d-%d" % (chrom, start, end)
                # check if chromosome needs to be filtered out or not
                if (options.chromPattern  != "" and not re.match(options.chromPattern, cols[0])):
                    chrom = ""
                    start = 0
                    end = 0

                else:
                    chrom = cols[0]
                    start = int(cols[1])
                    end = int(cols[2])
                counter = 0

            # check if fragment aggregation is fulfilled
            elif (counter >= options.fragmentAggregation):
                interval = Interval(chrom, start, end)
                intersect_tree.insert(interval, fragmentsCount)
                if (options.vverbose):
                    print >> sys.stdout,  "-- intervaltree.add %s:%d-%d" % (chrom, start, end)

                fragmentsMap[fragmentsCount] = tuple([chrom, start, end])
                start = int(cols[1])
                end = int(cols[2])
                counter = 0
                fragmentsCount += 1
            else:
                end = int(cols[2])

            # increment counter
            counter += 1

        except:
            if (options.verbose):
                print >> sys.stderr, 'skipping line in options.fragmentFile: %s' % (line)
            if (options.vverbose):
                traceback.print_exc()
                sys.exit(1)


    # handle last fragment
    if (end > 0):
        interval = Interval(chrom, start, end)
        intersect_tree.insert(interval, fragmentsCount)
        fragmentsMap[fragmentsCount] = tuple([chrom, start, end])
        fragmentsCount += 1
        fragmentsChrom[chrom] = tuple([fragmentsStart, fragmentsCount])

        if (options.vverbose):
            print >> sys.stdout, "-- intervaltree.add %s:%d-%d" % (chrom, start, end)

    if (options.verbose):
        print >> sys.stdout, "- %s FINISHED: intervaltree populated" % (timeStamp())

    return [fragmentsMap, intersect_tree, fragmentsCount, fragmentsChrom]

def getNext(iterator, options):
    '''
    get the next read and populate object
    '''

    try:
        return Read(iterator.next())
    except StopIteration:
        if (options.vverbose):
            traceback.print_exc()
        return Read("")


def findNextReadPair(samiter,options):
    '''
    find the next read pair
    '''

    readpair=[getNext(samiter, options), getNext(samiter, options)]

    while( re.sub('\/[1,2]$', '', readpair[0].qname) != re.sub('\/[1,2]$', '', readpair[1].qname)):
        if (options.verbose):
            print "[WARN] File: drop first from unpaired read %s %s" %(readpair[0].qname, readpair[1].qname)
        readpair.pop(0)
        readpair.append(getNext(samiter))

    return readpair

def find(interval, tree):
    ''' Returns a list with the overlapping intervals '''
    out = []
    tree.intersect( interval, lambda x: out.append(x) )
    return [ (x.start, x.end, x.linenum) for x in out ]


def getFragment(inputfile, read, lookup_structure, options):
    ''' When input is bam file '''

    fragmentID = None
    try:
        # get fragments for both reads
        rchrom = inputfile.getrname(read.tid)
        rstart = read.pos
        rend =  read.pos+read.alen

        if (options.vverbose):
            print >> sys.stdout, "- Check   : read %s %d %d" % (rchrom, rstart, rend )

        if (options.fragmentFile):
            interval = Interval(rchrom, rstart, rstart+1)
            # take first bin only
            fragmentID = find(interval, lookup_structure)[0][2]

        else:
            # for fixed bin sizes
            try:
                fragmentID = lookup_structure[tuple([rchrom, int(rstart/options.resolution)])]
            except:
                if (options.vverbose):
                    print >> sys.stderr, '[WARN] not in lookup : %s %d %d(skipping)' % (rchrom, rstart, int(rstart/options.resolution))
                return None

        if (fragmentID == None):
            return fragmentID

    except:
        if (options.verbose):
            print >> sys.stderr, '[WARN] problems with interval intersection: %s (skipping)' % (read)
            traceback.print_exc()
            sys.exit(1)
        if (options.vverbose):
            traceback.print_exc()
            sys.exit(1)

    return fragmentID

def mapFragment(rchrom, rstart, lookup_structure, options):
    ''' When input is fragment pair file '''

    fragmentID = None
    try:
        # get fragments for both reads

        if (options.vverbose):
            print >> sys.stdout, "- Check   : fragment %s %d " % (rchrom, rstart )

        if (options.fragmentFile):
            interval = Interval(rchrom, rstart, rstart+1)
            # take first bin only
            fragmentID = find(interval, lookup_structure)[0][2]

        else:
            # for fixed bin sizes
            try:
                fragmentID = lookup_structure[tuple([rchrom, int(rstart/options.resolution)])]
            except:
                if (options.vverbose):
                    print >> sys.stderr, '[WARN] not in lookup : %s %d (skipping)' % (rchrom, rstart)
                return None

        if (fragmentID == None):
            return fragmentID

    except:
        if (options.verbose):
            print >> sys.stderr, '[WARN] problems with interval intersection: %s %d (skipping)' % (rchrom, rstart)
            traceback.print_exc()
            sys.exit(1)
        if (options.vverbose):
            traceback.print_exc()
            sys.exit(1)

    return fragmentID

def populateFragmentPairs(fragmentPairs, x, count):
    f_tuple = tuple([x.minFragment, x.maxFragment])
    if (not fragmentPairs.has_key(f_tuple)):
        fragmentPairs[f_tuple] = 0
    fragmentPairs[f_tuple] += count

def countReadsPerFragmentSerial(lookup_structure, options, args):
    '''
        counts the reads per fragment and generates appropriate output files
    '''

    fragmentList = collections.defaultdict(int)
    fragmentPairs = collections.defaultdict(int)

    if (options.inputIsFragmentPairs):
        for fFile in xrange(len(args)):
            if (options.verbose):
                print >> sys.stdout, "- %s START   : processing fragment files: %s" % (timeStamp(), args[fFile])

            with gzip.open(args[fFile]) as infile:
                for line in infile:
                    (chr1, start1, chr2, start2, count) = line.strip().split("\t")

                    # skip trans contacts if focusing on cis only
                    if (options.onlycis and chr1 != chr2):
                        continue

                    fragmentID1 = mapFragment(chr1, int(start1), lookup_structure, options)
                    fragmentID2 = mapFragment(chr2, int(start2), lookup_structure, options)

                    if (fragmentID1 == None or fragmentID2 == None):
                        if (options.vverbose):
                            print >> sys.stdout, "-- one region does not co-occur with any fragment: %s %s" % (fragmentID1, fragmentID2)
                        continue

                    fragmentList[fragmentID1] += 1
                    fragmentList[fragmentID2] += 1

                    f_tuple = tuple([min(fragmentID1, fragmentID2), max(fragmentID1, fragmentID2)])
                    fragmentPairs[f_tuple] += int(count)

            if (options.verbose):
                print >> sys.stdout, "- %s FINISHED: getting counts form fragment files " % (timeStamp())

    elif (options.inputIsReadPairs != ""):
       # get column indexes
       chr_index = map(int,  options.inputIsReadPairs.split(",")[0:4])
       try:
            chr_prefix=options.inputIsReadPairs.split(",")[4]
       except:
            chr_prefix=""

       for fFile in xrange(len(args)):
            if (options.verbose):
                print >> sys.stdout, "- %s START   : processing read files: %s" % (timeStamp(), args[fFile])

            with gzip.open(args[fFile]) as infile:
                for line in infile:
                    cols   = line.strip().split()
                    chr1   = chr_prefix + cols[chr_index[0]]
                    start1 = int(cols[chr_index[1]])
                    chr2   = chr_prefix + cols[chr_index[2]]
                    start2 = int(cols[chr_index[3]])

                    # skip trans contacts if focusing on cis only
                    if (options.onlycis and chr1 != chr2):
                        continue

                    fragmentID1 = mapFragment(chr1, start1, lookup_structure, options)
                    fragmentID2 = mapFragment(chr2, start2, lookup_structure, options)

                    if (fragmentID1 == None or fragmentID2 == None):
                        if (options.vverbose):
                            print >> sys.stdout, "-- one region does not co-occur with any fragment: %s %s" % (fragmentID1, fragmentID2)
                        continue

                    fragmentList[fragmentID1] += 1
                    fragmentList[fragmentID2] += 1
                    f_tuple = tuple([min(fragmentID1, fragmentID2), max(fragmentID1, fragmentID2)])
                    fragmentPairs[f_tuple] += 1

            if (options.verbose):
                print >> sys.stdout, "- %s FINISHED: getting counts form read files " % (timeStamp())

    else:
        # lazy load
        import pysam
        for bamFile in xrange(len(args)):
            if (options.verbose):
                print >> sys.stdout, "- %s START   : processing reads from bam file: %s" % (timeStamp(), args[bamFile])

            samfile = pysam.Samfile(args[bamFile], "rb" )
            samiter = samfile.fetch(until_eof=True)
            readcounter = 0

            while(True):
                readpair = findNextReadPair(samiter,options)
                # if file contains any more reads, exit
                if (readpair[0].qname=="dummy"):
                    break

                # skip trans contacts if focusing on cis only
                if (options.onlycis and inputfile.getrname(readpair[0].tid) != inputfile.getrname(readpair[1].tid)):
                    continue

                fragmentID1 = getFragment(samfile, readpair[0], lookup_structure, options)
                fragmentID2 = getFragment(samfile, readpair[1], lookup_structure, options)

                if (fragmentID1 == None or fragmentID2 == None):
                    if (options.vverbose):
                        print >> sys.stdout, "-- one read does not co-occur with any fragment: %s %s" % (str(fragmentID1), str(fragmentID2))
                    continue
                fragmentList[fragmentID1] += 1
                fragmentList[fragmentID2] += 1
                f_tuple = tuple([min(fragmentID1, fragmentID2), max(fragmentID1, fragmentID2)])
                fragmentPairs[f_tuple] += 1
                readcounter+=1

                if (options.verbose and readcounter % 1000000 == 0 ):
                    print >> sys.stdout, "- %s         : %d read pairs processed" % (timeStamp(), readcounter)
            samfile.close()

            if (options.verbose):
                print >> sys.stdout, "- %s FINISHED: getting reads from bam file " % (timeStamp())

    return [ fragmentList.items(), fragmentPairs.items() ]

def countReadsPerFragmentParallel(iFile, lookup_structure, options):
    '''
        counts the reads per fragment and generates appropriate output files
    '''

    fragmentList = collections.defaultdict(int)
    fragmentPairs = collections.defaultdict(int)

    if (options.inputIsFragmentPairs):
        if (options.verbose):
            print >> sys.stdout, "- %s START   : processing fragment file: %s" % (timeStamp(), iFile)

        with gzip.open(iFile) as infile:
            for line in infile:
                (chr1, start1, chr2, start2, count) = line.strip().split("\t")

                # skip trans contacts if focusing on cis only
                if (options.onlycis and chr1 != chr2):
                    continue

                fragmentID1 = mapFragment(chr1, int(start1), lookup_structure, options)
                fragmentID2 = mapFragment(chr2, int(start2), lookup_structure, options)

                if (fragmentID1 == None or fragmentID2 == None):
                    if (options.vverbose):
                        print >> sys.stdout, "-- one region does not co-occur with any fragment: %s %s" % (fragmentID1, fragmentID2)
                    continue

                fragmentList[fragmentID1] += 1
                fragmentList[fragmentID2] += 1
                f_tuple = tuple([min(fragmentID1, fragmentID2), max(fragmentID1, fragmentID2)])
                fragmentPairs[f_tuple] += int(count)


        if (options.verbose):
            print >> sys.stdout, "- %s FINISHED: getting counts form fragment file " % (timeStamp())

    elif (options.inputIsReadPairs != ""):
        # get column indexes
        chr_index = map(int,  options.inputIsReadPairs.split(",")[0:4])
        try:
            chr_prefix=options.inputIsReadPairs.split(",")[4]
        except:
            chr_prefix=""

        if (options.verbose):
            print >> sys.stdout, "- %s START   : processing read files: %s" % (timeStamp(), iFile)

        with gzip.open(iFile) as infile:
            for line in infile:
                cols   = line.strip().split()
                chr1   = chr_prefix + cols[chr_index[0]]
                start1 = int(cols[chr_index[1]])
                chr2   = chr_prefix + cols[chr_index[2]]
                start2 = int(cols[chr_index[3]])

                # skip trans contacts if focusing on cis only
                if (options.onlycis and chr1 != chr2):
                    continue

                fragmentID1 = mapFragment(chr1, start1, lookup_structure, options)
                fragmentID2 = mapFragment(chr2, start2, lookup_structure, options)

                if (fragmentID1 == None or fragmentID2 == None):
                    if (options.vverbose):
                        print >> sys.stdout, "-- one region does not co-occur with any fragment: %s %s" % (fragmentID1, fragmentID2)
                    continue

                fragmentList[fragmentID1] += 1
                fragmentList[fragmentID2] += 1
                f_tuple = tuple([min(fragmentID1, fragmentID2), max(fragmentID1, fragmentID2)])
                fragmentPairs[f_tuple] += 1

        if (options.verbose):
            print >> sys.stdout, "- %s FINISHED: getting counts form read files " % (timeStamp())

    else:
        # lazy load
        import pysam
        if (options.verbose):
            print >> sys.stdout, "- %s START   : processing reads from bam file: %s" % (timeStamp(), iFile)

        samfile = pysam.Samfile(iFile, "rb" )
        samiter = samfile.fetch(until_eof=True)
        readcounter = 0

        while(True):
            readpair = findNextReadPair(samiter,options)
            # if file contains any more reads, exit
            if (readpair[0].qname=="dummy"):
                break

            # skip trans contacts if focusing on cis only
            if (options.onlycis and inputfile.getrname(readpair[0].tid) != inputfile.getrname(readpair[1].tid)):
                continue

            fragmentID1 = getFragment(samfile, readpair[0], lookup_structure, options)
            fragmentID2 = getFragment(samfile, readpair[1], lookup_structure, options)

            if (fragmentID1 == None or fragmentID2 == None):
                if (options.vverbose):
                    print >> sys.stdout, "-- one read does not co-occur with any fragment: %s %s" % (str(fragmentID1), str(fragmentID2))
                continue

            fragmentList[fragmentID1] += 1
            fragmentList[fragmentID2] += 1
            f_tuple = tuple([min(fragmentID1, fragmentID2), max(fragmentID1, fragmentID2)])
            fragmentPairs[f_tuple] += 1
            readcounter+=1

            if (options.verbose and readcounter % 1000000 == 0 ):
                print >> sys.stdout, "- %s         : %d read pairs processed" % (timeStamp(), readcounter)
        samfile.close()

        if (options.verbose):
            print >> sys.stdout, "- %s FINISHED: getting reads from bam file " % (timeStamp())

    # dump file to disc as big objects cannot be transfered through the multi-processing frameworks
    fn = iFile+".pickle"

    pickle.dump(tuple([fragmentList, fragmentPairs ]), open(fn, 'wb'))
    return fn


def countReadsPerFragment(lookup_structure, options, args):
    '''
        slurps in all input fils in parallel
        counts the reads per fragment and generates appropriate output files
    '''

    if (options.verbose):
        print >> sys.stdout, "- %s STARTED : reading input files : %s" % (timeStamp(), str(args))

    fragmentList = collections.defaultdict(int)
    fragmentPairs = collections.defaultdict(int)
    func = partial(countReadsPerFragmentParallel, lookup_structure=lookup_structure, options=options)
    results_iterator = imap_easy(func, args, n_jobs=8, chunksize=1)

    if (options.verbose):
        print >> sys.stdout, "- %s FINISHED: reading input files " % (timeStamp())
        print >> sys.stdout, "- %s STARTED : combining input files " % (timeStamp())


    # combine dictionaries
    for fn in results_iterator:
        if (options.verbose):
            print >> sys.stdout, "- reading pickled file %s" % (fn)

	(fl,fp) = pickle.load(open(fn,'rb'))
        for k, v in fl.iteritems():
            fragmentList[k]+=v

        for k, v in fp.iteritems():
            fragmentPairs[k]+=v
	# delete file
	os.remove(fn)

    if (options.verbose):
        print >> sys.stdout, "- %s FINISHED: combining input files " % (timeStamp())


    return [ fragmentList, fragmentPairs ]
