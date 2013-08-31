#!/bin/python
import sys
from optparse import OptionParser

RANGES = {
    'Sanger': (33, 73),
    'Solexa': (59, 104),
    'Illumina 1.3': (64, 104),
    'Illumina 1.5': (67, 104)
}


def get_qual_range(qual_str):
    """
    >>> get_qual_range("DLXYXXRXWYYTPMLUUQWTXTRSXSWMDMTRNDNSMJFJFFRMV")
    (68, 89)
    """

    vals = [ord(c) for c in qual_str]
    return min(vals), max(vals)

def process_lines():
    global options
    global args

    gmin, gmax  = 99, 0
    valid = []
    for i, line in enumerate(sys.stdin):
        lmin, lmax = get_qual_range(line.rstrip())
        if lmin < gmin or lmax > gmax:
            gmin, gmax = min(lmin, gmin), max(lmax, gmax)
            valid = get_encodings_in_range(gmin, gmax)
            if len(valid) == 0:
                print >>sys.stderr, "no encodings for range: %s" % str((gmin, gmax))
                sys.exit()
            if len(valid) == 1 and options.n == -1:
                print "\t".join(valid) + "\t" + str((gmin, gmax))
                sys.exit()

        if options.n > 0 and i > options.n:
            print "\t".join(valid) + "\t" + str((gmin, gmax))
            sys.exit()

    print "\t".join(valid) + "\t" + str((gmin, gmax))


def get_encodings_in_range(rmin, rmax, ranges=RANGES):
    valid_encodings = []
    for encoding, (emin, emax) in ranges.items():
        if rmin >= emin and rmax <= emax:
            valid_encodings.append(encoding)
    return valid_encodings

def main():
    global options
    global args
    usage = '''usage: %prog [options]
	guess the encoding of a stream of quality lines.

        Hint: extract quality lines of zipped fastq files with 
        zcat fastq.gz | awk 'NR % 4 == 0' | python %prog [options] 
        '''
    parser = OptionParser(usage)
    parser.add_option("-n", dest="n", help="number of qual lines to test default:-1"
                 " means test until end of file or until it it possible to "
                 " determine a single file-type",
                 type='int', default=-1)

    (options, args) = parser.parse_args()

    process_lines()

if __name__ == "__main__":
    main()
