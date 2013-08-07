#!/usr/bin/env

import argparse, sys, re

#
# Extract a specific (set of) feature(s) from a gtf file
# 
# Author Denis.Bauer@csiro.au


def walker(file,out,loc,delim,featuresKept,verbose,stdout):
    if (not(stdout)):
        out=open(out,"w")
    kept=0
    lines=0
    for i in open(file):
        lines+=1
        if i=="" or i[0]=="#":
            continue
        arr=re.split(delim,i)
        if (arr[loc].strip('";') in featuresKept):
            kept+=1
            if(not(stdout)):
                out.write(i)
            else:
                print i.strip("\n")
    if verbose:
        print "for feature %s %i were kept out of the %i processed lines" % (str(featuresKept),kept,lines)

def main():
    parser = argparse.ArgumentParser(description='Feature extractor')
    parser.add_argument('-f', type=str, help='GTF file input', required=True)
    parser.add_argument('-l', type=int, help='location of the feature name (GENCODE default 13)', default=13)
    parser.add_argument('-d', type=str, help='deliminator', default="[\t ]")
    parser.add_argument('-o', type=str, help='output (default stdout)', default=False)
    parser.add_argument('-v', action="store_true", help='verbose', default=False)
    parser.add_argument('--keep', type=str, nargs='+', help='element(s) to keep -- space separated list', required=True)
    args = parser.parse_args()

    stdout=True
    if args.o:
        stdout=False

    walker(args.f,args.o,args.l,args.d,set(args.keep),args.v,stdout)

if __name__ == '__main__': 
    main()
