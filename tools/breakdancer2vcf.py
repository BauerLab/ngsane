#!/usr/bin/env python

import csv, sys, optparse, re

parser = optparse.OptionParser()

parser.add_option("-i","--input", dest = "tsvfile", help = "input file")
parser.add_option("-o","--output", dest="vcffile", help = "vcf file")
(options, args) = parser.parse_args()

#16      27236526        2+0-    16      27237057        2+2-    DEL     568     79      2       RNA3kChr16rep.asd.bam|2  4.08
#16      27236526        .       .       .       .       PASS    PROGRAM=breakdancer;SVTYPE=DEL;SVLEN=-568;SVEND=27237057

outfile=open(options.vcffile,"w")

outfile.write("##fileformat=VCFv4.1\n")
outfile.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type\">\n")
outfile.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length\">\n")
outfile.write("##INFO=<ID=SVEND,Number=1,Type=Integer,Description=\"End\">\n")
outfile.write("##INFO=<ID=COV,Number=1,Type=String,Description=\"coverage\">\n")
outfile.write("#"+"\t".join(['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'])+"\n")

for i in open(options.tsvfile):
	if i.find("#Chr")>-1:
		samplename=i.split("\t")[-1]
		continue
	elif i.find("#")>-1:
		continue
	else:
		arr=re.split("[\t]",i)
		string="SVTYPE=%s;SVLEN=%i;SVEND=%s;COV=%s" % (arr[6],0-int(arr[7]),arr[4],arr[9])
		outa=[arr[0],arr[1],".",".",".",arr[9],"PASS",string]
        outfile.write("\t".join(outa)+"\n")
