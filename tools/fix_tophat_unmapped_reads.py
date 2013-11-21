#!/usr/bin/env python

# ----------------------------------------------------------------------------
#  "THE BEER-WARE LICENSE" (Revision 42):
#  Christian Brueffer <christian@brueffer.de> wrote this file. As long as you
#  retain this notice you can do whatever you want with this stuff. If we meet
#  some day, and you think this stuff is worth it, you can buy me a beer in
#  return.
# ----------------------------------------------------------------------------

# This script fixes the reads in a TopHat unmapped.bam to make them compatible
# with other tools, e.g. Picard and samtools.

import os
import pysam
import sys


def get_index_pos(index, unmapped_reads, read):
    """Returns the position of a read in the index or None."""

    if read.qname in index:
        return index[read.qname]
    else:
        return None


def main(path, outdir, mapped_file="accepted_hits.bam", unmapped_file="unmapped.bam"):
    bam_mapped = pysam.Samfile(os.path.join(path, mapped_file))
    bam_unmapped = pysam.Samfile(os.path.join(path, unmapped_file))
    unmapped_reads = list(bam_unmapped.fetch(until_eof=True))

    # Fix things that relate to all unmapped reads.
    unmapped_dict = {}
    unmapped_index = {}
    for i in range(len(unmapped_reads)):
        read = unmapped_reads[i]

        # remove /1 and /2 suffixes
        if read.qname.find("/") != -1:
            read.qname = read.qname[:-2]

        unmapped_index[read.qname] = i

        # work around "mate is unmapped" bug in TopHat
        if read.qname in unmapped_dict:
            unmapped_reads[unmapped_dict[read.qname]].mate_is_unmapped = True
            read.mate_is_unmapped = True
        else:
            unmapped_dict[read.qname] = i

        read.mapq = 0

        unmapped_reads[i] = read

    # Fix things that relate only to unmapped reads with a mapped mate.
    for mapped in bam_mapped:
        if mapped.mate_is_unmapped:
            i = get_index_pos(unmapped_index, unmapped_reads, mapped)
            if i is not None:
                unmapped = unmapped_reads[i]

                # map chromosome TIDs from mapped to unmapped file
                mapped_rname = bam_mapped.getrname(mapped.tid)
                unmapped_new_tid = bam_unmapped.gettid(mapped_rname)

                unmapped.tid = unmapped_new_tid
                unmapped.rnext = unmapped_new_tid
                unmapped.pos = mapped.pos
                unmapped.pnext = 0

                unmapped_reads[i] = unmapped

    bam_mapped.close()

    # for the output file, take the headers from the unmapped file
    base, ext = os.path.splitext(unmapped_file)
    out_filename = "".join([base, "_fixup", ext])
    bam_out = pysam.Samfile(os.path.join(outdir, out_filename), "wb",
                        template=bam_unmapped)

    bam_unmapped.close()

    for read in unmapped_reads:
        bam_out.write(read)

    bam_out.close()


def usage(scriptname):
    print "Usage:"
    print scriptname, "tophat_output_dir"
    sys.exit(1)


if __name__ == "__main__":
    if len(sys.argv) == 2:
        path = sys.argv[1]
        if os.path.exists(path) and os.path.isdir(path):
            # no tmpdir specified, use the bam dir
            main(path, path)
        else:
            usage(sys.argv[0])
    elif len(sys.argv) == 3:
        path = sys.argv[1]
        outdir = sys.argv[2]
        if os.path.exists(path) and os.path.isdir(path) and os.path.exists(outdir) and os.path.isdir(outdir):
            main(path, outdir)
        else:
            usage(sys.argv[0])
    else:
        usage(sys.argv[0])
