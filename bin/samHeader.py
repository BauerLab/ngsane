import pysam, sys

## date 3.August 2010
## author Denis C. Bauer



def usage():
    print """python samHeader.py -i <infile.bam> -o <outfile.bam> -r <reference.fai> [options]\n
    -a <assembly>  Assembly of the reference, default HG18\n
    -ID <name>     RG ID name, e.g. exeriment1, required\n
    -SM <name>     RG SM name, e.g. sample1, required\n
    -PL <name>     RG PL platform, default Illumna\n

    The script adds to an exising bam file the read group information.

    Note that infile needs to already know the reference sequence and needs \n
    to be sorted (by coordinate):
    \tsamtools view -bt reference.fai input.sam | samtools sort - input.sort"""
    sys.exit()

bam=""
out=""
reference=""
assembly="HG18"
ID=""
SM=""
PL="Illumina"

if(len(sys.argv)<6):
    print "too few parameters given"
    usage()

i=1
while(len(sys.argv)>i):
    if(sys.argv[i]=="-i"):
        i+=1
        bam=sys.argv[i]
    elif(sys.argv[i]=="-o"):
        i+=1
        out=sys.argv[i]
    elif(sys.argv[i]=="-r"):
        i+=1
        reference=sys.argv[i]
    elif(sys.argv[i]=="-a"):
        i+=1
        assembly=sys.argv[i]
    elif(sys.argv[i]=="-ID"):
        i+=1
        ID=sys.argv[i]
    elif(sys.argv[i]=="-SM"):
        i+=1
        SM=sys.argv[i]
    elif(sys.argv[i]=="-PL"):
        i+=1
        PL=sys.argv[i]
    else:
        print "don't understand "+sys.argv[i]
        usage()
    i+=1

if(bam=="" or out=="" or reference=="" or assembly=="" or ID=="" or SM=="" or PL==""):
    print "error: something is missing"
    print "bam "+bam
    print "out "+out
    print "assembly "+assembly
    print "ID "+ID
    print "SM "+SM
    print "PL "+PL
    usage()


# make the SQ tags for the different chromosomes in the .fai file
def makeSQTag(reference,assembly):
    ent=[]
    for f in open(reference):
        arr=f.split("\t")
        ent.append({'LN':int(arr[1]), 'SN':arr[0], 'AS':assembly})
    return ent

header = { 'HD': {'VN': '1.0', 'sorted':'coordinate'},
           'SQ': makeSQTag(reference,assembly),
           'RG': [{'ID': ID, 'SM': SM, 'PL':PL}] }

# open new bam file
outfile = pysam.Samfile( out, "wb", header = header )

# index old bam file -- in case it wasn't already and open it
pysam.index(bam)
infile = pysam.Samfile( bam, "rb" )

counter=0
# go through reads manipulate them and write them to the new bam file
# until_eof = True - important otherwise it just reads the mapped reads
for read in infile.fetch(until_eof = True):
    counter+=1
    # if they are unpapped create a read tag
    if(read.is_unmapped):
        read.tags=tuple([('RG',ID)])
    # if they are mapped just attach the RG info to it (RG:Z:ID)
    else:
        read.tags=tuple(list(read.tags)+[('RG',ID)])
    outfile.write(read)

print counter

infile.close()
outfile.close()
    
