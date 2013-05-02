from optparse import OptionParser
import sys, os, getopt

bundle=50

random="/clusterdata/hiseq_apps/hiSeqInf/tools/shufSTDIN.py"
picard="/clusterdata/apps/picard-tools-1.31/"
samtools="/clusterdata/hiseq_apps/tools/freeze001/samtools"
TMP="/clusterdata/denis/tmp/"

#qout="qout/downsample/"

# bundle the commands together in a file
def add (com, command, outdir, name, counter, last, opt):
    command.append(com)
    if (len(command)>bundle or last):
        out=open(outdir+name+str(counter),"w")
        out.write("\n".join(command))
        out.write("\nrm "+outdir+name+str(counter)+"\n")
        os.system("chmod -u=rwx "+outdir+name+str(counter))
        os.system("qsub "+opt+" -cwd -b y -j y -o "+options.qout+name+str(counter)+".out -N "+name+str(counter)+" ./"+outdir+name+str(counter))
        command=[]
        counter+=1
    return (command, counter)

# split in little windows
def split(options):
    
    b=options.begin
    e=options.end
    w=options.window

    outname=options.input.split("/")[-1].split(".bam")[0]
    collection=[]
    counter=0
    c2=0
    filelist=[]

    while( b+w<e ):
        outname="%s_%i_%i" % (options.input.split("/")[-1].split(".bam")[0],options.sample, b)
        com1="%s view %s %s:%i-%i | python %s -n %i > %s.body.tmp" % (samtools, options.input,options.chr,b,int(b+w-(w*0.25)),random,options.sample,options.tmp+outname)
        com2="%s view -H %s | cat - %s.body.tmp | samtools view -S -b - | %s sort - %s" % (samtools,options.input,options.tmp+outname,samtools,options.tmp+outname)
        com3="rm "+options.tmp+outname+".body.tmp"
        command= "\n".join([com1,com2,com3])
        #print command
        (collection,counter)=add(command, collection, options.tmp, "down"+options.input.split("/")[-1],counter, False, "")
        filelist.append(options.tmp+outname+".bam")
        c2+=1
        b=int(b+w-(w*0.25))

    if(b<e):
        sample=((e-b)/float(w))*options.sample
        print "%i %i %i" % (b,e,sample)
        outname="%s_%i_%i" % (options.input.split("/")[-1].split(".bam")[0],options.sample, b)
        com1="%s view %s %s:%i-%i | python %s -n %i > %s.body.tmp" % (samtools,options.input,options.chr,b,e,random,sample,options.tmp+outname)
        com2="%s view -H %s | cat - %s.body.tmp | samtools view -S -b - | %s sort - %s" % (samtools,options.input,options.tmp+outname,samtools,options.tmp+outname)
        com3="rm "+options.tmp+outname+".body.tmp"
        command= "\n".join([com1,com2,com3])
        #print command
        (collection,counter)=add(command, collection, options.tmp, "down"+options.input.split("/")[-1],counter, True, "")
        c2+=1
        filelist.append(options.tmp+outname+".bam")


    #print "find "+str(c2)+" files"
    return filelist

# merge windows at the end
def merge(options,filelist):
    output=options.output+options.input.replace(".bam","."+str(options.sample)+"x.bam").split("/")[-1]
    com1="java -Xmx10g -jar "+picard+"MergeSamFiles.jar OUTPUT="+output+" USE_THREADING=true"

    outname=options.input.split("/")[-1].split(".bam")[0]+"*.bam"
    for f in filelist:
        if(f!="" and f.find("x.bam")==-1):
            com1+=" INPUT="+f

    #remove duplicates
    com2="java -Xmx10g -jar "+picard+"MarkDuplicates.jar INPUT="+output+" OUTPUT="+output.replace("x.bam","d.bam")
    com2+=" METRICS_FILE="+output+".dupl  AS=true VALIDATION_STRINGENCY=LENIENT TMP_DIR="+TMP

    # make index
    com3=samtools+" index "+output.replace("x.bam","d.bam")

    #remove files
    com4="rm "+" ".join(filelist)
    
    #print command
    command= "\n".join([com1,com2,com3,com4])
    #print command
    n=options.input.split("/")[-1].split("bam")[0]+str(options.sample)+"x.bam"
    add(command, [], options.tmp, "merge"+n ,0, True, "-pe mpich2_mpd 2 -hold_jid down"+outname.strip(".bam"))

parser = OptionParser()
   
parser.add_option("-i","--input", dest = "input", help = "input bam file")
parser.add_option("-o","--output", dest = "output", help = "output bam file")
parser.add_option("-t","--tmp", dest = "tmp", help = "tempdir for window files")
parser.add_option("-c","--chr", dest = "chr", help = "chromosome") 
parser.add_option("-b","--begin", dest = "begin", help = "beginning of the downsampling", type="int")
parser.add_option("-e","--end", dest = "end", help = "end of downsampling", type="int")
parser.add_option("-r","--region", dest = "region", help = "region chrX:Y-Z")
parser.add_option("-w","--window", dest = "window", help = "window to sample <cov> reads from", type="int")
parser.add_option("-s","--sample", dest = "sample", help = "how many reads to sample", type="int")
parser.add_option("-q","--qout", dest = "qout", help = "qsub out")
(options, args) = parser.parse_args() 

if (options.region!=""):
    options.chr=options.region.split(":")[0]
    options.begin=int(options.region.split(":")[1].split("-")[0])
    options.end=int(options.region.split(":")[1].split("-")[1])

filelist=[]
filelist=split(options)
merge(options,filelist)
