#!/bin/python

import sys,os,math,re,traceback,datetime,glob

def removePrefix(text, prefix):
### remove prefix from a string ###
    return text[len(prefix):] if text.startswith(prefix) else text
   
def removeSuffix(text, suffix):
### remove prefix from a string ###
    return text[0:-len(suffix)] if text.endswith(suffix) else text
  
   
if (len(sys.argv)==1 or sys.argv[0].find("help")>-1):
    print "python2 times"
    die

dir=re.split("[ \n]",sys.argv[1])
ext=sys.argv[2]
type=sys.argv[3]
printing=True
percent=False
noSummary=False
noOverallSummary=False
overAll=True
link=""

if(dir[0]==""):
    dir.pop(0)


pseudocount=1e-20

#print "looking at "+str(dir)

i=3
# print minimal output only
while(len(sys.argv)>i):
    if(sys.argv[i]=="--essentials" or sys.argv[i]=="--e"):
        printing=False
    if(sys.argv[i]=="--percent" or sys.argv[i]=="--p"):
        percent=True
    if(sys.argv[i]=="--noOverAll" or sys.argv[i]=="--n"):
        overAll=False
    if(sys.argv[i]=="--link" or sys.argv[i]=="--l"):
		i+=1
		link=sys.argv[i]
    if(sys.argv[i]=="--noSummary" or sys.argv[i]=="--s"):
        noSummary=True
    if(sys.argv[i]=="--noOverallSummary" or sys.argv[i]=="--o"):
        noOverallSummary=True
    i+=1

if (dir == ["all"]):
    dir=[]
    for f in os.listdir('.'):
        if(f[-5::]=="_seqs"):
            if(type=="snpFilter"):
                dir.append(f+"/snp")
            else:
                dir.append(f+"/aln")
#if(printing):
#    print >>sys.stderr, dir

names=[]

def average(arr):
    sum=0.
    for a in arr:
        sum+=float(a)
    sum/=len(arr)
    return(sum)

# EXCEL comparable is =STDEV.P()
# because this is the population std
def std(arr,av):
    sum=0.
    for a in arr:
        sum+=(float(a)-float(av))*(float(a)-float(av))
    sum=math.sqrt(sum/(len(arr)))
    return(sum)

def ste(arr):
    av=average(arr)
    return(std(arr,av)/math.sqrt(len(arr)))

def per(max,arr):
    sum=0.
    for a in xrange(len(arr)):
        if (float(arr[a])!=0.0):
            sum+=float(arr[a])/(float(max[a])+pseudocount)*100
    if (sum==0.0):
        return 0
    sum/=len(arr)
    return sum


def printStats(arrV, arrN, arrS, noSummary, filestructure, filesuffix):

    if (len(arrV)==0):
        print " <i> - No result files detected</i>"
        return
    
    out=[[],[],[],[],[],[]]
    string=[]
    for c in xrange(len(arrV)):
        string+=["<div>%17s</div>" % arrN[c]]
        formatString="%17.2f"
        if (min(arrV[c]) > 0.0 and min(arrV[c])<0.009):
            formatString="%17.2e"
        out[0].append(formatString % (min(arrV[c])))
        out[1].append(formatString % (max(arrV[c])))
        out[2].append(formatString % (average(arrV[c])))
        out[3].append(formatString % (ste(arrV[c])))
        if (percent):
            out[4].append(formatString % (per(arrV[0],arrV[c])))
        out[5].append(formatString % (sum(arrV[c])))

    print "<table class='data'><thead><tr><th><div style='width:25px'><div></th><th>"+("</th><th>").join(string)+"</th><th class='left'>File</th></tr></thead>"
    if(printing and arrS!=0 ):       
        print "<tbody>"
        for l in arrS:
            resultPerS=[]
            for e in l[0]:
                formatString="%17.2f"
                if (e > 0.0 and e<0.009):
                    formatString="%17.2e"
                resultPerS+=[ formatString % e ]
            
            print "<tr><td></td><td>"+("</td><td>").join(resultPerS)+"</td><td class='left'>"+ removeSuffix(removePrefix(l[1], filestructure), filesuffix)+"</td></tr>"
        print "</tbody>"
            
    if(noSummary):
        print "</table>"
        return
        
    elif(arrS==0 or len(arrS)>1):

        print "<tfoot>"
        print "<tr><td class='left'>sum</td><td>"+("</td><td>").join(out[5])+"</td><td class='left'><i>aggregation</i></td></tr>"
        print "<tr><td class='left'>min</td><td>"+("</td><td>").join(out[0])+"</td><td></td></tr>"
        print "<tr><td class='left'>max</td><td>"+("</td><td>").join(out[1])+"</td><td></td></tr>"
        print "<tr><td class='left'>av</td><td>"+("</td><td>").join(out[2])+"</td><td></td></tr>"
        print "<tr><td class='left'>ste</td><td>"+("</td><td>").join(out[3])+"</td><td></td></tr>"
        if (percent):
            print "<tr><td>av%</td><td>"+"</td><td>".join(out[4])+"</td><td></td></tr>"
        print "</tfoot>"
    print "</table>"

# sam statiscis for initial aligment
def samstats(statsfile):
    names=["Total reads","QCfail","Duplicates","Duplicates %","Mapped","Mapped %","Paired", "Paired %", "Singletons"]
    values=[]
    st=re.split("[\n]+",open(statsfile).read())
    total= int(st[0].strip().split(" ")[0])
    QCfail = int(st[0].strip().split(" ")[2])
    dupl = int(st[1].strip().split(" ")[0]) 
    duplPercent = float(dupl)/(float(total)+pseudocount)*100
    mapped = int(st[2].strip().split(" ")[0])
    mappedPercent = float(mapped)/(float(total)+pseudocount)*100
    paired = int(st[6].strip().split(" ")[0])
    pairedPercent = float(paired)/(float(total)+pseudocount)*100
    singletons = int(st[8].strip().split(" ")[0])
    values = [total, QCfail, dupl, duplPercent, mapped, mappedPercent, paired, pairedPercent, singletons ]
    customRegion = open(statsfile).read().split("#custom region")
    if (len(customRegion) >= 2):
        names += ["Region mapped", "Region mapped %", "Region paired", "Region paired %"]
        st = customRegion[1].split("\n")
        regmapped = int(st[1].strip().split(" ")[0])
        regmappedPercent = float(regmapped)/(float(total)+pseudocount)*100
        regpaired = int(st[2].strip().split(" ")[0])
        regpairedPercent = float(regpaired)/(float(total)+pseudocount)*100
        values += [regmapped, regmappedPercent, regpaired, regpairedPercent]
#    if (len(st)>75):
#	print st
#        values.append(int(st[75])) # regmapped
#        values.append(float(values[-1])/values[0]*100) 
#        values.append(int(st[80]))
#        values.append(float(values[-1])/values[0]*100)
        
    return names,values


def fastqstats(statsfile):
	names,values=samstats(statsfile)
	return [names[0]],[values[0]]


# blue statiscis for initial aligment
def blue(statsfile):
    names=["Total reads","Reads OK","Reads OK%","healed","healed%","not healed", "not healed %","discarded", "discarded %","subs","dels","ins"]
    values=[]
    st=re.split("[\n]+",open(statsfile).read())
    values.append(int(st[2].split("\t")[0])) # total
    values.append(int(st[3].split("\t")[0])) # acepted
    values.append(float(values[-1])/float(values[0])*100) # %
    values.append(int(st[7].split("\t")[0])) # healed
    values.append(float(values[-1])/float(values[0])*100) # %
    values.append(int(st[11].split("\t")[0])) # not healed
    values.append(float(values[-1])/float(values[0])*100) # %
    values.append(int(st[17].split("\t")[0])) # discarded
    values.append(float(values[-1])/float(values[0])*100) # %
    values.append(int(st[20].split("\t")[0]))
    values.append(int(st[21].split("\t")[0]))
    values.append(int(st[22].split("\t")[0]))
        
    return names,values



# sam statiscis for initial aligment
def tophat(statsfile):
    names=["Total reads","Accepted","QCfail","Duplicates","Duplicates %","Mapped","Mapped %","Paired", "Paired %", "Singletons"]
    values=[]
    st=re.split("[ \n]+",open(statsfile).read())
    values.append(int(st[73])) # total
    values.append(int(st[0])) # acepted
    values.append(int(st[2])) # QCfail
    values.append(int(st[10])) # dupl
    values.append(float(values[-1])/float(values[1])*100) # dupl%
    values.append(int(st[14])) # mapped
    values.append(float(values[-1])/float(values[1])*100) # mapped %
    values.append(int(st[33])) # paired
    values.append(float(values[-1])/float(values[1])*100) # paired %
    values.append(int(st[47]))
    if (len(st)>76):
        names.append("Junction")
        names.append("Junction %")
        names.append("Jnct over GTF")
        names.append("Jnct over GTF %")
        values.append(int(st[76])) # junction reads
        #values.append(float(values[-1])/float(values[5])) # junction %
        values.append(float(values[-1])/float(values[1])*100) # junction %
        values.append(int(st[79])) # junction reads in ncbi genes
        values.append(float(values[-1])/float(values[10])*100) # junction reads in ncbi genes %
        
    return names,values

def cufflinksStats(logFile):
    
    names=["Transcripts", "Skipped", "Genes FPKM", "% OK", "% Lowdata", "Isoforms FPKM", "% OK", "% Lowdata"]
    values=[]
    file=open(logFile).read()
    # populate
    tmp=file.split("transcripts.gtf")[1].strip().split()[0]
    TS=float(tmp.strip())
    values.append(TS)

    tmp=file.split("skipped.gtf")[1].strip().split()[0]
    SK=float(tmp.strip())
    values.append(SK)

    tmp=file.split("genes.fpkm_tracking")[1].split(";")[0].strip().split()[0]
    GF=float(tmp.strip())
    values.append(GF)

    try:
        if not "OK;" in file.split("genes.fpkm_tracking")[1].split("\n")[0]: 
            raise error

        tmp=file.split("genes.fpkm_tracking")[1].split("OK;")[0].split(";")[-1].strip()
        GO=float(tmp.strip())
        values.append(100* GO / GF)
    except:
        values.append(0)
            

    try:
        if not "LOWDATA;" in file.split("genes.fpkm_tracking")[1].split("\n")[0]: 
            raise error
        tmp=file.split("genes.fpkm_tracking")[1].split("LOWDATA;")[0].split(";")[-1].strip()
        GL=float(tmp.strip())
        values.append(100* GL / GF)
    except:
        values.append(0)
    
    tmp=file.split("isoforms.fpkm_tracking")[1].split(";")[0].strip().split()[0]
    IF=float(tmp.strip())
    values.append(IF)

    try:    
        if not "OK;" in file.split("isoforms.fpkm_tracking")[1].split("\n")[0]: 
            raise error

        tmp=file.split("isoforms.fpkm_tracking")[1].split("OK;")[0].split(";")[-1].strip()
        IO=float(tmp.strip())
        values.append(100* IO / IF)
    except:
        values.append(0)
    
    try:
        if not "LOWDATA;" in file.split("isoforms.fpkm_tracking")[1].split("\n")[0]: 
            raise error

        tmp=file.split("isoforms.fpkm_tracking")[1].split("LOWDATA;")[0].split(";")[-1].strip()
        IL=float(tmp.strip())
        values.append(100 * IL / IF)
    except:
        values.append(0)
    
    
    return names, values


def htseqcountStats(logFile):
    
    names=[]
    values=[]
    lines=open(logFile).read().split("\n")
    for f in lines:
        cols = f.split(" ")
        if (len(cols)<6): 
            continue
        names+=[" ".join(cols[0:3])+" no_feature", " ".join(cols[0:3])+" ambiguous"]
        values+=[ float(cols[4]), float(cols[6]) ]
    return names, values


def onTarget(statsfile):
    names=["Total reads", "Total paired", "Total  Paired(%)" ,"OnTarget 100","(%)", "Paired on Target 100","(%)"]
    values=[]
    f=open(statsfile).read()
#    print f
    st=re.split("[ \n]+",f)
#    print st
#    print st[205]
    values.append(int(st[75])) # total
    values.append(int(st[85])) # paired total
    values.append(float(values[-1])/float(values[0])*100) # paired total %
    values.append(int(st[0])) # on target
    values.append(float(values[-1])/float(values[0])*100) # on taget %
    values.append(int(st[33])) # paired
    values.append(float(values[-1])/float(values[0])*100) # paired %
    if (f.find("# on target 0")>-1):
        names.append("pairedd oT 200")
        names.append("(%)")
        names.append("paired oT 0")
        names.append("(%)")
        values.append(int(st[205])) # on target 200
        values.append(float(values[-1])/float(values[0])*100) # on taget 200 %
        values.append(int(st[128])) # on target 0
        values.append(float(values[-1])/float(values[0])*100) # on taget 0 %
        
    return names,values


# picard RNAseq annotation
def annoStatsPicard(statsfile):
    names=["Total bases", "Aligned bases","(%)", "Ribosomal bases", "(%)" ,"Coding","(%)", "UTR","(%)","Intronic","(%)", "Intergenic", "(%)", "MEDIAN_CV_COVERAGE","MEDIAN_5PRIME_BIAS","MEDIAN_3PRIME_BIAS"]
    values=[]
    f=open(statsfile).read().split("\n")[7]
#    print f
    st=re.split("[ \t]+",f)
    for i in range(0,7):
        values.append(int(st[i]))
        if(i>0):
            values.append(float(values[-1]/float(values[0])*100))
    values.append(float(st[19]))
    values.append(float(st[20]))
    values.append(float(st[21]))
    return names,values

#def annoStats(statsfile):
#    names=["total", "genes", "(%)", "rRNA", "(%)", "tRNA", "(%)", "lincRNA", "(%)", "miRNA", "(%)", "snoRNA", "(%)", "snrna", "(%)", "miscRNA", "(%)", "polyA", "(%)", "other", "(%)", "HiSeq", "(%)", "ucsc_rRNA", "(%)", "SegDups", "(%)"]
#    values=[]
#    f=open(statsfile).read().split("\n")[1]
#    st=re.split("[ \t]+",f)
#    values.append(float(st[0]))
#    for i in range(1,14):
#        values.append(float(st[i]))
#        values.append(float(values[-1]/values[0]*100))
#    return names,values


def annoStats(statsfile):
	names=[]
	values=[]
	f=open(statsfile).read().split("\n")
	n=re.split("[ \t]+",f[0].strip())
	v=map(float,re.split("[ \t]+",f[2].strip())[1:])
	pairs=zip(n,v)
	pairs.sort()
	names = [ p[0] for p in pairs ]
	values = [ p[1] for p in pairs ]
	return names,values

def parsetime(string):
    arr=string.split(":")
    return (int(arr[0])*3600+int(arr[1])*60+int(arr[2]))/3600.0

def time(file):
    st=re.split("[ =,k\n]",open(file).read().split("cput")[1])
    real=parsetime(st[1]) #datetime.datetime.strptime(st[1], "%I:%M:%S")
    computer=parsetime(st[9])
    memory=int(st[6])*9.5367e-7
#    end=datetime.datetime.strptime(st.split(">>>>> enddate ")[1].split("\n")[0], "%a %b %d  %H:%M:%S EST %Y")
#    delta=end-start
    return ["cpu(h)", "real", "memory(GB)"],[real,computer,memory]


# sam statistics for recalibrated scores and realigned bams
def samstatsrecal(statsfile):
    # properties of the recalibrated file
    names,values = samstats (statsfile)

    # properties of the initial files
    n,valun = samstats(statsfile.replace(".recal.clean", ""))

    names+=["imprpaired"]
    values.append(valun[5]-values[5])
    if (len(valun)>9 and len(values)>9):
        names+=["imprregpaired"]
        values.append(valun[9]-values[9])
    return names,values
    

# sam statiscis for initial aligment
def bamDist(bamfile, col):

    names=[]
    values=[]
    pre=""
    if(col==1):
        pre="mapped"
    elif(col==5):
        pre="pqal20"

    st=os.popen3("python ../tools/bamchromdist.py "+bamfile)[1].read().split("\n")
    
    for s in st[2:-3]:
        arr=re.split("[ \t]+",s)
        if(arr[0]==""):
            arr.pop(0)
        names.append(pre+"_"+arr[0])
        values.append(int(arr[col]))
    if(printing):
        print str(values)+" "+bamfile
    return names,values


def variant(variantFile):
	names=["Total","known", "SNPdb Conc", "variantRatePerBp", "hetHomRatio", "novel", "variantRatePerBp","hetHomRatio"]
	values=[]
	file=open(variantFile).read()
#	CO=re.split("[ \t\n]+", file.split("CompOverlap")[3])
#	values.append(int(CO[5])) #total
	CO=re.split("[ \t\n]+", file.split("CountVariants")[3])
	values.append(int(CO[6])) #total
	CO=re.split("[ \t\n]+", file.split("CompOverlap")[4])
	values.append(int(CO[5])) #known
	values.append(float(CO[10])) #SNPconc

	CV=re.split("[ \t\n]+", file.split("CountVariants")[4])
	values.append(float(CV[10])) #variantRatePerBp
	values.append(float(CV[26])) #hetHomRatio
	
	CO=re.split("[ \t\n]+", file.split("CompOverlap")[5])
	values.append(int(CO[5])) #novel
#	values.append(float(CO[10])) #SNPconc

	CV=re.split("[ \t\n]+", file.split("CountVariants")[5])
	values.append(float(CV[10])) #variantRatePerBp
	values.append(float(CV[26])) #hetHomRatio

	return names,values


def trimgaloreStats(logFile):
    names=["Processed reads", "Trimmed reads","%","Too short reads","%","Too long reads","%", "Remaining reads","%"]
    values=[]
    file=open(logFile).read()
    # populate
    tmp=file.split("Processed reads:")[1].strip().split()[0]
    PR=float(tmp.strip())
    values.append(PR)
   
    tmp=file.split("Trimmed reads:")[1].split("(")[0]
    TR=float(tmp.strip())
    values.append(TR)
    values.append(100*TR/PR)

    tmp=file.split("Too short reads:")[1].split("(")[0]
    TS=float(tmp.strip())
    values.append(TS)
    values.append(100*TS/PR)

    tmp=file.split("Too long reads:")[1].split("(")[0]
    TL=float(tmp.strip())
    values.append(TL)
    values.append(100*TL/PR)

    tmp=file.split("remaining reads ")[1].split()[0]
    RM=float(tmp.strip())
    values.append(RM)
    values.append(100*RM/PR)

    return names, values

def trimmomaticStats(logFile):
    values=[]
    if "Input Read Pairs" in open(logFile).read():
        # paired library
        names=["Input Read Pairs", "Both Surviving","%","Forward Only","%","Reverse Only","%", "Dropped","%"]
        file=open(logFile).read()

        # populate
        tmp=file.split("Input Read Pairs:")[1].split("Both Surviving")[0].strip()
        PR=float(tmp.strip())
        values.append(PR)

        tmp=file.split("Both Surviving:")[1].split("(")[0]
        BSR=float(tmp.strip())
        values.append(BSR)
        values.append(100*BSR/PR)

        tmp=file.split("Forward Only Surviving:")[1].split("(")[0]
        FS=float(tmp.strip())
        values.append(FS)
        values.append(100*FS/PR)

        tmp=file.split("Reverse Only Surviving:")[1].split("(")[0]
        RS=float(tmp.strip())
        values.append(RS)
        values.append(100*RS/PR)

        tmp=file.split("Dropped:")[1].split("(")[0]
        DR=float(tmp.strip())
        values.append(DR)
        values.append(100*DR/PR)

    else:
        # single library
        names=["Input Reads", "Surviving","%","Forward Only","%","Reverse Only","%","Dropped","%"]
        file=open(logFile).read()

        # populate
        tmp=file.split("Input Reads:")[1].split("Surviving")[0].strip()
        PR=float(tmp.strip())
        values.append(PR)

        tmp=file.split("Surviving:")[1].split("(")[0]
        SR=float(tmp.strip())
        values.append(SR)
        values.append(100*SR/PR)
     
        values.append(0.)
        values.append(0.)
        values.append(0.)
        values.append(0.)

        tmp=file.split("Dropped:")[1].split("(")[0]
        DR=float(tmp.strip())
        values.append(DR)
        values.append(100*DR/PR)

    return names, values

def cutadaptStats(logFile):
    names=["Processed reads", "Trimmed reads","%","Too short reads","%","Too long reads","%", "Remaining reads","%"]
    values=[]
    file=open(logFile).read()
    # populate
    tmp=file.split("Processed reads:")[1].strip().split()[0]
    PR=float(tmp.strip())
    values.append(PR)

    tmp=file.split("Trimmed reads:")[1].split("(")[0]
    TR=float(tmp.strip())
    values.append(TR)
    values.append(100*TR/PR)

    tmp=file.split("Too short reads:")[1].split("(")[0]
    TS=float(tmp.strip())
    values.append(TS)
    values.append(100*TS/PR)

    tmp=file.split("Too long reads:")[1].split("(")[0]
    TL=float(tmp.strip())
    values.append(TL)
    values.append(100*TL/PR)

    tmp=file.split("remaining reads ")[1].split()[0]
    RM=float(tmp.strip())
    values.append(RM)
    values.append(100*RM/PR)

    return names, values

def hiclibStats(logFile):
    names=["Original reads", "Unmapped", "(%)", "Semi-dangling", "(%)", "Duplicates", "(%)", "Sm/Lg Fragments","(%)", "Extreme Fragments","(%)", "# Reads (final)", "(%)", "# Fragments (final)"]
    values=[]
    file=open(logFile).read()
    # populate
    tmp=file.split("Original reads:")[1].split("Number of reads changed")[1].split("Fragments number changed")[0].split("--->")
    OR=float(tmp[0].strip())
    values.append(OR)
    UM=OR-float(tmp[1].strip())
    values.append(UM)
    values.append(100*UM/OR)

    tmp=file.split("Semi-dangling end filter")[1].split("Number of reads changed")[1].split("Fragments number changed")[0].split("--->")
    SD=float(tmp[0].strip())-float(tmp[1].strip())
    values.append(SD)
    values.append(100*SD/OR)

    tmp=file.split("Filtering duplicates in DS reads")[1].split("Number of reads changed")[1].split("Fragments number changed")[0].split("--->")
    DR=float(tmp[0].strip())-float(tmp[1].strip())
    values.append(DR)
    values.append(100*DR/OR)

    tmp=file.split("Small/large fragments filter")[1].split("Number of reads changed")[1].split("Fragments number changed")[0].split("--->")
    SLF=float(tmp[0].strip())-float(tmp[1].strip())
    values.append(SLF)
    values.append(100*SLF/OR)

    tmp=file.split("Extreme fragments filter")[1].split("Number of reads changed")[1].split("Fragments number changed")[0].split("--->")
    EF=float(tmp[0].strip())-float(tmp[1].strip())
    values.append(EF)
    values.append(100*EF/OR)

    RF=float(tmp[1].strip())
    values.append(RF)
    values.append(100*RF/OR)

    tmp=file.split("Extreme fragments filter")[1].split("Fragments number changed")[1].split("--->")[1].split()
    FF=float(tmp[0].strip())
    values.append(FF)

    return names,values

def hicupStats(statsFile):
    names=[]
    values=[]
    file=open(statsFile).read()

    if (file.find("Average_length_truncated_sequence")>-1):
    	lines = file.split("\n")
    	p1 = lines[1].split("\t")
    	p2 = lines[2].split("\t")    	
    	names += ["P1 non-trunc reads", "%","P1 trunc reads", "%", "P2 non-trunc reads","%", "P2 trunc reads", "%"]
    	values += [float(p1[4]), float(p1[5]), float(p1[2]), float(p1[3]), float(p2[4]), float(p2[5]), float(p2[2]), float(p2[3])]
  
    if (file.find("Total_reads_processed")>-1):
    	lines = file.split("\n")
    	names += ["P1 reads", "P2 reads", "Unique Align P1", "%", "Unique Align P2", "%", "Multi mapper P1", "%", "Multi mapper P2", "%", "Nonaligned P1", "%", "Nonaligned P2", "%", "Paired P1", "%", "Paired P2", "%" ]
    	p1 = lines[1].split("\t")
    	p2 = lines[2].split("\t")
    	values += [float(p1[1]), float(p2[1]), float(p1[4]), float(p1[5]), float(p2[4]), float(p2[5]), float(p1[6]), float(p1[7]), float(p2[6]), float(p2[7]), float(p1[8]), float(p1[9]), float(p2[8]), float(p2[9]), float(p1[10]), float(p1[11]), float(p2[10]), float(p2[11])] 	

    if (file.find("Same_circularised")>-1):
    	lines = file.split("\n")
  	
    	sig = re.sub('[_]', ' ', lines[0]).split("\t")[1:]
    	val = lines[1].split("\t")[1:]
    	names += sig
    	values += [float(i) for i in val]

    if (file.find("Read_pairs_processed")>-1):
    	lines = file.split("\n")
    	sig = re.sub('[_]', ' ', lines[0]).split("\t")[1:]
    	val = lines[1].split("\t")[1:]
    	names += sig
    	values += [float(i) for i in val]

    return names,values

def fithicStats(statsFile):
    names=["Intra In Range Count", "Intra Out Of Range Count", "Intra Very Proximal Count","Inter Count" ]
    values=[]
    file=open(statsFile).read()

    # populate
    tmp=file.split("intraInRangeCount ")[1].split()[0]
    IR=float(tmp.strip())
    values.append(IR)

    tmp=file.split("intraOutOfRangeCount ")[1].split()[0]
    RC=float(tmp.strip())
    values.append(RC)

    tmp=file.split("intraVeryProximalCount ")[1].split()[0]
    PC=float(tmp.strip())
    values.append(PC)

    tmp=file.split("interCount ")[1].split()[0]
    IC=float(tmp.strip())
    values.append(IC)

    return names,values

def homerchipseqStats(logFile):
    names=["Est. fragment length","Autocorr. same strand", "Autocorr. diff. strand", "Autocorr. same/diff", "# Peaks","Peak size","IP efficiency (%)"]
    values=[]
    file=open(logFile).read()
    # populate
    tmp=file.split("# fragment length =")[1].strip().split()[0]
    FL=float(tmp.strip())
    values.append(FL)

    tmp=file.split("Same strand fold enrichment:")[1].strip().split()[0]
    AS=float(tmp.strip())
    values.append(AS)
    
    tmp=file.split("Diff strand fold enrichment:")[1].strip().split()[0]
    AD=float(tmp.strip())
    values.append(AD)

    tmp=file.split("Same / Diff fold enrichment:")[1].strip().split()[0]
    AF=float(tmp.strip())
    values.append(AF)
    
    tmp=file.split("# total peaks =")[1].strip().split()[0]
    TP=float(tmp.strip())
    values.append(TP)

    tmp=file.split("# peak size =")[1].strip().split()[0]
    PS=float(tmp.strip())
    values.append(PS)

    tmp=file.split("# Approximate IP efficiency =")[1].strip().split("%")[0]
    IP=float(tmp.strip())
    values.append(IP)

    return names,values

def peakrangerStats(logFile):
    names=["Estimated noise rate", "Total reads", "Unique reads","Library complexity","Peaks","Summits"]
    values=[]
    file=open(logFile).read()
    # populate
    tmp=file.split("Estimated noise rate:")[1].strip().split()[0]
    EN=float(tmp.strip())
    values.append(EN)

    tmp=file.split("Total reads:")[1].strip().split()[0]
    TR=float(tmp.strip())
    values.append(TR)

    tmp=file.split("Unique reads:")[1].strip().split()[0]
    UR=float(tmp.strip())
    values.append(UR)

    tmp=file.split("Library complexity:")[1].strip().split("%")[0]
    LC=float(tmp.strip())
    values.append(LC)

    tmp=file.split("Peaks:")[1].strip().split()[0]
    PE=float(tmp.strip())
    values.append(PE)

    return names, values

def macs2Stats(logFile):
    names=["Total IP tags", "IP (filtered)", "%","Control tags","Control (filtered)","%","Paired peaks","Fragment length","Refined peaks"]
    values=[]
    file=open(logFile).read()
    # populate
    tmp=file.split("#1  total tags in treatment:")[1].strip().split()[0]
    TT=float(tmp.strip())
    values.append(TT)

    tmp=file.split("#1  tags after filtering in treatment:")[1].strip().split()[0]
    TF=float(tmp.strip())
    values.append(TF)

    tmp=file.split("#1  Redundant rate of treatment:")[1].strip().split()[0]
    TR=float(tmp.strip())
    values.append(100.-TR)
    
    try:
        tmp=file.split("#1  total tags in control:")[1].strip().split()[0]
        CT=float(tmp.strip())
        values.append(CT)

        tmp=file.split("#1  tags after filtering in control:")[1].strip().split()[0]
        CF=float(tmp.strip())
        values.append(CF)

        tmp=file.split("#1  Redundant rate of control:")[1].strip().split()[0]
        CR=float(tmp.strip())
        values.append(100.-CR)
    except:
        values.append(0)
        values.append(0)
        values.append(0)

    try:
        tmp=file.split("#2 number of paired peaks:")[1].strip().split()[0]
        PP=float(tmp.strip())
        values.append(PP)
    except:
        values.append(0)

    tryalternative=False
    try:
        tmp=file.split("#2 predicted fragment length is")[1].strip().split("bps")[0]
        PF=float(tmp.strip())
        values.append(PF)
        
    except:
        tryalternative=True

    if (tryalternative):
        try:
            tmp=file.split(" as fragment length")[0].strip().split("#2 Use ")[1]
            PF=float(tmp.strip())
            values.append(PF)
    
        except:
            values.append(0)
     
    tmp=file.split("Final number of refined peaks:")[1].strip().split()[0]
    RP=float(tmp.strip())
    values.append(RP)
           
    return names, values

def chanceStats(logFile):
    names=[ "Cumulative % enrichment", "Input scaling factor", "Diff. Enrichment", "FDR (overall)", "FDR (TF normal)", "FDR (Histone normal)", "FDR (TF cancer)", "FDR (Histone cancer)"]
    values=[]
    file=open(logFile).read()
    # populate
    tmp=file.split("percent_genome_enriched,")[1].strip().split()[0]
    PE=float(tmp.strip())
    values.append(PE)

    tmp=file.split("input_scaling_factor,")[1].strip().split()[0]
    IS=float(tmp.strip())
    values.append(IS)

    tmp=file.split("differential_percentage_enrichment,")[1].strip().split()[0]
    DE=float(tmp.strip())
    values.append(DE)

    tmp=file.split("fdr,")[1].strip().split()[0]
    FD=float(tmp.strip())
    values.append(FD)
    
    tmp=file.split("tfbs_normal_fdr,")[1].strip().split()[0]
    TN=float(tmp.strip())
    values.append(TN)
      
    tmp=file.split("histone_normal_fdr,")[1].strip().split()[0]
    HN=float(tmp.strip())
    values.append(HN)
    
    tmp=file.split("tfbs_cancer_fdr,")[1].strip().split()[0]
    TC=float(tmp.strip())
    values.append(TC)
    
    tmp=file.split("histone_cancer_fdr,")[1].strip().split()[0]
    HC=float(tmp.strip())
    values.append(HC)
    
    return names, values


def fastqscreenStats(logFile):
    file=open(logFile).read()
    names=[]
    species=file.split("\n")[2:-3]
    for p in [ s.split()[0] for s in species ]:
        names+=[ p+" (unique hit)" ]
        names+=[ p+" (multi map)" ]
    
    names+=["Hit no libraries"]
    values=[]

    # populate
    for s in species:
        values.append(float(s.split()[2]))
        values.append(float(s.split()[3]))

    tmp=file.split("%Hit_no_libraries:")[1].strip().split()[0]
    
    TT=float(tmp.strip())
    values.append(TT)
    
    return names, values

def memechipStats(logFile):
    names=["Peak regions", "with strong sites","%", "w/o strong sites","%"]
    values=[]
    file=open(logFile).read()
    # populate
    tmp=file.split("Peak regions:")[1].strip().split()[0]
    PR=float(tmp.strip())
    values.append(PR)

    tmp=file.split("bound directely (strong site):")[1].strip().split()[0]
    SS=float(tmp.strip())
    values.append(SS)
    values.append(100. * SS/PR)

    tmp=file.split("bound indirectely (weak or no site):")[1].strip().split()[0]
    WS=float(tmp.strip())
    values.append(WS)
    values.append(100. * WS/PR)

    return names, values
 
 
def inchwormStats(logFile):
    names=[ "Read pairs", "Jellyfish kmers", "Inchworm fasta"]
    values=[]
    file=open(logFile).read()
    # populate

    tmp=file.split("both fasta:")[1].strip().split()[0]
    RP=float(tmp.strip())
    values.append(RP)
        
    tmp=file.split("jellyfish kmers:")[1].strip().split()[0]
    JK=float(tmp.strip())
    values.append(JK)

    tmp=file.split("inchworm fasta:")[1].strip().split()[0]
    IF=float(tmp.strip())
    values.append(IF)
    return names, values
    
def chrysalisStats(logFile):
    names=["Read count"]
    values=[]
    file=open(logFile).read()
    # populate
    tmp=file.split("Read count:")[1].strip().split()[0]
    RC=float(tmp.strip())
    values.append(RC)
    
    return names, values
    
def butterflyStats(logFile):
    names=[ "Transcripts"]
    values=[]
    file=open(logFile).read()
    # populate
    
    tmp=file.split("Transcripts:")[1].strip().split()[0]
    RP=float(tmp.strip())
    values.append(RP)
    
    return names, values

def bigwigStats(logFile):
    names=["Library Size", "Normalized to", "Scale factor"]
    values=[]
    file=open(logFile).read()
    # populate
    tmp=file.split("library size:")[1].strip().split()[0]
    LS=float(tmp.strip())
    values.append(LS)

    tmp=file.split("normalize to:")[1].strip().split()[0]
    NT=float(tmp.strip())
    values.append(NT)

    tmp=file.split("scale factor:")[1].strip().split()[0]
    SF=float(tmp.strip())
    values.append(SF)

    return names, values

#################33
# TEMP

def intersection(variantFile):
    #names=["Intersection", "Post", "Het/Hom", "Ti/Tv","known", "novel", "Pre", "Het/Hom", "Ti/Tv","known", "novel" ]
    names=["Total","Het/Hom","nHets","nHomRef","nHomVar","overlHapMap77", "concHapMap77", "overlHapMap78","concHapMap78", "overl1000G78", "conc1000G78", "overldbSNP", "concdbSNP", "Intersection","Het/Hom","nHets","nHomRef","nHomVar","overlHapMap77","concHapMap77","overlHapMap78","concHapMap78", "overl1000G78","conc1000G78","overldbSNP", "concdbSNP","Post","Het/Hom","nHets","nHomRef","nHomVar","overlHapMap77", "concHapMap77","overlHapMap78","concHapMap78","overl1000G78", "conc1000G78", "overldbSNP","concdbSNP","Pre", "Het/Hom","nHets","nHomRef","nHomVar","overlHapMap77","concHapMap77","overlHapMap78","concHapMap78","overl1000G78", "conc1000G78","overldbSNP","concdbSNP"]
    values=[]
    file=open(variantFile).read()
    CO=file.split("CompOverlap :")[1].split("\n")
    CV=file.split("CountVariants :")[1].split("\n")
    #TI=file.split("TiTvVariantEvaluator :")[1].split("\n")
    #print re.split("[ \t]+",CV[1])
    #print re.split("[ \t]+",CV[1])[25]
    #print re.split("[ \t]+",CV[3])[25]
    values.append(int(re.split("[ \t]+",CO[17])[5])) #total 
    values.append(float(re.split("[ \t]+",CV[17])[25])) #het/hom
    values.append(float(re.split("[ \t]+",CV[17])[18])) #nHets
    values.append(float(re.split("[ \t]+",CV[17])[19])) #nHomRef
    values.append(float(re.split("[ \t]+",CV[17])[20])) #nHomVar
    values.append(float(re.split("[ \t]+",CO[17+18])[9])) #1hapmap1
    values.append(float(re.split("[ \t]+",CO[17+18])[10])) #1hapmap1
    values.append(float(re.split("[ \t]+",CO[17+(18*2)])[9])) #1hapmap1
    values.append(float(re.split("[ \t]+",CO[17+(18*2)])[10])) #1hapmap1
    values.append(float(re.split("[ \t]+",CO[17])[9])) #1kg
    values.append(float(re.split("[ \t]+",CO[17])[10])) #1kg
    values.append(float(re.split("[ \t]+",CO[17+(18*3)])[9])) #snpdb
    values.append(float(re.split("[ \t]+",CO[17+(18*3)])[10])) #snpdb

    values.append(int(re.split("[ \t]+",CO[8])[5])) #Intersection
    values.append(float(re.split("[ \t]+",CV[8])[25])) #het/hom
    values.append(float(re.split("[ \t]+",CV[8])[18])) #nHets
    values.append(float(re.split("[ \t]+",CV[8])[19])) #nHomRef
    values.append(float(re.split("[ \t]+",CV[8])[20])) #nHomVar
    values.append(float(re.split("[ \t]+",CO[8+18])[9])) #1hapmap1
    values.append(float(re.split("[ \t]+",CO[8+18])[10])) #1hapmap1
    values.append(float(re.split("[ \t]+",CO[8+(18*2)])[9])) #1hapmap1
    values.append(float(re.split("[ \t]+",CO[8+(18*2)])[10])) #1hapmap1
    values.append(float(re.split("[ \t]+",CO[8])[9])) #1kg
    values.append(float(re.split("[ \t]+",CO[8])[10])) #1kg
    values.append(float(re.split("[ \t]+",CO[8+(18*3)])[9])) #snpdb  
    values.append(float(re.split("[ \t]+",CO[8+(18*3)])[10])) #snpdb

    values.append(int(re.split("[ \t]+",CO[2])[5])) #Post+filtered inPre
    values[-1]+=(int(re.split("[ \t]+",CO[11])[5])) #Post
    values.append(float(re.split("[ \t]+",CV[2])[25])) #Post+filtered inPre het/hom
    values[-1]+=(float(re.split("[ \t]+",CV[11])[25])) #Post het/hom
    values[-1]/=2
    values.append(float(re.split("[ \t]+",CV[2])[18])) #Post+fiPre nHets
    values[-1]+=(float(re.split("[ \t]+",CV[11])[18])) #Post nHets
    values.append(float(re.split("[ \t]+",CV[2])[19])) #Post+fiPre nHomRef
    values[-1]+=(float(re.split("[ \t]+",CV[11])[19])) #Post nHomRef
    values.append(float(re.split("[ \t]+",CV[2])[20])) #Post+fiPre nHomVar
    values[-1]+=(float(re.split("[ \t]+",CV[11])[20])) #Post nHomVar
    values.append(float(re.split("[ \t]+",CO[2+18])[9])) #1hapmap1
    values[-1]+=(float(re.split("[ \t]+",CO[11+18])[9])) #1hapmap1
    values.append(float(re.split("[ \t]+",CO[2+18])[10])) #1hapmap1
    values[-1]+=(float(re.split("[ \t]+",CO[11+18])[10])) #1hapmap1
    values[-1]/=2
    values.append(float(re.split("[ \t]+",CO[2+(18*2)])[9])) #1hapmap2
    values[-1]+=(float(re.split("[ \t]+",CO[11+(18*2)])[9])) #1hapmap2    
    values.append(float(re.split("[ \t]+",CO[2+(18*2)])[10])) #1hapmap2
    values[-1]+=(float(re.split("[ \t]+",CO[11+(18*2)])[10])) #1hapmap2
    values[-1]/=2
    values.append(float(re.split("[ \t]+",CO[2])[9])) #1kg
    values[-1]+=(float(re.split("[ \t]+",CO[11])[9])) #1kg
    values.append(float(re.split("[ \t]+",CO[2])[10])) #1kg
    values[-1]+=(float(re.split("[ \t]+",CO[11])[10])) #1kg
    values[-1]/=2
    values.append(float(re.split("[ \t]+",CO[2+(18*3)])[9])) #snpdb
    values[-1]+=(float(re.split("[ \t]+",CO[11+(18*3)])[9])) #snpdb
    values.append(float(re.split("[ \t]+",CO[2+(18*3)])[10])) #snpdb
    values[-1]+=(float(re.split("[ \t]+",CO[11+(18*3)])[10])) #snpdb
    values[-1]/=2

    values.append(int(re.split("[ \t]+",CO[5])[5])) #Pre+filtered inPost
    values[-1]+=(int(re.split("[ \t]+",CO[14])[5])) #Pre
    values.append(float(re.split("[ \t]+",CV[5])[25])) #Post+filtered inPre het/hom
    values[-1]+=(float(re.split("[ \t]+",CV[14])[25])) #Post het/hom
    values[-1]/=2
    values.append(float(re.split("[ \t]+",CV[5])[18])) #Post+fiPre nHets
    values[-1]+=(float(re.split("[ \t]+",CV[14])[18])) #Post nHets
    values.append(float(re.split("[ \t]+",CV[5])[19])) #Post+fiPre nHomRef
    values[-1]+=(float(re.split("[ \t]+",CV[14])[19])) #Post nHomRef
    values.append(float(re.split("[ \t]+",CV[5])[20])) #Post+fiPre nHomVar
    values[-1]+=(float(re.split("[ \t]+",CV[14])[20])) #Post nHomVar
    values.append(float(re.split("[ \t]+",CO[5+18])[9])) #1hapmap1
    values[-1]+=(float(re.split("[ \t]+",CO[14+18])[9])) #1hapmap1
    values.append(float(re.split("[ \t]+",CO[5+18])[10])) #1hapmap1
    values[-1]+=(float(re.split("[ \t]+",CO[14+18])[10])) #1hapmap1
    values[-1]/=2
    values.append(float(re.split("[ \t]+",CO[5+(18*2)])[9])) #1hapmap2
    values[-1]+=(float(re.split("[ \t]+",CO[14+(18*2)])[9])) #1hapmap2
    values.append(float(re.split("[ \t]+",CO[5+(18*2)])[10])) #1hapmap2
    values[-1]+=(float(re.split("[ \t]+",CO[14+(18*2)])[10])) #1hapmap2
    values[-1]/=2
    values.append(float(re.split("[ \t]+",CO[5])[9])) #1kg
    values[-1]+=(float(re.split("[ \t]+",CO[14])[9])) #
    values.append(float(re.split("[ \t]+",CO[5])[10])) #
    values[-1]+=(float(re.split("[ \t]+",CO[14])[10])) #
    values[-1]/=2
    values.append(float(re.split("[ \t]+",CO[5+(18*3)])[9])) #snpdb
    values[-1]+=(float(re.split("[ \t]+",CO[14+(18*3)])[9])) #
    values.append(float(re.split("[ \t]+",CO[5+(18*3)])[10])) #
    values[-1]+=(float(re.split("[ \t]+",CO[14+(18*3)])[10])) #
    values[-1]/=2


    return names,values


##########33
# Temp ende
##########

def coverage(file):
    cont=open(file).read().split("\n")
    names=cont[0].split("\t")[1::]
    values=[]
    for i in cont[1].split("\t")[1::]:
        values.append(float(i))
    #if(printing):
    #    print str(values)+" "+file
    return names,values                    



# adds an array [a3,b3,c3] as a column
# a1, a2, a3
# b1, b2, b3
# c1, c2, c3
def addValues(results,values):
    if(results==[]):
        for v in values:
            results.append([v])
    else:
        for v in range(0,len(values)):
            results[v].append(values[v])
    return results
        

#######3
# MAIN
#######3
oaresult=[]
for d in dir:
    result=[]
    psresult=[]
    name=glob.glob(d+'*'+ext)
    name.sort()
    for f in name:
        try:
            if (type=="samstats"):
                names,values=samstats(f)
            if (type=="fastqstats"):
                names,values=fastqstats(f)
            if (type=="samstatsrecal"):
                names,values=samstatsrecal(f)
            if (type=="bamdistMapped"):
                names,values=bamDist(f, 5)
            if (type=="coverage"):
                names,values=coverage(f)                    
            if (type=="variant"):
                names,values=variant(f)
            if (type=="tophat"):
                names,values=tophat(f)
            if (type=="cufflinks"):
                names,values=cufflinksStats(f)
            if (type=="htseqcount"):
                names,values=htseqcountStats(f)
            if (type=="times"):
                names,values=time(f)
            if (type=="target"):
                names,values=onTarget(f)
            if (type=="intersection"):
                names,values=intersection(f)     
            if (type=="annostats"):
                names,values=annoStats(f)
            if (type=="trimgalore"):
                names,values=trimgaloreStats(f)
            if (type=="trimmomatic"):
                names,values=trimmomaticStats(f)
            if (type=="cutadapt"):
                names,values=cutadaptStats(f)
            if (type=="hiclibMapping"):
                names,values=hiclibStats(f)
            if (type=="hicup"):
                names,values=hicupStats(f)
            if (type=="fithic"):
                names,values=fithicStats(f)
            if (type=="homerchipseq"):
                names,values=homerchipseqStats(f)
            if (type=="peakranger"):
                names,values=peakrangerStats(f)
            if (type=="macs2"):
                names,values=macs2Stats(f)
            if (type=="chance"):
                names,values=chanceStats(f)
            if (type=="memechip"):
                names,values=memechipStats(f)
            if (type=="fastqscreen"):
                names,values=fastqscreenStats(f)
            if (type=="blue"):
				names,values=blue(f)
            if (type=="trinity_inchworm"):
                names,values=inchwormStats(f)
            if (type=="trinity_chrysalis"):
                names,values=chrysalisStats(f)
            if (type=="trinity_butterfly"):
                names,values=butterflyStats(f)
            if (type=="bigwig"):
				names,values=bigwigStats(f)
            
            result=addValues(result,values)

            # only list file structure from current root
            filename="/".join(f.split("/")[-4::])
            if (link!=""):
                filename="<a href=\"%s\">%s</a>" % (link+"/"+filename, filename)
            psresult.append([values,filename])
            oaresult=addValues(oaresult,values)
                
        except:
            sys.stderr.write("error with "+f+"\n")
            traceback.print_exc()
            #sys.exit()

    filestructure="/".join(d.split("/")[-4::]) # only list file structure from current root
    print "<h3>"+ filestructure +"</h3>" 
    printStats(result,names,psresult,noSummary,filestructure, ext)

if (not noOverallSummary and overAll):
    print "<h3 class='overall'>Aggregation over all libraries</h3>"
    printStats(oaresult,names,0,noOverallSummary,"","")
