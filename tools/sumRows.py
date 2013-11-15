import sys, os, re, math, optparse


parser = optparse.OptionParser()

parser.add_option("-i","--input", dest = "fn", help = "input file")
parser.add_option("-s","--start", type="int",  dest = "start", help = "start collecting")
parser.add_option("-e","--end", type="int", dest = "end", help = "end collecting")
parser.add_option("-f","--header", action='store_true', dest="header", help = "header", default=False)
parser.add_option("-o","--oneline", action='store_true', dest="oneline", help = "oneline", default=False)
parser.add_option("-n","--names", dest = "names", help = "names", default="")
parser.add_option("-l","--location", type="int", dest = "location", help = "use location column's value", default=-1)
parser.add_option("--normalize", type="string", dest = "normalize", help = "normalize")
parser.add_option("--normalizeBy", type="float", dest = "normalizeBy", help = "normalize")
(options, args) = parser.parse_args()

names=options.names.split(",")
fn=options.fn
start=options.start
end=options.end
location=options.location
normalize=options.normalize
normalizeBy=options.normalizeBy

percent=False
printing=True

def nonZero(arr):
    sum=0
    for a in arr:
        if float(a)>0:
            sum+=1
    return sum


def average(arr):
    sum=0
    for a in arr:
        sum+=float(a)
    sum/=len(arr)
    return(sum)

def std(arr,av):
    sum=0
    for a in arr:
        sum+=(float(a)-float(av))*(float(a)-float(av))
    sum=math.sqrt(sum/(len(arr)))
    return(sum)

def ste(arr):
    av=average(arr)
    return(std(arr,av)/math.sqrt(len(arr)))

def per(max,arr):
    sum=0
    for a in range(0,len(arr)):
        if (float(arr[a])!=0.0):
            sum+=float(arr[a])/float(max[a])*100
    if (sum==0.0):
        return 0
    sum/=len(arr)
    return sum

def printStats(arrV, arrN, arrS):
    out=[[],[],[],[],[],[],[],[]]
    string="    "
    #normaliziation
    normsum=0
    normnZe=0

    for c in range(0,len(arrV)):
        string+="%17s " % arrN[c]
        formatString="%17.2f"
        if(min(arrV[c])<0.009):
            formatString="%17.2e"
            #formatString="%17.2f"     

        # normalize over all elements in the array (e.g. genes)
        normsum+=sum(arrV[c])
        normnZe+=nonZero(arrV[c])
                                                                                                                    
        out[0].append(formatString % (min(arrV[c])))
        out[1].append(formatString % (max(arrV[c])))
        out[2].append(formatString % (average(arrV[c])))
        out[3].append(formatString % (ste(arrV[c])))
        if (percent):
            out[4].append(formatString % (per(arrV[0],arrV[c])))
        out[5].append(formatString % (sum(arrV[c])))
        out[6].append(formatString % (nonZero(arrV[c])))

    #normalize
    if normalize:
        for i in range(0,len(out[5])):
            out[5][i]=str(float(out[5][i])/(normsum+1e-6))
        for i in range(0,len(out[6])):
            out[6][i]=str(float(out[6][i])/(normnZe+1e-6))
    if normalizeBy > 0:
        for i in range(0,len(out[5])):
            out[5][i]=str(float(out[5][i])/(normalizeBy+1e-6))
        for i in range(0,len(out[6])):
            out[6][i]=str(float(out[6][i])/(normalizeBy+1e-6))

    if(options.oneline):
		print "nZe %s sum %s min %s max %s av %s stde %s" % (" ".join(out[6])," ".join(out[5])," ".join(out[0])," ".join(out[1])," ".join(out[2])," ".join(out[3]))
		sys.exit(0)

    if(printing and arrS!=0 ):
        print string
        for l in arrS:
            resultPerS="    "
            for e in l[0]:
                formatString="%17.2f "
                if(e<0.009):
                    formatString="%17.2e "
                resultPerS+= formatString % e
            resultPerS+=" "+l[1]
            print resultPerS
        if(arrS==0 or len(arrS)>1):
            print "-----------------------------"
    if(arrS==0 or len(arrS)>1):
        print string
        print "nZe "+" ".join(out[6])
        print "sum "+" ".join(out[5])
        print "min "+" ".join(out[0])
        print "max "+" ".join(out[1])
        print " av "+" ".join(out[2])
        print "ste "+" ".join(out[3])
        if (percent):
            print "av% "+" ".join(out[4])


sumarr=[]
if names==[""]:
    name=[]
else:
    name=names
for i in range(start,end):
    sumarr+=[[0]]
    if (names==[""]):
        name.append("row"+str(i))

datastream=""
if (fn == "-"):
    datastream=sys.stdin
else:
	datastream=open(fn)

for i in datastream:
    if (i=="" or i[0]=="#" or options.header):
#        print >> sys.stderr, i.split("\t")
        options.header=False
        continue
    arr=re.split("\t",i)
#    print >> sys.stderr, arr
    count=0
    for j in range(start,end):
        # if count is nonzero add the information in column location (e.g. coverage)
        if (location!=-1 and int(arr[j])!=0):
#            print arr[j]+" "+arr[location]
            sumarr[count].append(float(arr[location]))
        else:
            sumarr[count].append(float(arr[j]))
        count=count+1
#    print >> sys.stderr, sumarr


printStats(sumarr,name,0)
    
