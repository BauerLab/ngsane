import argparse, sys, numpy, math, traceback
#import pandas as pd
#import matplotlib.pyplot as plt
#import numpy as np
#from rpy2.robjects.lib import ggplot2
#from rpy2.robjects.packages import importr
#import pandas.rpy.common as com


#from rpy2.robjects import pandas2ri
#pandas2ri.activate()

#def processStream(file, up, down, length):
#    matrix=numpy.empty(shape=(up+down))
#    counter=0
#    for f in open(file):
##        print "%i %i %i" % (counter%(up+down), counter, int(f.split("\t")[4]))
#         # take counter from file instead
#        if counter%(up+down)==0:
#            if counter!=0:
#                matrix=numpy.vstack((matrix, a))
#                #print matrix.shape
#            a=numpy.empty(shape=(up+down))
#        a[counter%(up+down)]=int(f.split("\t")[1])
#        counter+=1
#
#    matrix = numpy.delete(matrix, (0), axis=0)
#    return (numpy.mean(matrix, axis=0), numpy.std(matrix, axis=0))
                    

def processStream2(file, up, center, down, bin, length, args):
    matrix=numpy.empty(shape=(length,up+center+down))
    counter=0
    for f in open(file):
        try:
            #arr=map(float,f.split("\t"))
            arr=f.split("\t")
            #print "%i %s %s %i %s" % (counter/(up+down), arr[0], arr[1], len(arr),str(arr))
    #        if (len(arr)>2 and arr[2]=="-"):
    #       matrix[counter/(up+down),((up+down)-int(arr[0]))/bin]=float(arr[1])   
    #        else:
            matrix[counter/(up+center+down),(int(arr[0])-1)/bin]=float(arr[1])
            counter+=1
        except:
            print arr
            print traceback.format_exc()
        
    print "shape of matrix %s" % (str(matrix.shape))

    #ignore all genomic locations that have zero coverage
    if (args.ignore):
        z = numpy.all(matrix==0, axis=1)
        matrix = matrix.compress(numpy.logical_not(z), axis=0)
        print "shape of matrix %s after removal of zero coverage" % (str(matrix.shape))
    if (matrix.shape[0]==0):
        print("[ERROR] all TTS were removed for this mark, consider a less stringent filtering")
        sys.exit()

    #normalize by total reads in library
    if (args.normalize):
        # add pseudocount just in case
        if (args.normalize == 0):
            args.normalize = 0.001
        matrix=matrix*(1.0/args.normalize)
        print "normalize by %i" % (args.normalize)


    #remove outliers
    if (args.timestd):
        #std=numpy.std(matrix)
        #z = numpy.any(matrix>(std*args.TIMESSTD), axis=1)
        std=numpy.std(matrix,axis=0)
        mask=numpy.zeros(shape=(len(matrix)))
        for r in range(0,numpy.shape(matrix)[1]):
            loc=numpy.where(matrix[:,r]>(std[r]*args.timestd))
            mask[loc]=1
        matrix = matrix.compress(numpy.logical_not(mask), axis=0)
        print "shape of matrix %s after outlier removal" % (str(matrix.shape))

    #calculate median
    if (args.metric=="median"):
        print "median"
        return (numpy.median(matrix, axis=0), numpy.std(matrix, axis=0)*(1/math.sqrt(length)))
    if (args.metric=="sum"):
        print "sum"
        return (numpy.sum(matrix, axis=0), None)
    print "mean"
    return (numpy.mean(matrix, axis=0), numpy.std(matrix, axis=0)*(1/math.sqrt(length)))


#def processNumpy(file,up,down,length):
#    x=numpy.loadtxt(open(file,"rb"),delimiter="\t")
#    ndx = np.lexsort(keys=(data, id))
#    id, data = id[ndx], data[ndx]
#    df = DataFrame(dict(id=id, data=data))
#    print x 


#def plot():
#    print len(std)
#    print len(x)
#    df= pd.DataFrame({'mean':mean,'std':std,'x':x})
#    h = com.convert_to_r_dataframe(df)
#
#    grdevices = importr('grDevices')
#    print "plot "+args.image
#    grdevices.png(file=args.image, width=1300, height=1000)
#
#    p =ggplot2.ggplot(h) + \
#       ggplot2.aes_string(x = "x") + \
#       ggplot2.geom_line(ggplot2.aes_string(y="mean"),color='red') + \
#       ggplot2.geom_ribbon(ggplot2.aes_string(ymin="mean-std", ymax="mean+std"))
#    p.plot()


def writeR(metric,ste,xvalue,file,name,category,args):
    data=open(file+".txt","a")
    if (args.metric!="sum"):
        data.write("x\tmetric\tste\tmark\tcategory\n")
        for a,s,x in zip(metric,ste,xvalue):
            data.write("%i\t%g\t%g\t%s\t%s\n" % (x,a,s,name,category))
    else:
        data.write("x\tmetric\tmark\tcategory\n")
        for a,x in zip(metric,xvalue):
            data.write("%i\t%g\t%s\t%s\n" % (x,a,name,category))
        
    data.close()

def writeCode(file, label, args):
    geomribbon=""
    ddply=""
    anno=""
    if (args.metric!="sum"):
        ddply=", ste=mean(ste)"
        anno="-result$ste"
        geomribbon="geom_ribbon(aes(ymin=metric-ste, ymax=metric+ste, fill=mark), alpha=1/3) +\n"

    code=open(file+".R","w")
    if (args.center > 0 ):
        code.write(
        """
        library("ggplot2")
        library("plyr")
        result <- read.table('%s', header=T)
        result2<-ddply(.data=result, .(x,mark,category), summarize, metric=mean(metric)%s)
        pdf('%s.pdf', width=10, height=5*length(levels(result$category)))
        ggplot(result2, aes(x=x)) + %s
            geom_line(aes(y=metric, col=mark)) +
            labs(title = "Coverage Distribution", y = "%s coverage", x = "") +
            geom_vline(xintercept = 0, colour = "gray65") +
            geom_vline(xintercept = %d, colour = "gray65") +
            facet_grid(category ~ ., scales = "free_y") +
            annotate("text", x = %d, y = 0, label = "%s") + 
            scale_x_continuous(breaks=c(seq(%d, 0, 200),seq(%d,%d,200)), labels=c(seq(%d, 0, 200),seq(0,%d,200)))
        dev.off()
        """ % (file+".txt",ddply,file,geomribbon,args.metric, args.center,args.center/2, label,-args.upstream,args.center,args.center+args.downstream,-args.upstream,args.downstream))
    else:
        code.write(
        """
        library("ggplot2")
        library("plyr")
        result <- read.table('%s', header=T)
        result2<-ddply(.data=result, .(x,mark,category), summarize, metric=mean(metric)%s)
        pdf('%s.pdf', width=10, height=5*length(levels(result$category)))
        ggplot(result2, aes(x=x)) + %s
            geom_line(aes(y=metric, col=mark)) +
            labs(title = "Coverage Distribution", y = "%s coverage", x = "") +
            geom_vline(xintercept = 0, colour = "gray65") +
            facet_grid(category ~ ., scales = "free_y") +
            annotate("text", x = %d, y = 0, label = "%s")
        dev.off()
        """ % (file+".txt",ddply,file,geomribbon,args.metric, 0, label))    
    code.close()
    
#        annotate("text", x = %d, y = 1, label = "test %s") + 

def main(args):
    # Parse the user supplied arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', dest='inputfile', help='The coverage file')
    parser.add_argument('-n', dest='name', help='Name of feature')
    parser.add_argument('-C', dest='category', help='Name of feature category')
    parser.add_argument('-g', dest='genomic', help='Name of genomic feature')
    parser.add_argument('-u', type=int, dest='upstream', help='Downstream')
    parser.add_argument('-d', type=int, dest='downstream', help='Upstream')
    parser.add_argument('-c', type=int, dest='center', help='center', default=0)
    parser.add_argument('-l', type=int, dest='length', help='length')
    parser.add_argument('-i', dest='image', help='image name')
    parser.add_argument('-o', dest='rfile', help='R file')
    parser.add_argument('--normalize', type=int, dest='normalize', help='normalize by total numbers of reads')
    parser.add_argument('--ignoreUncovered', dest='ignore', action='store_true', help='ignore features that have zero covereage')
    parser.add_argument('--metric', dest='metric', default='mean', type=str, choices=["mean", "median", "sum"], help='metric')
    parser.add_argument('--removeoutlier', dest='timestd', type=int, help='remove outlier larger than <int>*std')
    parser.add_argument('--bin', dest='bin', type=int, help='binsize', default=1)
    args = parser.parse_args()

    upstream=args.upstream/args.bin
    downstream=(args.downstream+1)/args.bin

    #print "%i %i" % (upstream, downstream)

    #processNumpy(args.inputfile,args.upstream,args.downstream,args.length)
    #mean,std=processStream(args.inputfile,args.upstream,args.downstream,args.length)
    if (args.inputfile):
        metric,ste=processStream2(args.inputfile,upstream,args.center,downstream,args.bin,args.length, args)
        x=range(-(args.upstream),args.center+args.downstream,args.bin)
        writeR(metric,ste,x,args.rfile,args.name,args.category,args)
    if (args.image) :
        writeCode(args.rfile,args.genomic,args)
    #print mean
    #print std






#stderrd <- function(x) return (sqrt(var(x)/length(x)))
#
#result <- read.table("R/transcriptLengthFull.ggplot", header=T, quote="\"")
#result2<-ddply(.data=result, .(pos,method,status,expression), summarize, Q1=quantile(coverage,probs=c(0.25)), Q3=quantile(coverage,probs=c(0.75)),median=median(coverage))
#result2$expression <- factor(result2$expression, levels = c("high","medium","low"))
#head(result2)
#ggplot(result2, aes(x=pos)) + geom_ribbon(aes(ymin=Q1, ymax=Q3, fill=method), alpha=1/3) + 
#  geom_line(aes(y=median, col=method))+
#  facet_grid(expression~status, scales = "free_y") +
#  labs(title = "Transcript Length", y = "median coverage", x = " transcipt length (%)")

if __name__ == '__main__':
    main(sys.argv)
