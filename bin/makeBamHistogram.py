import sys, os, re, subprocess, commands
#
# Create the R plots for the categories 
#

dirname=re.split("[ \n]",sys.argv[1])
out=sys.argv[2]

setname=[["Pgenes"],["lincRNA"],["miRNA","snoRNA","snRNA", "miscRNA"], ["rRNA"], ["ucsc_rRNA"], ["PolyA"], ["tRNA","other"], ["HiSeq", "segDups"], ["unannotated"], ["unmapped"]]
printable=[]
nr=[]

content=[]
for d in dirname:
    if d=="":
        continue
    content+=commands.getoutput("head -n 2 "+d+"/*.merg.anno.stats").split("==>")

# locate value of annotation
annotations=0
for s in setname:
    pos=[]
    for e in s:
        counter=0
        annotations+=1
        for i in re.split("[ \n]+",content[1]):
            if i==e:
                pos.append(counter)
                break
            counter+=1
    nr.append(pos)

#print setname
#print nr

descript=""
printable=[0]*len(setname)
Rdata=open(out+"/distribution.ggplot","w")
Rdata.write("sample type feature number\n")
for line in content:
    arr=re.split("[ \n]", line)
    if line=="" or line[0]=="#":
        continue
    name=[re.split("/",arr[1])[-3],re.split("/",arr[1])[-1].split(".")[0]]
    descript=re.split("/",arr[1])[-4]
#    name=re.split("[/.]",a[0])
    for i in range(0,len(setname)):
        value=0
        for j in range(0,len(setname[i])):
#            print nr[i][j]
#            print arr[nr[i][j]+1+annotations]
#            print arr
            try:
                value+=int(arr[nr[i][j]+1+annotations])
            except:
                value+=0
#        print name
        Rdata.write("%s %s %s %i\n" % (name[1],name[0], "_".join(setname[i]), value))
        printable[i]="_".join(setname[i])

Rdata.close()



Rscript=open(out+"/distribution.ggplot.R","w")
Rimage=out+"/distribution.pdf"
Rscript.write('library("ggplot2")\nlibrary("reshape")\n')
Rscript.write('pdf(file = "'+Rimage+'",width=10,height=7)\n')
Rscript.write('distribution <- read.table("'+out+'/distribution.ggplot", header=T, quote="\\"")\n')
#Rscript.write('head(distribution)\n')
Rscript.write('distribution$feature <- factor(distribution$feature, levels = c("'+'","'.join(printable)+'"))\n')
Rscript.write('ggplot(distribution, aes(x = sample, y=number)) + geom_bar(stat="identity", aes(fill = feature), position = "fill") + scale_y_continuous("fraction") + opts(axis.text.x=theme_text(angle=-90, hjust=0), title="'+descript+'")\n')
Rscript.write('ggplot(distribution, aes(x = sample, y=number)) + geom_bar(stat="identity", aes(fill = feature)) + opts(axis.text.x=theme_text(angle=-90, hjust=0), title="'+descript+'")\n')
#+ ylim(0,17e+07)\n')
Rscript.write("dev.off()\n")
Rscript.close()


x=os.system("Rscript --vanilla "+out+"/distribution.ggplot.R")
x=os.system("convert "+Rimage+" "+Rimage.replace("pdf","jpg"))

print "</pre><h3>Annotation of mapped reads</h3><pre>"
print '<table><tr><td><a href="'+Rimage+'"><img src="'+Rimage.strip(".pdf")+'-0.jpg"></a></td>'
print '<td><a href="'+Rimage+'"><img src="'+Rimage.strip(".pdf")+'-1.jpg"></a></td></tr></table>'
