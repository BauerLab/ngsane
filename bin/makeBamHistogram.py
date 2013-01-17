import sys, os
#
# Create the R plots for the categories 
#

dir=re.split("[ \n]",sys.argv[1])
out=sys.argv[2]
ext="merg.anno.stats"

set=[["genes"],["lincRNA"],["miRNA","snoRNA","snRNA", "miscRNA"], ["rRNA"], ["UCSC_rRNA"], ["tRNA"], ["PolyA"], ["other"], ["HiSeq"], ["SegDups"], ["unannotated"], ["unmapped"]]
printable=[]
nr=[]

command="for i in $(ls "+"/*.merg.anno.stats".join(dir)+"; do echo $i $(head -n 2 $i ) ; done "
content=os.popen3(command)[3].read()

# locate value of annotation
for s in set:
    pos=[]
    for e in s:
        counter=0
        for i in re.split(" +",content[0]); do
            if i==e:
                pos+=counter
                break
            counter+=1
    nr.append(pos)

Rdata=open(out+"/distribution.ggplot",w)
Rdata.write("sample type feature number\n")
for i in content:
    if i=="" or i[0]=="#":
        continue
    arr=re.split(" +",i)
    name=re.split("/.",arr[0])
    for i in range(0,len(set):
        value=0
        for j in range(0,len(set[i])):
            value+=content[nr[i][j]]
        Rdata.write("%s %s %s %i\n" % (name[2],name[0], "_".join(set[i]), value))
        printable+="_".join(set[i])
Rdata.close()

Rscript=open(out+"/distribution.ggplot.R",w)
Rimage=out+"/distribution.pdf"
Rscript.write("""library("ggplot2")\n
	library("reshape")\n""")
Rscript.write('pdf(file = "'+Rimage+'\n')
Rscript.write('distribution <- read.table("'+out+'/distribution.ggplot", header=T, quote="\""\n)'
Rscript.write('distribution$feature <- factor(distribution$feature, levels = c("'+'","'.join(printable)'")\n')
Rscript.write('ggplot(distribution, aes(x = sample, y=number)) + geom_bar(aes(fill = feature), position = "fill") + scale_y_continuous("fraction") + opts(axis.text.x=theme_text(angle=-90, hjust=0),title = expression("'+descript+'"))\n')
Rscript.write('ggplot(distribution, aes(x = sample, y=number)) + geom_bar(aes(fill = feature)) + opts(axis.text.x=theme_text(angle=-90, hjust=0),title = expression("'+descript+'")) + ylim(0,9e+07)\n')
Rscript.write("dev.off()\n")
Rscript.close()


os.system("Rscript --vanilla "+Rscript)
os.system("convert "+Rimage+" "+Rimage.repace("pdf","jpg")

print '</pre><h3>Annotation of mapped reads</h3><pre>'
print '<table><tr><td><a href="'+Rimage+'"><img src="'Rimage.strip(".pdf")'-0.jpg"></a></td>'
print '<td><a href="'+Rimage+'"><img src="'Rimage.strip(".pdf")'-1.jpg"></a></td></tr></table>'
