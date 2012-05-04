# summarize all read quality plots on one page

SOURCE=$1
OUTFILE=$3
OUTDIR=$2

LATEX=$OUTDIR/${OUTFILE/.pdf/.tex}
echo $LATEX

number=`ls -a $SOURCE/*read1*/Images/per_base_quality.png | wc -l`
let size=$number*2

echo """\documentclass{article}
\usepackage[margin=0.3in, paperwidth=8.5in, paperheight=$size in]{geometry}
\usepackage{graphicx}


\begin{document}
""" >$LATEX


echo """\begin{figure}[!ht]
           \centering
           \begin{tabular}{|c|c|}  
           \hline""" >>$LATEX


for i in $( ls -a $SOURCE/*read1*/Images/per_base_quality.png); do
    name=${i/$SOURCE\//}
    name=${name/_fastqc\/Images\/per_base_quality.png/}
    name=${name//_/"\_"}
    if [ -e ${i/read1/read2} ]; then
	echo $i
	echo $name" & "${name/read1/read2}"\\\\">>$LATEX
	echo "\includegraphics[height=1.7in,width=2.3in,type=png,ext=.png,read=.png]{"${i/.png/}"}&" >>$LATEX
	i2=${i/read1/read2}
	echo "\includegraphics[height=1.7in,width=2.3in,type=png,ext=.png,read=.png]{"${i2/.png/}"}\\\\" >>$LATEX
    else
	echo $name "& \\\\" >>$LATEX
	echo "\includegraphics[height=1.5in,width=2in,type=png,ext=.png,read=.png]{"${i/.png/}"}& \\\\" >>$LATEX

    fi

done


echo """\hline
           \end{tabular} 
           \end{figure}""" >>$LATEX
echo "\end{document}" >>$LATEX


cd $OUTDIR
/old_tools/new_apps/texlive/2012/bin/x86_64-linux/pdflatex ${OUTFILE/pdf/tex}
rm ${OUTFILE/pdf/aux}
rm ${OUTFILE/pdf/log}
#rm ${OUTFILE/pdf/tex}
convert -density 1000 -geometry 800x20000 $OUTFILE ${OUTFILE/pdf/jpg} 
#-scale 2500x1500