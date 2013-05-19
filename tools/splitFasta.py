import sys

# Splits fasta entries into seperate files
# >f1
# ...
# >f2
# ..
# becomes two files f1.fa and f2.fa
# author Denis C. Bauer


file=open(sys.argv[1])
dir=sys.argv[2]

outfile=""

while 1:
    f = file.readline()
    if not f:
        break
    pass
    if (f[0]==">"):
        name=f.strip(">").strip("\n")
        print name
        if(outfile==""):
            outfile=open(dir+name+".fa","w")
        else:
            outfile.close()
            outfile=open(dir+name+".fa","w")
    outfile.write(f)


print "done"
            
            
