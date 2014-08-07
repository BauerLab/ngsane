import sys

filename=sys.argv[1]
outfile=open(sys.argv[2],"w")

countall=0
count=0 


for i in open(filename):
    if i.find("#")>-1:
        outfile.write(i)
    else:
        countall+=1
        content=i.split("\t")
        start=int(content[1])
        end=int(i.split("END=")[1].split(";")[0])
        length=len(content[3])
        if (start+length-1!=end):
            count+=1
            #print "update length %i %i : %s" % (start+length-1, end, i)
            content[7]=content[7].replace("END="+str(end), "END="+str(start+length-1))
        outfile.write("\t".join(content))

print "%i updated out of %i" % (count, countall)

        
        
