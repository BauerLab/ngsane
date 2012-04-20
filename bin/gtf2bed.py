import os,sys,re

sections=open(sys.argv[1]).read().split("transcript\t")

for s in sections[1::]:
    line=s.split("\n")
    arr=re.split("[ \t]",line[0])
    start=arr[0]
    stop=arr[1]
    name=arr[6].split('"')[1]
    score=int(float(arr[10].split('"')[1])*100)
    starts=[]
    lengths=[]
    for l in line[1:-1]:
        arr=re.split("[ \t]",l)
        starts.append(str(int(arr[3])-int(start)))
        lengths.append(str(int(arr[4])-int(arr[3])))
        chr=arr[0]

    print "%s\t%s\t%s\t%s\t%i\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (chr,start,stop,name,score,".",start,stop,"255,0,0",len(line)-2,",".join(lengths),",".join(starts))

        
    

