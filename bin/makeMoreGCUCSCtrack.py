

file="/clusterdata/hiseq_apps/bin/freeze001/annovar/humandb/hg19_gc5Base.txt"

merge=5

def merging(prev,current):
    new = prev
    if (prev[1]==current[1]):
        new=[prev[0],prev[1],prev[2],current[3],current[4]]
        new.append(str(int(prev[5])+int(current[5])))
        new.append(str(int(prev[6])+int(current[6])))
        new.append(prev[7])
        new.append(prev[8])
        new.append(str(min(int(prev[9]),int(current[9]))))
        new.append("100")
        new.append(str(int(prev[11])+int(current[11])))
        new.append(str(float(prev[12])+float(current[12])))
        new.append(str(float(prev[13])+float(current[13])))
    
    return new


file=open(file).read().split("\n")
i=1
while(i<len(file)-merge):
    if (file[i]==""):
        i+=1
        continue
    stop=i+merge
    new=file[i-1].split("\t")
    while (i < stop):
        current=file[i].split("\t")
        new=merging(new,current)
        i=i+1
        
    print "\t".join(new)
    i=i+1
