import os,sys,re

def split(string):
    result=[]
    i=0
    arr=string.split(",")
    join=False
    while i<len(arr):
        if (join):
            if (arr[i][-1]=="\""):
                 join=False
                 #add rest of the quoted but comma separated stuff
                 result[-1]+=","+(arr[i])
            else:
                result[-1]+=","+(arr[i])
        else:
            # start a run if the string has not ended in this element
            if (arr[i]!="" and arr[i][0]=="\"" and not(arr[i][-1]=="\"")):
                join=True
            result.append(arr[i])
        i+=1
    return result

#"0/1:11,0:11:1.48:0,4,18"
def format(str):
    string="'./.',"
    if (str.find("./.")==-1):
        arr=str.split(":")
        string="'"+arr[0][1::]+"',\""+":".join(arr[1::])
    return string

file=open(sys.argv[1]).read().split("\n")
names=sys.argv[2::]

header=file[0].split("Otherinfo")[0]+",".join(names[0:9])
for n in names[9::]:
    header+=","+n+","+n+"_details"
print header

for r in file[1:-1]:
    arr=split(r)
    line=",".join(arr[0:29])
    for a in arr[29::]:
        line+=","+format(a)
    print line
