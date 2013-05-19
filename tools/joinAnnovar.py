import os,sys,re

#"0/0:23,0:23:66.21:0,66,817"
def format(str):
    st=str.split(":")
    if (st[0]=="./.\",\""):
        st[0]="./."
    return "'"+st[0]+"'\",\""+":".join(st[1::])

#Robyn.filter.snps.Robyns_1.sum.exome_summary.csv
def format2(str):
    return str.split(".")[3]

i=1
files=[]
while (len(sys.argv)>i):
    files.append(sys.argv[i])
    i+=1

command="paste -d $',' "+" ".join(files)
result=os.popen3(command)[1].read().split("\n")

header=result[0].split("Otherinfo")[0]+"CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,"
for f in files:
    header+=format2(f)+","
    header+=format2(f)+"_details,"
print header[0:-1]

for r in result[1:-1]:
    arr=re.split("GT:AD:DP",r)

    line=arr[0]#",".join(arr[0].split(","))[1::]
    for i in range(0,len(files)):
        line+=format(arr[i+1].split("\",\"")[1]+"\",\"")
    print line[0:-3]
