import os, re, sys

def parser(l,r):
    rn=re.findall("[0-9.]+",r)
    ln=re.findall("[0-9.,]+",l)
#    print rn
#    print ln
    for i,j in zip(rn,ln):
        if (i[0]=="." or j[0]==","):
            continue
        elif(float(i)*0.8>float(j) or float(i)*1.2<float(j)):
            print "Too different %s %s" % (i,j)
            return -1
    return 0

line=sys.stdin
local=[]
remote=[]
content=""
for line in sys.stdin:
    content+=line
    if line.find("Last modi")>-1 or line.find("NGSANE")>-1 or line.find("[NOTE]")>-1 :
        continue
#| grep -v "\-\-" | grep -v
# "NGSANE" |       grep -v "[NOTE]" |       grep -v "[0-9]c[0-9]")
    if line[0]=="<" :
        local.append(line)
    if line[0]==">" :
        remote.append(line)

#print len(local)
#print len(remote)

for i, j in zip(local, remote):
    if re.match("[0-9]",i) or re.match("[0-9]",j) or parser(i,j)==-1:
        print content
        sys.exit(1)
