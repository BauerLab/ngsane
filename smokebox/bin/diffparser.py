import os, re, sys

def parser(l,r):
    rn=re.findall("[ \"'>][0-9]+[0-9.]*[ \"'<]",r)
    ln=re.findall("[ \"'>][0-9]+[0-9.]*[ \"'<]",l)
#    print "ln %s" % (ln)
#    print "rn %s" % (rn)
    for i,j in zip(rn,ln):
        if (i[0]=="." or j[0]==","):
#            print "OK %s vs %s with context\n%svs\n%s-----REST------" % (i,j,l,r)
            continue
        elif(float(i)*0.8>float(j) or float(i)*1.2<float(j)):
            print "Too different %s vs %s with context\n%svs\n%s-----REST------" % (i,j,l,r)
            return -1
    return 0

line=sys.stdin
local=[]
remote=[]
content=""
for line in sys.stdin:
    content+=line
    if line.find("Last modi")>-1 or line.find("NGSANE")>-1 or line.find("[NOTE]")>-1 or line.find("class='citation'")>-1:
        continue
    if line[0]=="<" :
        local.append(line)
    if line[0]==">" :
        remote.append(line)

if len(local) != len(remote):
    print "local and remote changed lines differ (%d != %d)" % (len(local), len(remote))
    print content
    sys.exit(1)

for i, j in zip(local, remote):
    if re.match("[0-9]",i) or re.match("[0-9]",j) or parser(i,j)==-1:
        print content
        sys.exit(1)
