import sys, os,re

def comment(conf,setFile):
    os.system("cp "+conf+" "+conf+".old")
    if (open(conf).read().find("<--")>-1):
        remove(conf)
    file=open(conf).read().split("\n")
    tileSet=set(open(setFile).read().split(">,")[1].split("\n")[0].split(","))
    out=open(conf,"w")
    print "Mask: "+str(tileSet)
    w=1
    lane=0
    while (w <len(file)):
        comment=False
        tile=-1
        if (file[w].find("Lane Index")>-1):
            lane=file[w].split("\"")[1]
        if (file[w].find("<Tile>")>-1):
            tile="%04i" % int(re.split("[<>]",file[w])[2])
            if ("s_"+lane+"_"+tile in tileSet):
                comment=True
        if (comment):
            out.write("    <!--%s-->" % file[w].strip()+"\n")
        else:
            out.write(file[w].strip("\M")+"\n")
        w+=1
    out.close



def remove(conf):
    os.system("cp "+conf+" "+conf+".old")
    file=open(conf).read().split("\n")
    out=open(conf,"w")
    for f in file:
        if (f.find("<!--")>-1):
            f="       "+f.split("<!--")[1].split("-->")[0]
        out.write(f+"\n")
    out.close()


config=sys.argv[1]
task=sys.argv[2]
if(task=="--set"):
    comment(config,sys.argv[3])
elif (task=="--remove"):
    remove(config)
else:
    print "dont understand "+task

    
