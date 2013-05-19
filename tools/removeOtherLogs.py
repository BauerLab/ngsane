import sys,os

dir=sys.argv[1]
objects=os.listdir(dir)

for o in objects:
    print o
    f=open(dir+"/"+o).read()
    file=f.split(">>>>> enddate")
    if(len(file)>1):
        content="\n".join(file[-2].split("\n")[1::])+">>>>> enddate"+file[-1]
    else:
        content=f
    out=open(dir+"/"+o,"w")
    out.write(content)
    out.close()
