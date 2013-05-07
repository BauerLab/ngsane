import sys,random

# randomizing stdin with sampling (optional)

lim=-1
if(len(sys.argv)>1):
    lim=int(sys.argv[2])
    
file=sys.stdin.read().split("\n")
random.shuffle(file)

if(lim>-1):
    print "\n".join(file[0:lim])
else:
    print "\n".join(file)
