import sys

# add a column indicating whehter the element is not annotated by bedanno
# c s    e     l c gene rib miRNA *unannotated*
# 1 1111 11112 3 4    1   1     1            0
# 1 1112 11113 3 4    0   0     0            1

c=0
for i in sys.stdin:
    i=i.strip("\n")
    i=i.rstrip("\t")
 #   if i.find("gene")>-1:
    if i[0]=="#":
#    if (c==0):
        print i+"\tunannotated"
#        c+=1
        continue
    arr=i.split("\t")
#    print arr[5::]
    # from 6th column
    if sum(map(float,arr[5::]))==0:
        arr.append("1")
    else:
        arr.append("0")
    print "\t".join(arr)
