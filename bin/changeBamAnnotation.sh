###########
# change the chromosome annotation from chr -> 1 
#
# ./changeAnnotation.sh BAM_IN BAM_OUT 
#
#  HWI-ST212:243:D0W6RACXX:3:2316:17108:100013     161     _chr_X    71492464        37      100M    _chr_10   32391455        0       CTGCCACAATATTTTTAATTACGTACAAAGATCTGACATGTCACCCAGGGACCCATTTCACCCACTGCTCTGTTTGGCCGCCAGTCGTTTGTCTCTCTCT B@CFFFDFFHHHHJJJI9FHIIIGJ;FHIJIJJJJJJJJJJJJJ@HHIIJJJIJJGHJIJJJJJJJJJJHHHHHHFEFFD?/9=@5',8?BDCD3>@CDD RG:Z:OmicsRNA   XT:A:U  NM:i:1  SM:i:37 AM:i:23 X0:i:1  X1:i:0  XM:i:1  XO:i:0  XG:i:0  MD:Z:86T13


samtools view -h $1 |\
gawk  'OFS="\t" {if ($1~"@SQ"){split($2,arr,":");$2=arr[1]":chr"arr[2]} 
                  else{
                  if (!($1~"@")){
                  if (!(($3=="*" || $3=="="))){$3="chr"$3}
                  if (!(($7=="*" || $7=="="))){$7="chr"$7}}}; 
                  print $0}' | samtools view -b -S - >$2 