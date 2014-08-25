NGSANE := $(shell trigger.sh -v 2>/dev/null)

all: checkNGSANE result/SummaryCHIPSEQ.html result/SummaryTOPHATCUFFHTSEQ.html  result/SummaryVARCALLS.html result/SummaryHIC.html 
# result/SummaryRNASEQ.html

checkNGSANE:
ifdef NGSANE
	mkdir -p tmp
	@echo "Found version $(NGSANE)"
else
	@echo NGSANE Not found
	exit -2
endif


result/SummaryBWA.html:
	./bin/testRUNBWA.qsub

result/SummaryRNASEQ.html:
	@echo "##############################################"
	./bin/testRUNRNASEQ.qsub

result/SummaryTOPHATCUFFHTSEQ.html:
	@echo "##############################################"
	./bin/testRUNTOPHATCUFFHTSEQ.qsub

result/SummaryBOWTIE.html:
	./bin/testRUNBOWTIE.qsub

result/SummaryBOWTIE2.html:
	./bin/testRUNBOWTIE2.qsub

result/SummaryVARCALLS.html:
	 @echo "##############################################"
	./bin/testRUNVARCALLS.qsub

result/SummaryCHIPSEQ.html:
	@echo "##############################################"
	./bin/testRUNCHIPSEQ.qsub

result/SummaryHIC.html:
	@echo "##############################################"
	./bin/testRUNHIC.qsub
	
test:
	@chmod a+r result/*
	@chmod 777 qout/*
	@chmod a+r qout/*/*
	@wc -l result/*.txt | grep -v total | gawk '{if ($$1==0) {print "PASS "$$2} else {print "FAIL "$$2}}'

clean:
	rm -f result/*.html
	rm -f result/diff*.txt
	rm -f result/*.tar.gz
	rm -fr Transcript ChIPseq_TF ChIPseq_Histone ChIPseq_input HiC
	rm -fr variant
	rm -fr qout
	rm -fr tmp/*
	rm -fr runStats
	rm -f igv.*
	rm -f core*
	rm -fr common
	rm -rf includes
	rm -f referenceData/chr16.*.bwt
	rm -f referenceData/chr16.*.ebwt
	rm -f referenceData/chr16.*.bt2
	rm -f referenceData/chr16_HiCcopy.*
	rm -f referenceData/chr16_RNAseqcopy.*
	rm -f hs_err_pid*
	rm -rf A1 A2 A3 B1 B2 B3 RNA3kChr16 RNA3kChr16rep 
