#all: checkNGSANE result/SummaryVARCALL.html

NGSANE := $(shell trigger.sh -v 2>/dev/null)

all: checkNGSANE result/SummaryCHIPSEQ.html  result/SummaryTOPHATCUFFHTSEQ.html  result/SummaryVARCALL.html

checkNGSANE:
ifdef NGSANE
	@echo "Found version $(NGSANE)"
	mkdir -p tmp
else
	@echo NGSANE Not found
	exit -2
endif

result/SummaryBWA.html:
	./bin/testRUNBWA.qsub

result/SummaryTOPHATCUFFHTSEQ.html:
	@echo "##############################################"
	./bin/testRUNTOPHATCUFFHTSEQ.qsub

result/SummaryBOWTIE.html:
	./bin/testRUNBOWTIE.qsub

result/SummaryBOWTIE2.html:
	./bin/testRUNBOWTIE2.qsub

result/SummaryVARCALL.html:
	 @echo "##############################################"
	./bin/testRUNVARCALL.qsub

result/SummaryCHIPSEQ.html:
	@echo "##############################################"
	./bin/testRUNCHIPSEQ.qsub

test:
	@chmod a+r result/*
	@chmod 777 qout/*
	@chmod a+r qout/*/*
	@wc -l result/*.txt | head -n -1 | gawk '{if ($$1==0) {print "PASS "$$2} else {print "FAIL "$$2}}'

clean:
	rm -f result/*.html
	rm -f result/diff*.txt
	rm -fr Transcript ChIPseq_TF ChIPseq_Histone ChIPseq_input
	rm -fr variant
	rm -fr qout
	rm -fr tmp/*
	rm -fr runStats
