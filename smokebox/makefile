#all: result/SummaryTOPHATCUFF.html result/SummaryBWA.html result/SummaryBOWTIE.html result/SummaryBOWTIE2.html result/SummaryVARCALL.html

NGSANE := $(shell trigger.sh -v 2>/dev/null)

all: checkNGSANE result/SummaryTOPHATCUFF.html result/SummaryCHIPSEQ.html #result/SummaryVARCALL.html

checkNGSANE:
ifdef NGSANE
	@echo "Found version $(NGSANE)"
else
	@echo NGSANE Not found
endif

result/SummaryBWA.html:
	./bin/testRUNBWA.qsub

result/SummaryTOPHATCUFF.html:
#	[ -n $$(hash trigger.sh) ] && echo "[ERROR] NGSANE not loaded" && exit -1
	./bin/testRUNTOPHATCUFF.qsub

result/SummaryBOWTIE.html:
	./bin/testRUNBOWTIE.qsub

result/SummaryBOWTIE2.html:
	./bin/testRUNBOWTIE2.qsub

result/SummaryVARCALL.html:
	./bin/testRUNVARCALL.qsub

result/SummaryCHIPSEQ.html:
	./bin/testRUNCHIPSEQ.qsub

test:
	chmod a+r result/*
	chmod 777 qout/*
	chmod a+r qout/*/*
	wc -l result/*.txt | head -n -1 | gawk '{if ($$1==0) {print "PASS "$$2} else {print "FAIL "$$2}}'

clean:
	rm -f result/*.html
	rm -f result/diff*.txt
	rm -fr Transcript ChIPseq ChIPseq_input
	rm -fr variant
	rm -fr qout
	rm -f tmp/*
