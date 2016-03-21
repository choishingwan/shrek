SHELL := /bin/bash
#HOME := $(shell pwd)
HOME := /home/sam/workspace/inheritance_estimate/simulation/WTCCC
RELATED := 0.05
THREADS := 10
GENO := 0.01
MAF := 0.05
MIND := 0.95
HWE := 1e-5
NSET:= 10
SET :=0
TYPE :=0
PLINK :=$(HOME)/script/plink
GCTA :=$(HOME)/script/gcta64
SHREK := $(HOME)/script/shrek
PERM := 2  
CHR := 22
NSAMPLE := 1000
RSAMPLE := 500

.PHONY: help getSample processSample simulation groupSample cleanInter run
.DEFAULT_GOAL:= help
genProgramme: ## Compile the C++ programmes in the script folder
	@cd $(HOME)/script;\
	make -s all;

getSample: ## Get the required samples for the simulation
	@mkdir -p $(HOME)/samples;
	@cd samples && find /psychipc01/disk2/references/WTCCC/ -maxdepth 1 -name "[^W]*[f|b][i|a|e][d|m]" -exec echo ln -sf {} . \; | bash;\
	for i in `ls /psychipc01/disk2/references/WTCCC/| egrep "^[^W]*fam"`; do \
	disease=($${i//./ });\
	awk -v fam=$$disease '{print fam" "$$2" "$$3" "$$4" "$$5" 1"}' $$disease.fam > $$disease.temp.fam;\
	rm $$disease.fam;\
	mv $$disease.temp.fam $$disease.fam;\
	done;
processSample: QCLog genProgramme getSample ## Preprocess the samples such that they are clean before performing the analysis
	@related=$(RELATED);\
	cd $(HOME)/samples;\
	for i in `ls /psychipc01/disk2/references/WTCCC/| egrep "^[^W]*fam"`; do\
	    disease=($${i//./ });\
	    echo "Working on $$disease";\
	    $(PLINK) --bfile $$disease --geno $(GENO) --mind $(MIND) --maf $(MAF) --make-bed --out $$disease.temp --threads $(THREADS) --indep 50 5 2 > /dev/null 2>&1;\
	    echo "Prune the samples for subsequent processing";\
	    $(PLINK) --bfile $$disease.temp --extract $$disease.temp.prune.in --ibc --out $$disease.temp --threads $(THREADS) >/dev/null 2>&1;\
	    Rscript $(HOME)/script/heteroCheck.R $$disease.temp.ibc $$disease.temp.remove >/dev/null 2>&1;\
	    temp=$$(wc -l $$disease.temp.remove | cut -f 1 -d " ");\
	    echo "$${temp[0]} samples removed after heterozygousity test";\
	    $(PLINK) --bfile $$disease.temp --remove $$disease.temp.remove --extract $$disease.temp.prune.in --genome --min $$related --threads $(THREADS) --out $$disease.temp >/dev/null 2>&1;\
	    $(HOME)/script/GreedyRelativeness $$disease.temp.genome >> $$disease.temp.remove;\
	    temp2=$$(wc -l $$disease.temp.remove | cut -f 1 -d " ");\
	    echo "$$((temp2-temp)) samples removed after considering genetic relationship";\
	    $(PLINK) --bfile $$disease.temp --remove $$disease.temp.remove --extract $$disease.temp.prune.in --threads $(THREADS) --pca --out $$disease.temp >/dev/null 2>&1;\
	    Rscript $(HOME)/script/plotPCA.R $$disease.temp.eigenvec $$disease.temp >/dev/null 2>&1;\
	    temp=$$(wc -l $$disease.temp.pc | cut -f 1 -d " ");\
	    echo "$${temp[0]} samples seems to be outlier considering PCA plot";\
	    NUMOFLINES=$${temp[0]};\
	    while [ $$NUMOFLINES -gt 0 ]; do \
		cat $$disease.temp.remove $$disease.temp.pc > $$disease.temp.remove.temp;\
		mv $$disease.temp.remove.temp $$disease.temp.remove;\
		rm $$disease.temp.pc;\
		$(PLINK) --bfile $$disease.temp --remove $$disease.temp.remove --extract $$disease.temp.prune.in --threads $(THREADS) --pca --out $$disease.temp>/dev/null 2>&1;\
		Rscript $(HOME)/script/plotPCA.R $$disease.temp.eigenvec $$disease.temp >/dev/null 2>&1;\
		temp=$$(wc -l $$disease.temp.pc | cut -f 1 -d " ");\
		NUMOFLINES=$${temp[0]};\
		echo "Addition $$NUMOFLINES samples seems to be outlier considering PCA plot";\
	    done;\
	    $(PLINK) --bfile $$disease.temp --remove $$disease.temp.remove --geno $(GENO) --mind $(MIND) --maf $(MAF) --hardy --out $$disease.temp.filter --make-bed >/dev/null 2>&1;\
            temp=$$(wc -l $$disease.temp.filter.fam | cut -f 1 -d " ");\
	    temp2=$$(wc -l $$disease.temp.filter.bim | cut -f 1 -d " ");\
            echo "In total, $$temp samples and $$temp2 SNPs were retained";\
	    rm -f *temp.[^fp][^i]*;\
	    rm -f *temp.?[^in]*;\
	done;\
	
groupSample: genProgramme ## Group the processed samples to form the final sample
	@count=1;\
	cd $(HOME)/samples;\
	related=$(RELATED);\
	firstDisease="";\
	rm -f sampleList;\
	for i in `ls *.temp.filter.fam`; do \
	disease=($${i//./ });\
	if [ $$count -eq 1 ]; then \
	    firstDisease=$$disease;\
	    count=2;\
	fi;\
	echo $$disease.temp.filter.bed $$disease.temp.filter.bim $$disease.temp.filter.fam >> sampleList;\
	done;\
	tail -n +1 sampleList > samples;\
	rm sampleList;\
	echo "Merging samples";\
	$(PLINK) --bfile $$firstDisease.temp.filter --make-bed --out WTCCC17000.filter --merge-list samples --geno $(GENO) --mind $(MIND) --maf $(MAF) --hwe $(HWE) --indep 50 5 2 --threads $(THREADS) >/dev/null 2>&1; \
	echo "Start checking for relatedness in the samples";\
	$(PLINK) --bfile WTCCC17000.filter --extract WTCCC17000.filter.prune.in --genome --min $$related --out WTCCC17000.filter --threads $(THREADS) >/dev/null 2>&1; \
	$(HOME)/script/GreedyRelativeness WTCCC17000.filter.genome > WTCCC17000.filter.remove; \
	echo "Generate the final sample";\
	$(PLINK) --bfile WTCCC17000.filter --remove WTCCC17000.filter.remove --threads $(THREADS) --geno $(GENO) --mind $(MIND) --maf $(MAF) --hwe $(HWE) --make-bed --out WTCCC17000.final >/dev/null 2>&1;\
	temp=$$(wc -l WTCCC17000.final.fam | cut -f 1 -d " ");\
	temp2=$$(wc -l WTCCC17000.final.bim | cut -f 1 -d " ");\
	echo "In total, $$temp samples and $$temp2 SNPs were retained";\
	cp WTCCC17000.final.fam WTCCC17000.final.bk.fam;\
	#plink --bfile WTCCC17000.final --pca --extract WTCCC17000.filter.prune.in --out WTCCC17000.final --threads $(THREADS);\
	#Rscript $(HOME)/script/plotPCA.R WTCCC17000.final.eigenvec WTCCC17000.final;\
	#echo "Ready for simulation";

cleanInter: ## Cleaning the intermediate samples. 
	@cd samples;\
	rm -f *temp.[^p]*;\
	rm -f *.fi[^n]*;\
	rm -f $$(ls * | grep -v final | grep -v png);
	@echo "Intermediate files all clean";

simulation: genProgramme ## Generate the GRM and also other required files for the simulation. Will automatically continue to next step.
	@echo "Start running the simulation";
	@herit=($$(seq -s ' ' 0 0.1 0.9));\
	heritIndex=($$(seq -s ' ' 0 1 9));\
	causal=( 10 50 100 500 1000 );\
	prevalence=( 0.5 0.25 0.1);\
	observePrev=0.5;\
	echo "Simulation Parameters:";\
	echo "Heritability:          $${herit[@]}";\
	echo "Number of causal SNPs: $${causal[@]}";\
	echo "Population Prevalence: $${prevalence[@]}";\
	echo "Observed Prevalence:   $$observePrev";\
	echo "Sample Size:           $(NSAMPLE)";\
	echo "Size of reference:     $(RSAMPLE)";\
	echo "Sets of Simulation:    $(NSET)";\
	echo "Number of permutation: $(PERM)";\
	echo "Chromosome:            $(CHR)";\
	echo "===================================";\
	echo "Start Phenotype Simulation";\
	mkdir -p simulated;\
	$(PLINK) --bfile $(HOME)/samples/WTCCC17000.final --chr $(CHR) --make-bed --out $(HOME)/simulated/WTCCC17000.final.chr$(CHR)>/dev/null 2>&1;\
	echo $(HOME)/script/SimulatePheno $(HOME)/simulated/WTCCC17000.final.chr$(CHR) $(NSET) \"$${causal[*]}\" | bash ;\
	echo "Phenotype Simulated";\
	echo "Submit job to caclulate GRM";\
	grmJob=$$(echo "$(GCTA) --bfile $(HOME)/simulated/WTCCC17000.final.chr$(CHR) --make-grm-bin --autosome --maf $(MAF) --out $(HOME)/simulated/WTCCC17000.final.chr$(CHR).gcta --thread-num $(THREADS)" | qsub -l nodes=1:ppn=$(THREADS) -l mem=10gb -l walltime=12:00:00 -N GRM -q default -V -d $(HOME) );\
	echo "Initialize Simulation";\
	cd $(HOME);\
	for i in $${heritIndex[@]}; do \
	echo "make -s preparePhenotype SET=$$i" | qsub -W depend=afterok:$$grmJob -l nodes=1:ppn=1 -l mem=8gb -l walltime=4:00:00 -N sim_set$$i -V -d $(HOME) > /dev/null 2>&1;\
	(crontab -u sam -l 2>/dev/null ; echo "*/5 * * * * cd $(HOME) && make plinkWait SET=$$i ") | crontab -u sam -;\
	(crontab -u sam -l 2>/dev/null ; echo "*/10 * * * * cd $(HOME) && make heritWait SET=$$i ") | crontab -u sam  -;\
	done;\
	echo "Start $${#herit[@]} jobs";\

preparePhenotype: ## Simulate all the phenotypes and run the plink association afterward
	@herit=($$(seq -s ' ' 0 0.1 0.9));\
	causal=( 10 50 100 500 1000 );\
	prevalence=( 0.5 0.25 0.1);\
	observePrev=0.5;\
	echo "Start performing simulation";\
	echo "A total of $(PERM) repeats will be performed";\
	echo "Generating the pheno file for the simulation"; \
	mkdir -p $(HOME)/simulated/set$(SET);\
	echo "Rscript $(HOME)/script/phenoProcess.R $(HOME)/simulated/WTCCC17000.final.chr$(CHR).fam $(HOME)/simulated/WTCCC17000.final.chr$(CHR).pheno $(HOME)/simulated/set$(SET)/WTCCC17000.final.chr$(CHR) $(SET) $(PERM) $(NSAMPLE) $(RSAMPLE) \"$${observePrev[@]}\" \"$${prevalence[@]}\" \"$${herit[@]}\"" \"$${causal[@]}\" | bash > /dev/null 2>&1;\
	echo "Extract the reference sample for LDSC and SHREK";\
	$(PLINK) --bfile $(HOME)/simulated/WTCCC17000.final.chr$(CHR) --keep $(HOME)/simulated/set$(SET)/WTCCC17000.final.chr$(CHR).ref --make-bed --out $(HOME)/simulated/set$(SET)/WTCCC17000.final.chr$(CHR).ref > /dev/null 2>&1;\
	echo "Calculate the LD Score";\
	python $(HOME)/script/ldsc.py --bfile $(HOME)/simulated/set$(SET)/WTCCC17000.final.chr$(CHR).ref --l2 --ld-wind-kb 1000 --out $(HOME)/simulated/set$(SET)/WTCCC17000.final.chr$(CHR).ldsc > /dev/null 2>&1;\
	echo "Ready for association testing";\
	mkdir -p $(HOME)/simulated/set$(SET)/QT;\
	mkdir -p $(HOME)/simulated/set$(SET)/CC;\
	mkdir -p $(HOME)/simulated/set$(SET)/EX;\
	touch $(HOME)/.set$(SET)pheno;	

plinkWait: ## Process to distribute the plink association jobs.
	@#First, check if there are empty locations for job to be submitted
	@#It is important to make sure we do not access the crontab together. So use the lock file trick
	@#If there is no space to run, then just wait for 5 minutes and we will rerun this script
	@if [ ! -e $(HOME)/.set$(SET)pheno ]; then \
	    exit;\
	else \
	    while ! (set -o noclobber ; echo > $(HOME)/.lock) 2>/dev/null; do \
		sleep 10;\
	    done; \
	    echo "Plink wait $(SET)"; \
	    currentJob=$$( qstat -an1 | awk '{if(($$10=="R" || $$10=="Q" || $$10=="H")&&$$3=="small") print $$0}' |wc -l);\
	    if [ $$currentJob -lt 27 ]; then \
		echo "Can initialize job";\
		(crontab -l 2> /dev/null | grep -v "plinkWait SET=$$SET ") | crontab -; \
		echo "make -s plinkRun SET=$(SET) TYPE=1" | qsub -l nodes=1:ppn=1 -l mem=8gb -l walltime=4:00:00 -N sim_qt_s$(SET) -q small -V -d $(HOME) > /dev/null 2>&1 ;\
		echo "make -s plinkRun SET=$(SET) TYPE=2" | qsub -l nodes=1:ppn=1 -l mem=8gb -l walltime=4:00:00 -N sim_cc_s$(SET) -q small -V -d $(HOME) > /dev/null 2>&1;\
		echo "make -s plinkRun SET=$(SET) TYPE=3" | qsub -l nodes=1:ppn=1 -l mem=8gb -l walltime=4:00:00 -N sim_ex_s$(SET) -q small -V -d $(HOME) > /dev/null 2>&1;\
	    fi;\
	    rm $(HOME)/.lock;\
	fi;\

plinkRun: ## Run the plink association analysis.
	@#Type=1 -> QT; Type=2->CC; Type=3 -> EX;
	@prevalence=( 0.5 0.25 0.1);\
	observePrev=0.5;\
	if [[ $(TYPE) -eq 1 ]] ; then \
	set +C;\
	$(PLINK) --bfile $(HOME)/simulated/WTCCC17000.final.chr$(CHR) --assoc --pheno $(HOME)/simulated/set$(SET)/WTCCC17000.final.chr$(CHR).qt.pheno --all-pheno --out $(HOME)/simulated/set$(SET)/QT/QT ;\
	elif [[ $(TYPE) -eq 2 ]]; then \
	for i in $${prevalence[@]}; do \
	set +C;\
	$(PLINK) --bfile $(HOME)/simulated/WTCCC17000.final.chr$(CHR) --assoc --pheno $(HOME)/simulated/set$(SET)/WTCCC17000.final.chr$(CHR)_$${observePrev}_$$i.cc.pheno --all-pheno --out $(HOME)/simulated/set$(SET)/CC/CC_$${i}_ ;\
	done; \
	elif [[ $(TYPE) -eq 3 ]]; then \
	for i in $${prevalence[@]}; do \
	set +C;\
	$(PLINK) --bfile $(HOME)/simulated/WTCCC17000.final.chr$(CHR) --assoc --pheno $(HOME)/simulated/set$(SET)/WTCCC17000.final.chr$(CHR)_$$i.ex.pheno --all-pheno --out $(HOME)/simulated/set$(SET)/EX/EX_$${i}_ ;\
	done; \
	fi;\
	touch $(HOME)/.set$(SET)plink;	

heritWait: ## Process to distribute the heritability estimation jobs
	@if [ ! -e $(HOME)/.set$(SET)plink ]; then \
	    exit;\
	else \
	    while ! (set -o noclobber ; echo > $(HOME)/.lock) 2>/dev/null; do \
		Sleep 10;\
	    done;\
	    currentJob=$$( qstat -an1 | awk '{if(($$10=="R" || $$10=="Q" || $$10=="H")&&$$3=="medium") print $$0}' |wc -l);\
	    if [ currentJob -lt 17 ]; then \
		echo "Can initialize job";\
		(crontab -l 2> /dev/null | grep -v "heritWait SET=$$SET ") | crontab -; \
		echo "make -s heritRun SET=$(SET) TYPE=1" | qsub -l nodes=1:ppn=1 -l mem=20gb -l walltime=12:00:00 -N sim_sk_s$(SET) -q medium -V -d $(HOME) ;\
		echo "make -s heritRun SET=$(SET) TYPE=2" | qsub -l nodes=1:ppn=1 -l mem=20gb -l walltime=12:00:00 -N sim_ld_s$(SET) -q medium -V -d $(HOME) ;\
		echo "make -s heritRun SET=$(SET) TYPE=3" | qsub -l nodes=1:ppn=1 -l mem=20gb -l walltime=12:00:00 -N sim_gc_s$(SET) -q medium -V -d $(HOME) ;\
	    fi;\
	    rm $(HOME)/.lock; \
	fi;
 
heritRun: ## Running the heritability estimation
	@exit;
	@#1 = shrek 2 = ldsc 3 = gcta
	@if [[ $(TYPE) -eq 1 ]]; then \
	    make heritShrek SET=$(SET);\
	elif [[ $(TYPE) -eq 2 ]]; then \
	    make heritLDSC SET=$(SET);\
	elif [[ $(TYPE) -eq 3 ]]; then \
	    make heritGCTA SET=$(SET);\
	fi;\

heritShrek: ## Use Shrek to perform the heritability estimation
	
heritLDSC: ## Use LDSC to perform the heritability estimation

heritGCTA: ## Use Gcta to perform the heritability estimation
	@herit=($$(seq -s ' ' 0 0.1 0.9));\
	causal=( 10 50 100 500 1000 );\
	prevalence=( 0.5 0.25 0.1);\
	observePrev=0.5;\
	count=1;\
        causeIndex=0;\
        heritIndex=0;\
        perm=1;\
	echo "Quantitative Trait Estimation (GCTA)";\
	rm -f $(HOME)/simulated/set$(SET)/QT/gcta.res;\
	echo "Perm Causal Herit Estimate SE" > $(HOME)/simulated/set$(SET)/QT/gcta.res;\
	for i in `ls $(HOME)/simulated/set$(SET)/QT/*.qassoc`; do \
	    temp=($${i//./ });\
	    name=$${temp[1]};\
	    $(GCTA) --grm-bin $(HOME)/simulated/WTCCC17000.final.chr$(CHR).gcta --reml-no-constrain --pheno $(HOME)/simulated/set$(SET)/WTCCC17000.final.chr$(CHR).qt.pheno --mpheno $$count --reml --out $(HOME)/simulated/set$(SET)/QT/WTCCC17000.final.chr$(CHR).gcta --thread-num $(THREADS) > /dev/null 2>&1;\
	    gctaRes=($$(grep "V(G)/Vp" $(HOME)/simulated/set$(SET)/QT/WTCCC17000.final.chr$(CHR).gcta.hsq));\
	    echo "$$p $${causal[$$cause]} $${herit[$$h]} $${gctaRes[1]} $${gctaRes[2]}" >> $(HOME)/simulated/set$(SET)/QT/gcta.res;\
	    ((heritIndex+=1));\
	    if [[ $$heritIndex -eq 10 ]]; then \
		heritIndex=0; \
		((causeIndex+=1));\
		if [[ $$cause -eq 5 ]]; then \
		    causeIndex=0;\
		    ((p+=1));\
		fi;\
	    fi;\
	    ((count+=1));\
	done;\
	perm=1;\
	count=1;\
	causeIndex=0;\
	heritIndex=0;\
	echo "Binary Trait Estimation (GCTA)";\
	rm -f $(HOME)/simulated/set$(SET)/CC/gcta.res;\
	echo "Perm Causal Prevalence Herit Estimate SE" > $(HOME)/simulated/set$(SET)/CC/gcta.res;\
	for i in `ls $(HOME)/simulated/set$(SET)/CC/*.assoc`; do \
	    temp=($${i//./ });\
	    name=$${temp[2]};\
	    temp=($${i//_/ });\
	    prev=$${temp[2]};\
	    $(GCTA) --grm-bin $(HOME)/simulated/WTCCC17000.final.chr$(CHR).gcta --reml-no-constrain --pheno $(HOME)/simulated/set$(SET)/WTCCC17000.final.chr$(CHR)_$${observePrev}_${prev}.cc.pheno --mpheno $$count --reml --out $(HOME)/simulated/set$(SET)/CC/WTCCC17000.final.chr$(CHR).gcta --thread-num $(THREADS) --prevalence $$prev > /dev/null 2>&1;\
	    gctaRes=($$(grep "V(G)/Vp_L" $(HOME)/simulated/set$(SET)/CC/WTCCC17000.final.chr$(CHR).gcta.hsq));\
	    echo "$$p $${causal[$$cause]} $$prev $${herit[$$h]} $${gctaRes[1]} $${gctaRes[2]}" >> $(HOME)/simulated/set$(SET)/CC/gcta.res;\
	    ((heritIndex+=1));\
	    if [[ $$heritIndex -eq 10 ]]; then \
		heritIndex=0; \
		((causeIndex+=1));\
		if [[ $$cause -eq 5 ]]; then \
		    causeIndex=0;\
		    ((p+=1));\
		fi;\
	    fi;\
	    ((count+=1));\
	done;\
	perm=1;\
	count=1;\
	causeIndex=0;\
	heritIndex=0;\
	echo "Extreme Selection Estimation (GCTA)";\
	rm -f $(HOME)/simulated/set$(SET)/EX/gcta.res;\
	echo "Perm Ratio Prevalence Herit Estimate SE" > $(HOME)/simulated/set$(SET)/EX/gcta.res;\
	for i in `ls $(HOME)/simulated/set$(SET)/EX/*.qassoc`; do \
	    temp=($${i//./ });\
	    name=$${temp[1]};\
	    temp=($${i//_/ });\
	    ratio=$${temp[1]};\
	    $(GCTA) --grm-bin $(HOME)/simulated/WTCCC17000.final.chr$(CHR).gcta --reml-no-constrain --pheno $(HOME)/simulated/set$(SET)/WTCCC17000.final.chr$(CHR)_$${ratio}.ex.pheno --mpheno $$count --reml --out $(HOME)/simulated/set$(SET)/CC/WTCCC17000.final.chr$(CHR).gcta --thread-num $(THREADS) > /dev/null 2>&1;\
	    gctaRes=($$(grep "V(G)/Vp" $(HOME)/simulated/set$(SET)/EX/WTCCC17000.final.chr$(CHR).gcta.hsq));\
	    echo "$$p $${causal[$$cause]} $$ratio $${herit[$$h]} $${gctaRes[1]} $${gctaRes[2]}" >> $(HOME)/simulated/set$(SET)/EX/gcta.res;\
	    ((heritIndex+=1));\
	    if [[ $$heritIndex -eq 10 ]]; then \
		heritIndex=0; \
		((causeIndex+=1));\
		if [[ $$cause -eq 5 ]]; then \
		    causeIndex=0;\
		    ((p+=1));\
		fi;\
	    fi;\
	    ((count+=1));\
	done;\	
	
cleanLog: ## Delete all the qsub records
	@rm -f *.e2*;
	@rm -f *.o2*;
	@rm -rf .lock;
	@rm -rf .set*;

cleanSim: cleanLog ## Delete all previous simulation
	@rm -rf simulated;
	@rm -rf .lock;
	@rm -rf .set*;
	@echo " " | crontab -;
QCLog: ## Print the default parameters for the whole procedure
	@echo "Parmeters for the preprocessing are:";
	@echo "Relationship threshold = $(RELATED)";
	@echo "GENO threshold = $(GENO)";
	@echo "MIND threshold = $(MIND)";
	@echo "MAF threshold = $(MAF)";
	@echo "HWE threshold = $(HWE)";
	@echo "Number of thread used = $(THREADS)";

clear: cleanSim cleanLog ## Delete everything and start over
	@rm -rf samples
	@rm -rf simulated
help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) |  awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
