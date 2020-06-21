configfile: "/home/x_seoju/ScRNACloneEvaluation/Tirosh/config.yaml"
configfile: "/home/x_seoju/ScRNACloneEvaluation/Tirosh/samples.yaml"

include: "ichorCNA.snakefile"
include: "getAlleleCounts.snakefile"
import os.path

CLUST = {1:[1], 2:[1,2], 3:[1,2,3], 4:[1,2,3,4], 5:[1,2,3,4,5], 6:[1,2,3,4,5,6], 7:[1,2,3,4,5,6,7], 8:[1,2,3,4,5,6,7,8], 9:[1,2,3,4,5,6,7,8,9], 10:[1,2,3,4,5,6,7,8,9,10]}
PLOIDY = {2:[2], 3:[2,3], 4:[2,3,4]}


rule all:
	input: 
		expand("/proj/sc_ml/Tirosh/bulk-wes/CY79/cna/results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.txt", tumor=config["pairings"], clustNum=CLUST[config["TitanCNA_maxNumClonalClusters"]], ploidy=PLOIDY[config["TitanCNA_maxPloidy"]]),
		expand("/proj/sc_ml/Tirosh/bulk-wes/CY79/cna/results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.ichor.seg.txt", tumor=config["pairings"], clustNum=CLUST[config["TitanCNA_maxNumClonalClusters"]], ploidy=PLOIDY[config["TitanCNA_maxPloidy"]]),
		expand("/proj/sc_ml/Tirosh/bulk-wes/CY79/cna/results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.ichor.cna.txt", tumor=config["pairings"], clustNum=CLUST[config["TitanCNA_maxNumClonalClusters"]], ploidy=PLOIDY[config["TitanCNA_maxPloidy"]]),
		"/proj/sc_ml/Tirosh/bulk-wes/CY79/cna/results/titan/hmm/optimalClusterSolution.txt",
		"/proj/sc_ml/Tirosh/bulk-wes/CY79/cna/results/titan/hmm/optimalClusterSolution/"
		
rule runTitanCNA:
	input:
		alleleCounts="/proj/sc_ml/Tirosh/bulk-wes/CY79/cna/results/titan/tumCounts/{tumor}.tumCounts.txt",
		corrDepth="/proj/sc_ml/Tirosh/bulk-wes/CY79/cna/results/ichorCNA/{tumor}/{tumor}.correctedDepth.txt"		
	output:		
		titan="/proj/sc_ml/Tirosh/bulk-wes/CY79/cna/results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.txt",
		param="/proj/sc_ml/Tirosh/bulk-wes/CY79/cna/results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.params.txt",
		segTxt="/proj/sc_ml/Tirosh/bulk-wes/CY79/cna/results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.segs.txt",
		seg="/proj/sc_ml/Tirosh/bulk-wes/CY79/cna/results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.seg"
	params:
		outRoot="/proj/sc_ml/Tirosh/bulk-wes/CY79/cna/results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}/",
		titanRscript=config["TitanCNA_rscript"],
		libdir=config["TitanCNA_libdir"],
		numCores=config["TitanCNA_numCores"],
		normal=config["TitanCNA_normalInit"],
		chrs=config["TitanCNA_chrs"],
		sex=config["sex"],
		genomeStyle=config["genomeStyle"],
		genomeBuild=config["genomeBuild"],
		cytobandFile=config["cytobandFile"],
		estimatePloidy=config["TitanCNA_estimatePloidy"],
		estimateClonality=config["TitanCNA_estimateClonality"],
		estimateNormal=config["TitanCNA_estimateNormal"],
		centromere=config["centromere"],
		alphaK=config["TitanCNA_alphaK"],
		#alphaR=config["TitanCNA_alphaR"],
		#alleleModel=config["TitanCNA_alleleModel"],
		txnExpLen=config["TitanCNA_txnExpLen"],
		plotYlim=config["TitanCNA_plotYlim"]
	log:
		"/proj/sc_ml/Tirosh/bulk-wes/CY79/cna/logs/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.log"
	shell:
		"Rscript {params.titanRscript} --hetFile {input.alleleCounts} --cnFile {input.corrDepth} --outFile {output.titan} --outSeg {output.segTxt} --outParam {output.param} --outIGV {output.seg} --outPlotDir {params.outRoot} --libdir {params.libdir} --id {wildcards.tumor} --numClusters {wildcards.clustNum} --numCores {params.numCores} --normal_0 {params.normal} --ploidy_0 {wildcards.ploidy} --genomeStyle {params.genomeStyle} --genomeBuild {params.genomeBuild} --cytobandFile {params.cytobandFile} --chrs \"{params.chrs}\" --gender {params.sex} --estimateNormal {params.estimateNormal} --estimatePloidy {params.estimatePloidy} --estimateClonality {params.estimateClonality}  --centromere {params.centromere} --alphaK {params.alphaK} --txnExpLen {params.txnExpLen} --plotYlim \"{params.plotYlim}\" > {log} 2> {log}"
	
rule combineTitanAndIchorCNA:
	input:
		titanSeg="/proj/sc_ml/Tirosh/bulk-wes/CY79/cna/results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.segs.txt", 
		titanBin="/proj/sc_ml/Tirosh/bulk-wes/CY79/cna/results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.txt",
		titanParam="/proj/sc_ml/Tirosh/bulk-wes/CY79/cna/results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.params.txt",
		ichorSeg="/proj/sc_ml/Tirosh/bulk-wes/CY79/cna/results/ichorCNA/{tumor}/{tumor}.seg.txt",
		ichorBin="/proj/sc_ml/Tirosh/bulk-wes/CY79/cna/results/ichorCNA/{tumor}/{tumor}.cna.seg",
		ichorParam="/proj/sc_ml/Tirosh/bulk-wes/CY79/cna/results/ichorCNA/{tumor}/{tumor}.params.txt"
	output:
		segFile="/proj/sc_ml/Tirosh/bulk-wes/CY79/cna/results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.ichor.seg.txt",
		binFile="/proj/sc_ml/Tirosh/bulk-wes/CY79/cna/results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.ichor.cna.txt",
	params:
		combineScript=config["TitanCNA_combineTitanIchorCNA"],
		libdir=config["TitanCNA_libdir"],
		centromere=config["centromere"],
		sex=config["sex"]
	log:
		"/proj/sc_ml/Tirosh/bulk-wes/CY79/cna/logs/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.combineTitanIchorCNA.log"
	shell:
		"Rscript {params.combineScript} --libdir {params.libdir} --titanSeg {input.titanSeg} --titanBin {input.titanBin} --titanParam {input.titanParam} --ichorSeg {input.ichorSeg} --ichorBin {input.ichorBin} --ichorParam {input.ichorParam} --sex {params.sex} --outSegFile {output.segFile} --outBinFile {output.binFile} --centromere {params.centromere} > {log} 2> {log}"	
	
rule selectSolution:
	input:
		#ploidyDirs=expand("/proj/sc_ml/Tirosh/bulk-wes/CY79/cna/results/titan/hmm/titanCNA_ploidy{ploidy}/", ploidy=PLOIDY[config["TitanCNA_maxPloidy"]]),
		resultFiles=expand("/proj/sc_ml/Tirosh/bulk-wes/CY79/cna/results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.txt", tumor=config["pairings"], clustNum=CLUST[config["TitanCNA_maxNumClonalClusters"]], ploidy=PLOIDY[config["TitanCNA_maxPloidy"]])
	output:
		"/proj/sc_ml/Tirosh/bulk-wes/CY79/cna/results/titan/hmm/optimalClusterSolution.txt"
	params:
		solutionRscript=config["TitanCNA_selectSolutionRscript"],
		threshold=config["TitanCNA_solutionThreshold"]
	log:
		"/proj/sc_ml/Tirosh/bulk-wes/CY79/cna/logs/titan/selectSolution.log"
	shell:
		"""
		ploidyRun2=/proj/sc_ml/Tirosh/bulk-wes/CY79/cna/results/titan/hmm/titanCNA_ploidy2/
		if [ -d /proj/sc_ml/Tirosh/bulk-wes/CY79/cna/results/titan/hmm/titanCNA_ploidy3/ ]; then
			ploidyRun3=/proj/sc_ml/Tirosh/bulk-wes/CY79/cna/results/titan/hmm/titanCNA_ploidy3/
		else
			ploidyRun3=NULL
		fi
		if [ -d /proj/sc_ml/Tirosh/bulk-wes/CY79/cna/results/titan/hmm/titanCNA_ploidy4/ ]; then
			ploidyRun4=/proj/sc_ml/Tirosh/bulk-wes/CY79/cna/results/titan/hmm/titanCNA_ploidy4/
		else
			ploidyRun4=NULL
		fi
		Rscript {params.solutionRscript} --ploidyRun2 $ploidyRun2 --ploidyRun3 $ploidyRun3 --ploidyRun4 $ploidyRun4 --threshold {params.threshold} --outFile {output} > {log} 2> {log}
		"""
		
rule copyOptSolution:
	input:
		"/proj/sc_ml/Tirosh/bulk-wes/CY79/cna/results/titan/hmm/optimalClusterSolution.txt"
	output:
		directory("/proj/sc_ml/Tirosh/bulk-wes/CY79/cna/results/titan/hmm/optimalClusterSolution/")
	params:
	log:
		"/proj/sc_ml/Tirosh/bulk-wes/CY79/cna/logs/titan/hmm/optSolution/copyOptSolution.log"
	shell:
		"""
		curDir=`pwd`
		for i in `cut -f11 {input} | grep -v "path"`;
		do
			echo -e "Copying ${{i}} to {output}"
			cp -r ${{i}}* {output}
		done		
		"""

	
