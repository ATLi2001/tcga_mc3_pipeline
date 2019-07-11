#!/bin/bash
set -e
JOBS=$1
if [ ! -e $1 ];then
	echo Job file not found!
	exit 1
fi

#/home/austintianli/tcga_phylowgs
VCFDIR=$2

# || true to skip errors
for SAMPLE in $(cat $JOBS); do
	echo ${SAMPLE}
    if [ -e phylowgs_input/${SAMPLE}/cnv_data.txt ];then
        echo Skip ${SAMPLE}
        continue
    fi
	mkdir -p phylowgs_input/${SAMPLE}
    gsutil cp gs://bucket1q2w/phylowgs_mc3/non_pan12/titancna_res/${SAMPLE}/optimalClusterSolution.txt phylowgs_input/${SAMPLE}/ || true
	cellularity=$(tail -n +2 phylowgs_input/${SAMPLE}/optimalClusterSolution.txt | cut -f 6) || true
	echo "Cellularity = ${cellularity}"
	optimal=$(tail -n +2 phylowgs_input/${SAMPLE}/optimalClusterSolution.txt | cut -f 11 | sed "s/\/\//\//" | sed "s/results\/titan\/hmm\///").segs.txt || true
	echo "Optimal Solution: ${optimal}"

	gsutil cp gs://bucket1q2w/phylowgs_mc3/non_pan12/titancna_res/${optimal} phylowgs_input/${SAMPLE}/ || true

	solution=phylowgs_input/${SAMPLE}/$(basename ${optimal}) || true
	python /home/austintianli/phylowgs_titancna/parser/parse_cnvs.py -f titan -c ${cellularity} \
	--cnv-output phylowgs_input/${SAMPLE}/cnv.txt ${solution} || true

	python /home/austintianli/phylowgs_titancna/parser/create_phylowgs_inputs.py --cnv sample1=phylowgs_input/${SAMPLE}/cnv.txt \
	--vcf-type sample1=mutect_smchet sample1=${VCFDIR}/${SAMPLE:0:12}.filtered.vcf.gz --output-cnvs phylowgs_input/${SAMPLE}/cnv_data.txt --output-variants phylowgs_input/${SAMPLE}/ssm_data.txt || true
	echo $(wc -l phylowgs_input/${SAMPLE}/cnv_data.txt) 
done
