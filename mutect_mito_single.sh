#!/usr/bin/bash
#$ -cwd
#$ -v PYENV_ROOT
#$ -v PATH
#$ -v HTTP_PROXY
#$ -v LLVM_CONFIG
#$ -q r.q
#$ -pe smp 12

### mutect_mito_single.sh
### qsub script for mutect + EAGLE validation

# $1 to $3 are as specified from companion R script as inputs
F_TUMOR="$1" 
F_OUTPUT="$2"
F_SEQUENCE="$3"
F_EAGLE="$4"

echo "#   BAM INPUT : $F_TUMOR" >&2
echo "#   VCF OUTPUT: $F_OUTPUT" >&2
echo "#   REF GENOME: $F_SEQUENCE" >&2

### BAM indexing as necessary
echo "### QSUB START: $(date)" >&2

if [[ -f "$F_TUMOR" ]]
then
        FT_BAI=`echo $F_TUMOR | sed 's#\.bam#\.bai#'`
        if [[ ! -f "$F_TUMOR.bai" ]] || [[ ! -f "$FT_BAI" ]]
        then
                echo "## QSUB MSG: INDEXING TUMOR BAM" >&2
                /opt/bin/samtools index $F_TUMOR
        fi
        /opt/gatk/gatk GetSampleName -I $F_TUMOR -O $F_TUMOR.smid
        N_TUMOR=`cat $F_TUMOR.smid`
else
        echo "### QSUB ERR: TUMOR SAMPLE BAM (F_TUMOR) DOES NOT EXIST" >&2
        exit 102
fi

RS0=0
if [[ -f "$F_SEQUENCE" ]]
then
	F_SEQ_DICT=`echo $F_SEQUENCE | sed 's#\.fa*$#\.dict#'`
	if [[ ! -f "$F_SEQ_DICT" ]]
	then
		if [[ "$(dirname $F_SEQUENCE)" != "$(pwd)" ]]
		then
			RS0=26
			echo "## QSUB MSG: COPYING REFERENCE" >&2
			cp $F_SEQUENCE $(pwd)/
			F_SEQUENCE=`basename $F_SEQUENCE`
		fi
		echo "## QSUB MSG: CREATING SEQUENCE DICTIONARY" >&2
		/opt/gatk/gatk --java-options "-Xmx16G" CreateSequenceDictionary -R $F_SEQUENCE
	fi
else
        echo "### QSUB ERR: REFERENCE GENOME SEQUENCE FILE (F_SEQUENCE) DOES NOT EXIST" >&2
        exit 103
fi

echo "### QSUB MSG: BEGINNING MUTECT CALLS" >&2

F_HEADER_TUMOR=`basename $(echo $F_TUMOR | sed 's#\.bam.*$##')`

# Mutect2 mitochondrial mode
/opt/gatk/gatk --java-options "-Xmx16G" Mutect2 -R $F_SEQUENCE -I $F_TUMOR --output $F_OUTPUT --mitochondria-mode --disable-read-filter NotDuplicateReadFilter

# EAGLE
echo "### QSUB MSG: BEGINNING EAGLE FOR VALIDATION CALLS" >&2
echo "### QSUB MSG: OUTPUT: $(pwd)/$F_HEADER_TUMOR.egl.txt" >&2
/opt/bin/eagle -t 12 -n 0 --hetbias=0 -o $F_EAGLE  -v $F_OUTPUT -a $F_TUMOR  -r $F_SEQUENCE

if [[ "$RS0" -eq "26" ]]
then
	echo "## QSUB MSG: REMOVING COPIED REFERENCE"
	rm $F_SEQUENCE
	rm -rf $F_NORMAL.smid
fi

echo "### QSUB FINISHED: $(date)" >&2