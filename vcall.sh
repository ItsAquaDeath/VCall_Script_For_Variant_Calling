#!/bin/bash

# Usage page of the program
help(){
cat << EOF

Usage: bash $0 -g [genome] SRR_Acc_List.txt

    This script allows bioinformatic analysis of genomic variants from fastq
    files available from the Sequence Read Archive (SRA) inputing a properly
    formatted SRR_Acc_List.txt file (obtained from SRA Run Selector)

    ATTENTION: It is NECESSARY to locate most of the programs
               in a directory within $\PATH variable. This includes
               fastqc, cutadapt, bwa, samtools, bcftools, freebayes,
               and snpEff.

               For picard and gatk, add the path to the program in the
               first lines in the script.

       mandatory arguments:
       -g [reference genome in fasta format]
       -b [Technology target BED file (only mandatory for freebayes protocol)]
       -p [Ploidy file (only mandatory with samtools/bcftools protocol)]

       optional arguments:
       -o [outdir]
       -I [input dir with fastq files]
       -m [aligner: bwa (default) / bowtie2]
       -r [remove duplicates: yes (default) / no]
       -t [type: freebayes (default) / bcftools / gatk]
       -f [Fisher Strand, a measure of strand bias of the allele in GATK protocol, 60 (default)]
          [In bctools and freebayes protocol works as fraction of alternative/reference, 0.3 (default)]
       -q [quality by allele depth in GATK protocol,  2 (default)]
          [In bcftools and freebayes protocol works as read quality threshold, 20 (default)]
       -c [perform quality check with fastqc and trimming with cutadapt: no (default) / yes]
       -d [ReadsPosRankSum Test value in GATK protocol, -8 (default)]
          [In bcftools and freebayes protocol works as depth threshold, 10 (default)]
       -a [Name of the annotation database. GRCh38.p13 (human) is default]
       -e [Error solving. This solves "the sample list cannot be null or empty" error: no (default) / yes]

EOF
exit
}

#######################
#### PROGRAM PATHS ####
#######################
picard="java -jar /opt/picard.jar"
gatk=""

###########################
#### ARGUMENT PARSING #####
###########################

# Fetch SRR_Acc_List.txt filename
filelist=$(echo "$@" | tr " " "\n" | grep -E "SRR(.)*.txt")


while getopts "m:b:g:t:r:d:q:c:o:I:f:p:a:e:h" option;
do
	case $option in
		m) mapper=$OPTARG;;
		b) bed=$OPTARG;;
		g) refgenome=$OPTARG;;
		t) type=$OPTARG;;
		r) rmdup=$OPTARG;;
		d) rprs=$OPTARG;;
		q) qd=$OPTARG;;
		c) qcheck=$OPTARG;;
		o) outdir=$OPTARG;;
		I) indir=$OPTARG;;
		f) fs=$OPTARG;;
		p) ploidy=$OPTARG;;
		a) annotation=$OPTARG;;
		e) err=$OPTARG;;
		h) help;;
		?) help;;
	esac
done

# Some default values if unfilled
if [[ -z $type ]];
then
	type="freebayes"
fi
if [[ -z $annotation ]];
then
	annotation="GRCh38.p13"
fi

#############################
#### MANDATORY ARGUMENTS ####
#############################
if [[ -z $refgenome ]];
then
	help
fi
if [[ $type == "bcftools" ]] && [[ -z $ploidy ]];
then
	echo ""
	echo "Argument -p [ploidy] is mandatory for samtools/bcftools protocol"
	help
fi
if [[ $type == "freebayes" ]] && [[ -z $bed ]];
then
	echo ""
	echo "Argument -b [BED] is mandatory for freebayes protocol"
	help
fi

###################################################
#### DEFAULT VALUES IF ARGUMENT IS NOT DEFINED ####
###################################################
if [[ -z $qd ]] && [[ $type == "gatk" ]];
then
	qd=2.0
elif [[ -z $qd ]] && [[ $type != "gatk" ]];
then
	qd=20
fi

if [[ -z $rprs ]] && [[ $type == "gatk" ]];
then
	rprs=-8.0
elif [[ -z $rprs ]] && [[ $type != "gatk" ]];
then
	rprs=10
fi

if [[ -z $fs ]] && [[ $type == "gatk" ]];
then
	fs=60.0
elif [[ -z $fs ]] && [[ $type != "gatk" ]];
then
	fs=0.3
fi

#########################
#### INITIAL MESSAGE ####
#########################
echo "*****************************************************************"
echo "*******************   Vcall initialized!  ***********************"
echo "*****************************************************************"
echo "*** - Author: Álvaro Martínez Martínez ***"
echo "*** - Call: vcall $* ***"
echo "*****************************************************************"

# Create output dir if necessary
if [[ -z $outdir ]];
then
	outdir="output"
	mkdir -p $outdir
	echo "$outdir/ has been created to save output!"
else
	outdir=$(basename $outdir)
	echo "$outdir/ will be taken to dump program output files!"
fi


###############################
#### DOWNLOAD FILES CHECK #####
###############################
if [[ -z $indir ]];
then
	# Download from SRA
	prefetch -O $outdir -p --option-file $filelist
	echo ""
	echo "*****************************************************************"
	echo Download of files has finished!
	echo Proceeding with conversion from sra to fastq...
	# Order newly downloaded files
	ls $outdir | while read -r dir; 
	do 
		mv ./$outdir/$dir/*.sra ./$outdir;
		rmdir ./$outdir/$dir;
	done
	# Perform fasterq-dump
	ls $outdir | while read -r file;
	do
		fasterq-dump -e 24 $outdir/$file -O $outdir;
		rm $outdir/$file;
	done
	echo ""
	echo "*****************************************************************"
	echo Files have been successfully converted to fastq!
	echo "*****************************************************************"
else
	# If fastq are allocated locally, list them and check they are
	# compatible with SRR.text file
	indir=$(basename $indir)
	ls $indir | while read -r file;
	do
		basename $file | cut -d "." -f1 | cut -d "_" -f1 >> check.temp
	done
	# Remove duplicated lines due to possible _1 _2 fastq files
	awk -i inplace '!seen[$0]++' check.temp
	# Scan differences with diff
	difference=$(diff $filelist check.temp)
	rm check.temp
	if [[ -z $difference ]];
	then
		echo All files in file list are already downloaded in $indir/!
		echo Script will proceed with local allocated files...
	else
		echo $filelist file indicates not all fastq files are downloaded in $indir/!
		echo Please, make sure this file matches with EVERY fastq file located in $indir/ before intiating the script!
		exit 1
	fi
fi

####################################
#### QUALITY CHECK PRE-TRIMMING ####
####################################
if [[ $qcheck == "yes" ]];
then
	echo "*****************************************************************"
	echo "Checking quality of sequencing data..."
	mkdir $outdir/quality_pre
	cat $filelist | while read -r file;
	do
		if [[ -z $indir ]];
		then
			fastqc -t 24 -o $outdir/quality_pre $outdir/$file.fastq
		else
			fastqc -t 24 -o $outdir/quality_pre $indir/$file.fastq
		fi
	done
	echo "Quality reports pre-trimming saved at $outdir/quality_pre/"
	echo "*****************************************************************"
fi


######################
#### TRIM LIBRARY ####
######################
if [[ $qcheck == "yes" ]];
then
	cat $filelist | while read -r file;
	do
		echo "*****************************************************************"
		echo "Performing trimming of files $file..."
		if [[ -z $indir ]];
		then
			cutadapt --max-n 0.1 -q 10 -j 24 -o $outdir/${file}_trimmed.fastq $outdir/$file.fastq
			mv $outdir/${file}_trimmed.fastq $outdir/$file.fastq
		else
			cutadapt --max-n 0.1 -q 10 -j 24 -o $indir/${file}_trimmed.fastq $indir/$file.fastq
			mv $indir/${file}_trimmed.fastq $outdir/$file.fastq
		fi
		echo "$file trimmed and saved at $outdir/$file.fastq!"
	done
fi


#####################################
#### QUALITY CHECK POST-TRIMMING ####
#####################################
if [[ $qcheck == "yes" ]];
then
	echo "*****************************************************************"
	echo "Checking quality of sequencing data..."
	mkdir $outdir/quality_post
	cat $filelist | while read -r file;
	do
			fastqc -t 24 -o $outdir/quality_post $outdir/$file.fastq
	done
	echo "Quality reports post-trimming saved at $outdir/quality_post/"
	echo "*****************************************************************"
fi


##########################
#### MAPPER SELECTION ####
##########################
echo "*****************************************************************"
echo "Performing alignments with reference genome..."
while read -r fastq;
do
	echo "*****************************************************************"
	echo "Alignment started with $fastq.fastq file..."
	if [[ -z $mapper ]] || [[ $mapper == "bwa" ]];
	then
		if [[ -z $indir ]];
		then
			bwa mem -M -t 24 -o $outdir/$fastq.sam $refgenome $outdir/$fastq.fastq
		else
			bwa mem -M -t 24 -o $outdir/$fastq.sam $refgenome $indir/$fastq.fastq
		fi
	elif [[ $mapper == "bowtie" ]]; 
	then
		if [[ -z $indir ]];
		then
			bowtie2 -q -p 24 -x $refgenome -1 $outdir/$fastq.fastq -S $outdir/$fastq
		else
			bowtie2 -q -p 24 -x $refgenome -1 $indir/$fastq.fastq -S $outdir/$fastq
		fi
	fi
	# Check aligment result
	if [[ $? -eq 0 ]];
	then
		echo "Alignment with file $fastq.fastq: SUCCESS!"
		echo "Now converting into bam format..."
		samtools view -bSh $outdir/$fastq.sam -o $outdir/$fastq.bam
		rm $outdir/$fastq.sam
	else
		echo "Alignment with file $fastq.fastq: FAIL..."
		exit 1
	fi
	# Sort the bam file and remove the unsorted to save disk space
	echo "Sorting $fastq.bam file..."
	samtools sort --threads 24 $outdir/$fastq.bam -o $outdir/$fastq.sorted.bam
	rm $outdir/$fastq.bam
	echo "Mapping done with $fastq.fastq!"
	echo "*****************************************************************"	
done < $filelist
echo "Aligment phase has concluded with no issues!"
echo "*****************************************************************"


#########################
#### MARK DUPLICATES ####
#########################
echo "*****************************************************************"
echo "Marking duplicated reads for all bam files"
while read -r name;
do
	if [[ -z $rmdup ]] || [[ $rmdup == "yes" ]];
	echo "Marking duplicates for file $name.sorted.bam..."
	then
		if [[ $type != "bcftools" ]];
		then
			$picard MarkDuplicates INPUT=$outdir/$name.sorted.bam OUTPUT=$outdir/$name.sorted.rmDup.bam METRICS_FILE=$outdir/$name.sorted.rmDup.metrics REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=tmp/ 2>/dev/null
		else
			samtools rmdup $outdir/$name.sorted.bam $outdir/$name.sorted.rmDup.bam 2>/dev/null
		fi
		if [[ $? -eq 0 ]];then
			echo "Mark duplicates: SUCCESS!"
		else
			echo "Mark duplicates: FAIL..."
			exit 1
		fi
	else
		mv $outdir/$name.sorted.bam $outdir/$name.sorted.rmDup.bam
	fi
	samtools index $outdir/$name.sorted.rmDup.bam
done < $filelist
echo "All files have been indexed for variant calling!"
echo "*****************************************************************"

#######################################
#### SOLVE POSSIBLE ERROR WITH @RG ####
#######################################
# If there are no @RG defined in bam files an error pops with GATK
# This code solve that by adding @RG dummy groups to the files
if [[ $err == "yes" ]];
then
	echo "*****************************************************************"
	echo "Filling dummy @RG in bam files to avoid gatk error..."
	count=1
	cat $filelist | while read -r name;
	do
		$picard AddOrReplaceReadGroups I=$outdir/$name.sorted.rmDup.bam O=$outdir/$name.sorted.RG.bam RGLB=$name RGPL=illumina RGPU=unit1 RGSM=Sample$count CREATE_INDEX=True
		rm $outdir/$name.sorted.rmDup.bam
		rm $outdir/$name.sorted.rmDup.bam.bai
		mv $outdir/$name.sorted.RG.bam $outdir/$name.sorted.rmDup.bam
		mv $outdir/$name.sorted.RG.bai $outdir/$name.sorted.rmDup.bam.bai
		echo "Finished with file $name.sorted.bam!"
		count=$((count + 1))
	done
fi
echo "All bam files have been formated for downstream analysis!"
echo "*****************************************************************"

########################
#### VARIANT CALLER ####
########################
while read -r name;
do
	echo "*****************************************************************"
	echo "Starting Variant Calling of file $name.sorted.rmDup.bam..."
	if [[ $type == "freebayes" ]];
	then
		freebayes --min-base-quality 3 -F 0.15 -f $refgenome -t $bed -b $outdir/$name.sorted.rmDup.bam -v $outdir/$name.vcf
		# Check last command status
		if [[ $? -eq 0 ]];then
			echo "Variant calling of $name: SUCCESS"
			echo "Variants located in $name.vcf.gz"
		else
			echo "Variant calling of $name: FAIL"
			exit 1
		fi
		bgzip $outdir/$name.vcf
	fi
	if [[ $type == "bcftools" ]];
	then
		bcftools mpileup -Ou -f $refgenome $outdir/$name.sorted.rmDup.bam 2>/dev/null | bcftools call -mv -Ob --ploidy-file $ploidy -o $outdir/$name.bcf 2>/dev/null
		# Check last command status
		if [[ $? -eq 0 ]];then
			echo "Variant calling of $name: SUCCESS"
			echo "Variants located in $name.vcf.gz"
		else
			echo "Variant calling of $name: FAIL"
			exit 1
		fi
		bcftools view $outdir/$name.bcf > $outdir/$name.vcf
		bgzip $outdir/$name.vcf
	fi
	if [[ $type == "gatk" ]];
	then
		# Create a dict for the reference genome first
		dict=$(echo $refgenome | cut -d'.' -f1)
		if [[ ! -f $dict.dict ]];
		then
			$picard CreateSequenceDictionary R=$refgenome O=$dict.dict 2>/dev/null
		fi
		samtools faidx $refgenome
		# Then apply gatk
		$gatk --java-options "-Xmx4g" HaplotypeCaller -R $refgenome -I $outdir/$name.sorted.rmDup.bam -O $outdir/$name.vcf.gz 2>/dev/null
		# Check last command status
		if [[ $? -eq 0 ]];then
			echo "Variant calling of $name: SUCCESS"
			echo "Variants located in $name.vcf.gz"
		else
			echo "Variant calling of $name: FAIL"
			exit 1
		fi
	fi
	# Elaborate an index for the resulting vcf file
	tabix -f -p vcf $outdir/$name.vcf.gz
done < $filelist
echo "Variant calling completed!"
echo "*****************************************************************"

##########################
#### FILTERS TO APPLY ####
##########################
# Filter are different depending on the variant calling protocol used
if [[ $type == "freebayes" ]];
then
	echo "*****************************************************************"
	echo "Filtering variants according to the filters of the script:"
	printf "\t - Read quality (Q): %s.\n" "$qd"
	printf "\t - Fraction of alternative alleles: %s.\n" "$fs"
	printf "\t - Coverage (ReadDepth): %s.\n" "$rprs"
	while read -r name;
	do
		echo bcftools filter -i "'QUAL>=$qd && (AB>=$fs || AB==0) && INFO/DP>=$rprs'" -O z -o $outdir/$name.filtered.vcf.gz $outdir/$name.vcf.gz | sh
		if [[ $? -eq 0 ]];
		then
			echo "Variant filter of file $name.vcf.gz: SUCCESS"
		else
			echo "Variant filter of file $name.vcf.gz: FAIL"
			exit 1
		fi
	done < $filelist
fi

if [[ $type == "bcftools" ]];
then
	echo "*****************************************************************"
	echo "Filtering variants according to the filters of the script:"
	printf "\t - Read quality (Q): %s.\n" "$qd"
	printf "\t - Fraction of alternative alleles: %s.\n" "$fs"
	printf "\t - Coverage (ReadDepth): %s.\n" "$rprs"
	while read -r name;
	do
		echo bcftools filter -i "'QUAL>=$qd && (AF>=$fs) && INFO/DP>=$rprs'" -O z -o $outdir/$name.filtered.vcf.gz $outdir/$name.vcf.gz | sh
		if [[ $? -eq 0 ]];
		then
			echo "Variant filter of file $name.vcf.gz: SUCCESS"
		else
			echo "Variant filter of file $name.vcf.gz: FAIL"
			exit 1
		fi
	done < $filelist
fi

if [[ $type == "gatk" ]];
then
	echo "*****************************************************************"
	echo "Filtering variants according to the filters of the script:"
	printf "\t - QualByDepth (QD): %s.\n" "$qd"
	printf "\t - FisherStrand (FS): %s.\n" "$fs"
	printf "\t - ReadPosRankSum (RPRS): %s.\n" "$rprs"
	while read -r name;
	do
		$gatk VariantFiltration -R $refgenome -V $outdir/$name.vcf.gz -O $outdir/$name.filtered.vcf.gz -filter "QD<$qd" --filter-name "QD" -filter "FS>$fs" --filter-name "FS" -filter "ReadPosRankSum<$rprs" --filter-name "ReadPosRankSum"
		if [[ $? -eq 0 ]];
		then
			echo "Variant filter of file $name.vcf.gz: SUCCESS"
		else
			echo "Variant filter of file $name.vcf.gz: FAIL"
			exit
		fi
	done < $filelist
fi
echo "Variant filtering completed!"
echo "*****************************************************************"

############################
#### VARIANT ANNOTATION ####
############################
echo "*****************************************************************"
echo "Now annotating variants with snpEff..."
while read -r name;
do
	snpEff -v $annotation -s $outdir/"$name"_snpEff.html $outdir/$name.filtered.vcf.gz > $outdir/$name.annotated.vcf
	grep -v "^GL" $outdir/$name.annotated.vcf > $outdir/$name.tmp; mv $outdir/$name.tmp $outdir/$name.annotated.vcf

	if [[ $? -eq 0 ]];then
		echo "Variant Annotation of $name.filtered.vcf.gz: SUCCESS"
		echo "*****************************************************************"
	else
		echo "Variant annotation of $name.filtered.vcf.gz: FAIL"
		exit 1
	fi
done < $filelist
echo "Variant annotation completed!"
echo "*****************************************************************"