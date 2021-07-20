#!/bin/bash -ue


## This script is meant to be run in a directory with a directory of fastq files, samples.txt, comparisons.txt, and a genome fasta file, and gtf file.
## example:
## Project
##  |--RNA.Seq.bash.2.sh
##  |--genome.fa
##  |--genes.gtf
##  |--sammples.txt
##  |--comparisons.txt
##  |--Raw_fastq
##	|--Sample1_fastq.gz
##	|--Sample2.fastq.gz		


## Run script using nohup. First set as executable, then run using nohup.
## chmod +x RNA-Seq.hisat.bash.sh
## screen -S RNAseq
## ./RNA-Seq.hisat.bash.sh > log.txt &


## variables to change:
## Definitely check that the samples, comparisons, genome, gtf fastq, se, strand and length variables have beed modified to suit your data.

samples=samples_test.txt
##The samples sheet is a tab-delimited text file with the following columns:
##id sampleName treatment fastq1 fastq2
## where fastq1 and fastq2 are the fastq file names, if data is single-end, only fastq1 column is required.
## when copying text files from windows and macs, check the file format of the sample sheet with vim first.
## make sure there are separate lines, and you don't see "[noeol]" at the bottom of the screen.
threads=10
length=150
##length is length of reads
genome=genome.fa
## create a symbolic link the annotation_db directory to your current directory.
## example for mouse genome.fa, (the space and period ath the end is part of the command):
## ln -s /lower_bay/local_storage/annotation_db/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa .
gtf=genes.gtf
## same for the genes.gtf
##ln -s /lower_bay/local_storage/annotation_db/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf .
deg=true
splicing=false
fastq=Raw_data
## The directory that you have your fastq files in		
comparisons=comparisons_test.txt
## tab delimited text file  with two columns:
## treat control

strand=""
##strandness:
## "" -unstranded
##rf -reversely stranded - use for illumina directional libraries, will say "Illumina Truseq stranded"
##fr stranded
se=false
##if data is paired end, set se to false
outDir=RNA-Seq
bigwigs=true
##select true to generate bigwigs
annotation=mm10
#####

#################
##testing if options worked
##########################################





if [ -z "${samples:-}" ]
then
	echo "please specify a samples file"
	exit
elif [ ! -e $samples ]
then
        echo "sample file $samples not found"
        exit 1

fi



if [ -z "${comparisons:-}" ] && [ "$deg" = "true" ]
then
	echo "please specify a comparisons file"
	exit 1

elif [ ! -e $comparisons ] && [ "$deg" =  "true" ]
then
        echo "comparisons file $comparisons not found"
        exit 1

fi



if [ -z "${genome:-}" ]
then
	echo "please specify a reference genome fasta file"
	exit 1

elif [ ! -e $genome ]
then
	echo "Genome fasta file not found"
	exit 1

fi


if [ -z "${gtf:-}" ]
then
	echo "please specify a gtf file"
elif [ ! -e $gtf ]
then

	echo "gtf file not found"
	exit 1

fi


if [ -z "${length:-}" ]
then
	echo "please specify the length of your reads"
	exit 1
fi


if [ -z "${threads:-}" ]
then

	threads=1
fi


if [ -z "${outDir:-}" ]
then
	now=$(date +"%m_%d_%Y")
	outDir=RNA-Seq_out_${now}
	echo "No output directory specified. Using output directory $outDir"

fi



if [ -z "${splicing:-}" ]
then

        splicing=false
        echo "Splicing Analysis will not be done"
fi


if [ -z "${deg:-}" ]
then

        deg=false
        echo "Differential gene expression analysis will not be done"
fi


if [ "$deg" != "true" ] && [ "$splicing" != "true" ]
then
        echo "Only alignment with STAR will be done"
fi



if [ -z "${strand:-}" ]
then
	echo "strandness not specified. setting strandness to unstranded"
	strand=""	
fi

if [ -z "${se:-}" ]
then

	echo "data is single-end: $se "
	se=false
fi 


if [ -z "${fastq:-}" ]
then
	echo "Please specify the location of your fastq files"
	exit 1
fi


#####################
##set up directories
#can be changed


fastqc=($outDir/FastQC)
trimmed_fastq=($outDir/Trimmed_fastq)
trimmed_fastqc=($outDir/Trimmed_fastQC)
star1=($outDir/STAR1pass)
star1index=($outDir/STAR1index)
cuffnorm=($outDir/Cuffnorm)
counts=($outDir/Counts)
assembly=($outDir/Assembly)
ballgown=($outDir/Ballgown)
deseq=($outDir/DESeq2)
bw_dir=($outDir/BigWigs)
multiqc=($outDir/Multiqc)
enrich=($outDir/Cluster_Analysis)

## add hisat
hisat=($outDir/Hisat2_alignments)

##############step1. FastQC


if [ ! -d $fastqc ]
then
	mkdir -p $fastqc
fi



fastqc -t 10 -o $fastqc $fastq/*gz


######## Step2. Remove adaptors


if [ ! -d $trimmed_fastq ]
then
	mkdir -p $trimmed_fastq
fi


while read -r coreNumber sampleName group r1 r2
do
	if [ $se != "true" ]
	then

	read1=($fastq/${r1})
	read2=($fastq/${r2})
	i=$coreNumber
	echo trimming adaptors from $i
	 java -jar /usr/share/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads $threads -phred33 -trimlog $trimmed_fastq/Trimmomatic_${i}.log \
         $read1 $read2 $trimmed_fastq/${i}_1_at.fq.gz $trimmed_fastq/${i}_1_unpaired.fq.gz $trimmed_fastq/${i}_2_at.fq.gz $trimmed_fastq/${i}_2_unpaired.fq.gz \
         ILLUMINACLIP:/usr/share/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10:8:TRUE LEADING:3 TRAILING:20 MINLEN:12 2> $trimmed_fastq/${i}_trimStats.txt
	
	gzip $trimmed_fastq/Trimmomatic_${i}.log
	else
	
	read1=($fastq/${r1})
        i=$coreNumber
        echo trimming adaptors from $i
         java -jar /usr/share/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads $threads -phred33 -trimlog $trimmed_fastq/Trimmomatic_${i}.log \
         $read1  $trimmed_fastq/${i}_1_at.fq.gz \
         ILLUMINACLIP:/usr/share/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10:8:TRUE LEADING:3 TRAILING:20 MINLEN:12 2> $trimmed_fastq/${i}_trimStats.txt
	
	gzip $trimmed_fastq/Trimmomatic_${i}.log
	fi

done < $samples




##### FastQC again



fastqc -t 10 -o $trimmed_fastq $trimmed_fastq/*at.fq.gz


if [ ! -d $hisat ]
then

	mkdir $hisat
fi



for i in $(ls $trimmed_fastq/*[12]_at.fq.gz | cut -f3 -d "/" | sed 's/_[12]_at.fq.gz//g' | sort | uniq)
do

	echo $i

        r1=$(ls $trimmed_fastq/${i}_1_at.fq.gz)
        r2=$(ls $trimmed_fastq/${i}_2_at.fq.gz)

	echo $r1
	echo $r2
        echo aligning $i using hisat

        hisat2 -p $threads --new-summary --secondary -x /lower_bay/local_storage/annotation_db/Mus_musculus/mm10_hisat2_index/genome -1 $r1 -2 $r2 2> $hisat/${i}_hisat.out.txt | samtools view -hb - > $hisat/${i}_hisat.bam

	samtools sort -o $hisat/${i}_sorted.bam $hisat/${i}_hisat.bam

	samtools index $hisat/${i}_sorted.bam

	rm $hisat/${i}_hisat.bam
done


if [ "$bigwigs" = "true" ]
then
	if [ ! -d $bw_dir ]
	then
		mkdir -p "$bw_dir"
	fi


	for file in $hisat/*_sorted.bam
	do

		sample=$(echo $file | cut -f3 -d "/" | cut -f1 -d "_")
		
		bamCoverage -b $file --normalizeUsing CPM -of bigwig -o $bw_dir/${sample}.bw
	done

fi



##get expression using cuffnorm

if [ ! -d $cuffnorm ]
then

	mkdir $cuffnorm
fi

lib=""



if [ "$strand" = "rf" ]
then
lib="--library-type fr-firststrand"

elif [ "$strand" = "fr" ]
then

lib="--library-type fr-secondstrand"

fi


files=$(ls $hisat/*sorted.bam)
cuffnorm -o $cuffnorm -p $threads $lib $gtf $files



if [ ! -d $counts ]
then

	mkdir -p $counts
fi

fc_strand=""
if [ "$strand" = "rf" ]
then
	fc_strand="-s 2"

elif [ "$strand" = "fr" ]
then

	fc_strand="-s 1"

fi


if [ "$se" != "true" ]
then            
        
	featureCounts -a $gtf $fc_strand -p -o $counts/featureCounts_gene_counts.txt $hisat/*sorted.bam
        
else            
	featureCounts -a $gtf $fc_strand -o $counts/featureCounts_gene_counts.txt $hisat/*sorted.bam
        
fi              


if [ "$deg" = "true" ]
then


                        
       ##DESeq2         
                        
       ##add a prepfeatureCounts
                        
       if [ ! -d $deseq ]
       then             
                        
       	mkdir -p $deseq 
       fi               
                        
	prepfeatureCounts.R --counts=${counts}/featureCounts_gene_counts.txt --out=${counts}
	                
                        
	DESeq2.R --counts=${counts}/featureCounts_for_DESeq2.csv --annotation=${counts}/gene_lengths.csv --species=${annotation} --fpkms=${cuffnorm}/genes.fpkm_table --samples=${samples} --comparisons=${comparisons} --out=${deseq}

	if [ ! -d $enrich ]
	then

		mkdir -p $enrich

	fi

        ClusterProfiler_no_wikipathways.R --DE=${deseq} --org=${annotation} --out=${enrich}
fi                      
                        
                        
if [ "$splicing" = "true" ]
then                    
                        
                        
	if [ ! -d $assembly ]
	then            
		mkdir $assembly
                        
	fi              
	                
                        
	if [ $se != "true" ]
	then


	for i in $hisat/*sorted.bam #change
	do
		sample=$(basename $i _sorted.bam) #change
		
		stringtie -p 14 -G $gtf -o $assembly/${sample}.gtf $i
	done

	else
	for i in $hisat/*sorted.bam
        do
                sample=$(basename $i sorted.bam) #changed
                
                stringtie -p 14 -G $gtf -o $assembly/${sample}.gtf $i
        done
	fi

	

	ls $assembly/*gtf | tr " " "\n" > $assembly/gtf_files.txt

	gtf_filenames=$assembly/gtf_files.txt

	stringtie --merge -p 14 -G $gtf -o $assembly/merged.gtf $gtf_filenames

	if [ $se != "true" ]
	then


	for i in $hisat/*sorted.bam
	do
		sample=$(basename $i _sorted.bam)
		
		
		stringtie_dir=""
		if [ $strand = "fr" ]
		then

			strintie_dir="--fr"
		else
			stringtie_dir="--rf"

		fi
	if [ ! -d $ballgown ]
	then
		mkdir -p $ballgown
	fi
	
		stringtie -p $threads -e -B -o $ballgown/${sample}/${sample}.gtf -G $gtf $i

	done
	
	else

	for i in $hisat/*sorted.bam
        do
                sample=$(basename $i sorted.bam)
                

		if [ $strand = "fr" ]
                then

                        strintie_dir="--fr"
                else
                        stringtie_dir "--rf"

                fi
                
                stringtie -p $threads -e -B -o ballgown/${sample}/${sample}.gtf -G $gtf $i

        done

	fi
	python /usr/share/pathing/prepDE.py3 -i $ballgown -g $counts/stringtie_gene_count_matrix.csv -t $counts/transcript_count_matrix.csv
fi




#####################
if [ ! -d $multiqc ]
then
	mkdir $multiqc
fi

cd $multiqc
ln -s ../../$trimmed_fastq/*trimStats.txt .
ln -s ../../$star1/*final.out .
ln -s ../../$fastqc/*  .
ln -s ../../$trimmed_fastq/*fastqc* .
ln -s ../../$assembly/*tsv .
ln -s ../../$counts/*summary .
ln -s ../../$deseq/*mq* .
ln -s ../../$deseq/DE_summary.txt


multiqc .
cd -

#get software versions
fastqc --version &> $outDir/softwareVersions.txt
which trimmomatic-0.39.jar | cut -f5 -d "/" | cut -f1,2 -d "." &>> $outDir/softwareVersions.txt
STAR --version &>> $outDir/softwareVersions.txt
samtools --version | grep "samtools" &>> $outDir/softwareVersions.txt
echo stringtie $(stringtie --version) &>> $outDir/softwareVersions.txt
gffcompare --version &>> $outDir/softwareVersions.txt
