# VCall_Script_For_Variant_Calling
## Author: ItsAquaDeath
This project aimed to design a fully functional **bash** script to perform all basic steps on a variant calling pipeline starting from knowing the Sequence Read Archive (SRA) ID of the files. 
The script can be executed in bash so long as the following programs are installed and added to the $PATH variable in your distribution:

- [SRA-Toolkit](https://github.com/ncbi/sra-tools)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [cutadapt](https://cutadapt.readthedocs.io/en/stable/)
- [Burrows-Wheeler Aligner (BWA)](https://bio-bwa.sourceforge.net/)
- [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [Samtools and BCFtools](https://www.htslib.org/)
- [Freebayes](https://github.com/freebayes/freebayes)
- [Genome Analysis Toolkit (GATK)](https://gatk.broadinstitute.org/hc/en-us)
- [snpEff](https://pcingola.github.io/SnpEff/#snpeff)

I don't own any rights regarding the forementioned programs, and all credit goes to their respective authors.

The shell script is a result of an individual project within the MsC in Bioinformatics from Universitat de València (UV), Spain. With this last sentence you might be thinking: "Woah, a script made by a newbie bioinformatician... It's going to be a total mess! I can surely do something much better than this!" And you'd be totally right about that! I'm still a bit of a noob with shell scripting, so there are surely a ton of areas within the script that can be further optimized or better programmed.

I'm currently not intending to upgrade this version of the script. Maybe in the future? Who knows...

Anyway, if this helped you in your personal project I'm more than happy! :D
