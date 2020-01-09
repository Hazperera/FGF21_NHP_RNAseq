# Variant calling and Bam Matcher pipeline documentation
# Sara Murray 
# November 2017


##### Installation #####

#official description and installation guide
# https://bitbucket.org/sacgf/bam-matcher



### Things to install for Bam Matcher
# python v2.7
    # python libraries:  (instructions for installing these in bam matcher installation wiki)
        #PyVCF
        #ConfigParser
        #Cheetah
        #pysam (requires python-dev and zlib1g-dev)
        #fisher (requires numpy and python-dev)

# varscan.v2.4.0 (java)        http://dkoboldt.github.io/varscan/
# samtools-1.4.1 (cmd line)    http://samtools.sourceforge.net/.
# bamMatcher (python)          https://bitbucket.org/sacgf/bam-matcher



### Things to install for variant calling (ie, to make and filter vcf file)
# picard (java)
# htslib-1.4.1 (cmd line)
# bcftools-1.4.1 (cmd line)




##### Make a config file for Bam Matcher #####
## chunk below is entire config file as written/commented, saved as Macaque.conf in Code/bam-matcher/


################################################


# BAM-matcher configuration file

# All parameters are in the format:
# KEYWORD:  VALUE

## DO NOT REMOVE SECTION HEADERS (e.g. [Variantcallers])
## DO NOT REMOVE PARAMETER KEYWORDS

# If not setting a specific parameter, just leave it blank, rather than deleting or commenting out the line
# Missing parameter keywords will generate errors

[VariantCallers]
# file paths to variant callers and other binaries

# This is the default caller to use (gatk, freebayes, or varscan)
# note:  if using varscan, only need to install varscan
caller:    varscan

# These are paths (or commands) to the caller executables
# full paths is always required for *.jar files (GATK and VarScan2)
# sometime you may need to specify full path to the binary (for freebayes, samtools and java)
GATK:      
freebayes: 
samtools:  /usr/bin/samtools
varscan:   /home/smurray/Code/VarScan.v2.4.0.jar
java:      java

[ScriptOptions]
DP_threshold:   15
number_of_SNPs:
  
  # fast_freebayes enables --targets option for Freebayes, faster but more prone to Freebayes errors
  # set to False will use --region, each variant is called separately
  fast_freebayes: True

# This is the file containing variant positions to use for genotype calling
# The format of the variant genomic positions must match the default reference (REFERENCE)
# not the alternate reference (REF_ALTERNATE)
# note: This default can be changed when calling bam-matcher with --vcf option
VCF_file: /home/smurray/Code/FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf


[VariantCallerParameters]
# GATK memory usage in GB
GATK_MEM: 4

# GATK threads (-nt)
GATK_nt:  1

# VarScan memory usage in GB
VARSCAN_MEM: 4

[GenomeReference]
# default reference fasta file
REFERENCE: /home/smurray/Code/Genomes/MacaM_Rhesus_Genome_v7.fasta
REF_ALTERNATE:
  # CHROM_MAP is required if using two different genome references that have different (but compatible) chromosome names
  # this is mainly to deal with the hg19 "chr" issue
  CHROM_MAP:
  
  [BatchOperations]
# you MUST specify a cache directory, this directory should be read/write-able by all BAM-matcher users
CACHE_DIR:/home/smurray/Code/BamMatcherCache

[Miscellaneous]

### end config file ###
#######################################################













#####  Test: Run Bam Matcher using Varscan option  #####

# Can do test run with test data provided by bam matcher.  Tiles and instructions come with bam matcher download.

# Try on two FGF21 samples with a trimmed down vcf file - SNPs from Chr01 only (should be positive result)
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.Chr01.vcf -B1 FGF21_BAM/FGF21_treated/US-1580283/US-1580283.bam -B2 FGF21_BAM/FGF21_treated/US-1580277/US-1580277.bam;
# this second test run should be a negative result (or at least lower fraction in common!)
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.Chr01.vcf -B1 FGF21_BAM/FGF21_treated/US-1580283/US-1580283.bam -B2 FGF21_BAM/FGF21_treated/US-1577042/US-1577042.bam;



# Try on those same FGF21 samples with full vcf file
# positive control:
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580283/US-1580283.bam -B2 FGF21_BAM/FGF21_treated/US-1580277/US-1580277.bam;
# negative control:
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580283/US-1580283.bam -B2 FGF21_BAM/FGF21_treated/US-1577042/US-1577042.bam;








###### make a list of samples to run from the fgf21 data #####
# can paste this into cmd line to run (as below), or figure out a way to source this list to run bamMatcher.

toRun<- read.csv("/home/smurray/Code/bam-matcher/bamMatcherSamplesToRun.csv")
toRun$codeB1 <- paste0(paste(paste0("-B1 FGF21_BAM/FGF21_treated/US-", toRun$unknown), toRun$unknown, sep="/US-"), ".bam")
toRun$codeB1[2]
toRun$codeB2 <- paste0(paste(paste0("-B2 FGF21_BAM/FGF21_treated/US-", toRun$known), toRun$known, sep="/US-"), ".bam --verbose")
toRun$codeB2[2]
toRun$codeToRun <- paste(paste(toRun$codeBeginning, toRun$codeB1, sep=" "), toRun$codeB2, sep=" ")
toRun$codeToRun[2]

write.csv(toRun, "/home/smurray/Code/bam-matcher/bamMatcherSamplesToRunFinal.csv")







##### For Reals: run bam matcher with varscan using full vcf file  #####

# chunk below this one contains all samples vs all samples.
# once each sample is run, it's variants will be cached and that cache will be used the next time the sample is run.


# ran with chr01 vcf file earlier, now running with full, depth filtered vcf file.
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580283/US-1580283.bam -B2 FGF21_BAM/FGF21_treated/US-1580277/US-1580277.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580266/US-1580266.bam -B2 FGF21_BAM/FGF21_treated/US-1580277/US-1580277.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1577093/US-1577093.bam -B2 FGF21_BAM/FGF21_treated/US-1580277/US-1580277.bam

# 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1577042/US-1577042.bam -B2 FGF21_BAM/FGF21_treated/US-1577093/US-1577093.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580279/US-1580279.bam -B2 FGF21_BAM/FGF21_treated/US-1580277/US-1580277.bam;

bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580283/US-1580283.bam -B2 FGF21_BAM/FGF21_treated/US-1580279/US-1580279.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1577083/US-1577083.bam -B2 FGF21_BAM/FGF21_treated/US-1580279/US-1580279.bam
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1577087/US-1577087.bam -B2 FGF21_BAM/FGF21_treated/US-1576044/US-1576044.bam

bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1576044/US-1576044.bam -B2 FGF21_BAM/FGF21_treated/US-1580277/US-1580277.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1577042/US-1577042.bam -B2 FGF21_BAM/FGF21_treated/US-1580277/US-1580277.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1577056/US-1577056.bam -B2 FGF21_BAM/FGF21_treated/US-1580277/US-1580277.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1577066/US-1577066.bam -B2 FGF21_BAM/FGF21_treated/US-1580277/US-1580277.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1577083/US-1577083.bam -B2 FGF21_BAM/FGF21_treated/US-1580277/US-1580277.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1577086/US-1577086.bam -B2 FGF21_BAM/FGF21_treated/US-1580277/US-1580277.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1577087/US-1577087.bam -B2 FGF21_BAM/FGF21_treated/US-1580277/US-1580277.bam;


bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580304/US-1580304.bam -B2 FGF21_BAM/FGF21_treated/US-1581502/US-1581502.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580308/US-1580308.bam -B2 FGF21_BAM/FGF21_treated/US-1581502/US-1581502.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580311/US-1580311.bam -B2 FGF21_BAM/FGF21_treated/US-1581502/US-1581502.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580315/US-1580315.bam -B2 FGF21_BAM/FGF21_treated/US-1581502/US-1581502.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581535/US-1581535.bam -B2 FGF21_BAM/FGF21_treated/US-1581502/US-1581502.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581538/US-1581538.bam -B2 FGF21_BAM/FGF21_treated/US-1581502/US-1581502.bam;

bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580304/US-1580304.bam -B2 FGF21_BAM/FGF21_treated/US-1577066/US-1577066.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580308/US-1580308.bam -B2 FGF21_BAM/FGF21_treated/US-1577066/US-1577066.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580311/US-1580311.bam -B2 FGF21_BAM/FGF21_treated/US-1577066/US-1577066.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580315/US-1580315.bam -B2 FGF21_BAM/FGF21_treated/US-1577066/US-1577066.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581535/US-1581535.bam -B2 FGF21_BAM/FGF21_treated/US-1577066/US-1577066.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581538/US-1581538.bam -B2 FGF21_BAM/FGF21_treated/US-1577066/US-1577066.bam;

bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580304/US-1580304.bam -B2 FGF21_BAM/FGF21_treated/US-1580266/US-1580266.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580308/US-1580308.bam -B2 FGF21_BAM/FGF21_treated/US-1580266/US-1580266.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580311/US-1580311.bam -B2 FGF21_BAM/FGF21_treated/US-1580266/US-1580266.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580315/US-1580315.bam -B2 FGF21_BAM/FGF21_treated/US-1580266/US-1580266.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581535/US-1581535.bam -B2 FGF21_BAM/FGF21_treated/US-1580266/US-1580266.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581538/US-1581538.bam -B2 FGF21_BAM/FGF21_treated/US-1580266/US-1580266.bam;

bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580304/US-1580304.bam -B2 FGF21_BAM/FGF21_treated/US-1577093/US-1577093.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580308/US-1580308.bam -B2 FGF21_BAM/FGF21_treated/US-1577093/US-1577093.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580311/US-1580311.bam -B2 FGF21_BAM/FGF21_treated/US-1577093/US-1577093.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580315/US-1580315.bam -B2 FGF21_BAM/FGF21_treated/US-1577093/US-1577093.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581535/US-1581535.bam -B2 FGF21_BAM/FGF21_treated/US-1577093/US-1577093.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581538/US-1581538.bam -B2 FGF21_BAM/FGF21_treated/US-1577093/US-1577093.bam;

bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580304/US-1580304.bam -B2 FGF21_BAM/FGF21_treated/US-1577087/US-1577087.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580308/US-1580308.bam -B2 FGF21_BAM/FGF21_treated/US-1577087/US-1577087.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580311/US-1580311.bam -B2 FGF21_BAM/FGF21_treated/US-1577087/US-1577087.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580315/US-1580315.bam -B2 FGF21_BAM/FGF21_treated/US-1577087/US-1577087.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581535/US-1581535.bam -B2 FGF21_BAM/FGF21_treated/US-1577087/US-1577087.bam;
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581538/US-1581538.bam -B2 FGF21_BAM/FGF21_treated/US-1577087/US-1577087.bam;

bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580304/US-1580304.bam -B2 FGF21_BAM/FGF21_treated/US-1580311/US-1580311.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580304/US-1580304.bam -B2 FGF21_BAM/FGF21_treated/US-1580308/US-1580308.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580304/US-1580304.bam -B2 FGF21_BAM/FGF21_treated/US-1580315/US-1580315.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580308/US-1580308.bam -B2 FGF21_BAM/FGF21_treated/US-1580311/US-1580311.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580308/US-1580308.bam -B2 FGF21_BAM/FGF21_treated/US-1580315/US-1580315.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580304/US-1580304.bam -B2 FGF21_BAM/FGF21_treated/US-1581538/US-1581538.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580304/US-1580304.bam -B2 FGF21_BAM/FGF21_treated/US-1581535/US-1581535.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580308/US-1580308.bam -B2 FGF21_BAM/FGF21_treated/US-1581535/US-1581535.bam 


# last two pesky samples - run them each against one sample from each animal (different reference sample from above)
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581535/US-1581535.bam -B2 FGF21_BAM/FGF21_treated/US-1577056/US-1577056.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581535/US-1581535.bam -B2 FGF21_BAM/FGF21_treated/US-1581492/US-1581492.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581535/US-1581535.bam -B2 FGF21_BAM/FGF21_treated/US-1580283/US-1580283.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581535/US-1581535.bam -B2 FGF21_BAM/FGF21_treated/US-1577042/US-1577042.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581535/US-1581535.bam -B2 FGF21_BAM/FGF21_treated/US-1581509/US-1581509.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581535/US-1581535.bam -B2 FGF21_BAM/FGF21_treated/US-1581507/US-1581507.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581535/US-1581535.bam -B2 FGF21_BAM/FGF21_treated/US-1581497/US-1581497.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581535/US-1581535.bam -B2 FGF21_BAM/FGF21_treated/US-1581514/US-1581514.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581535/US-1581535.bam -B2 FGF21_BAM/FGF21_treated/US-1581493/US-1581493.bam 

bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581535/US-1581535.bam -B2 FGF21_BAM/FGF21_treated/US-1580319/US-1580319.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581535/US-1581535.bam -B2 FGF21_BAM/FGF21_treated/US-1581492/US-1581492.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581535/US-1581535.bam -B2 FGF21_BAM/FGF21_treated/US-1581517/US-1581517.bam 

bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581538/US-1581538.bam -B2 FGF21_BAM/FGF21_treated/US-1577056/US-1577056.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581538/US-1581538.bam -B2 FGF21_BAM/FGF21_treated/US-1581492/US-1581492.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581538/US-1581538.bam -B2 FGF21_BAM/FGF21_treated/US-1580283/US-1580283.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581538/US-1581538.bam -B2 FGF21_BAM/FGF21_treated/US-1577042/US-1577042.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581538/US-1581538.bam -B2 FGF21_BAM/FGF21_treated/US-1581509/US-1581509.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581538/US-1581538.bam -B2 FGF21_BAM/FGF21_treated/US-1581507/US-1581507.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581538/US-1581538.bam -B2 FGF21_BAM/FGF21_treated/US-1581497/US-1581497.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581538/US-1581538.bam -B2 FGF21_BAM/FGF21_treated/US-1581514/US-1581514.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581538/US-1581538.bam -B2 FGF21_BAM/FGF21_treated/US-1581493/US-1581493.bam 

bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581476/US-1581476.bam -B2 FGF21_BAM/FGF21_treated/US-1577056/US-1577056.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581476/US-1581476.bam -B2 FGF21_BAM/FGF21_treated/US-1581492/US-1581492.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581476/US-1581476.bam -B2 FGF21_BAM/FGF21_treated/US-1580283/US-1580283.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581476/US-1581476.bam -B2 FGF21_BAM/FGF21_treated/US-1577042/US-1577042.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581476/US-1581476.bam -B2 FGF21_BAM/FGF21_treated/US-1581509/US-1581509.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581476/US-1581476.bam -B2 FGF21_BAM/FGF21_treated/US-1581507/US-1581507.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581476/US-1581476.bam -B2 FGF21_BAM/FGF21_treated/US-1581497/US-1581497.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581476/US-1581476.bam -B2 FGF21_BAM/FGF21_treated/US-1581514/US-1581514.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581476/US-1581476.bam -B2 FGF21_BAM/FGF21_treated/US-1581493/US-1581493.bam 

bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580266/US-1580266.bam -B2 FGF21_BAM/FGF21_treated/US-1580283/US-1580283.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1577086/US-1577086.bam -B2 FGF21_BAM/FGF21_treated/US-1580315/US-1580315.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580283/US-1580283.bam -B2 FGF21_BAM/FGF21_treated/US-1580315/US-1580315.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580277/US-1580277.bam -B2 FGF21_BAM/FGF21_treated/US-1580315/US-1580315.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580315/US-1580315.bam -B2 FGF21_BAM/FGF21_treated/US-1580283/US-1580283.bam 

bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1577093/US-1577093.bam -B2 FGF21_BAM/FGF21_treated/US-1581524/US-1581524.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1577042/US-1577042.bam -B2 FGF21_BAM/FGF21_treated/US-1580325/US-1580325.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580311/US-1580311.bam -B2 FGF21_BAM/FGF21_treated/US-1577093/US-1577093.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580325/US-1580325.bam -B2 FGF21_BAM/FGF21_treated/US-1577093/US-1577093.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1581524/US-1581524.bam -B2 FGF21_BAM/FGF21_treated/US-1577042/US-1577042.bam 

bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1577083/US-1577083.bam -B2 FGF21_BAM/FGF21_treated/US-1581487/US-1581487.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580279/US-1580279.bam -B2 FGF21_BAM/FGF21_treated/US-1581521/US-1581521.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1577083/US-1577083.bam -B2 FGF21_BAM/FGF21_treated/US-1581496/US-1581496.bam 
bam-matcher.py --config bam-matcher/Macaque.conf --vcf FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf -B1 FGF21_BAM/FGF21_treated/US-1580279/US-1580279.bam -B2 FGF21_BAM/FGF21_treated/US-1581497/US-1581497.bam 




##### Prep steps done to create necessary files #####

##### Prepare reference genome, etc #####

## Prep fasta file to use as reference:  need .dict and .fai files.
# unclear if this is necessary for varscan or just GATK (which I'm not using anymore)

# change directory 
cd Genomes

# make .dict file for reference genome
java -jar picard.jar CreateSequenceDictionary R= MacaM_Rhesus_Genome_v7.fasta O= MacaM_Rhesus_Genome_v7.dict

# make a fasta index file (.fai) for reference genome
samtools faidx MacaM_Rhesus_Genome_v7.fasta




###### Index Bam files and make a vcf file with our reference genome #####


# (these are the full-length and shortened vcf files used above in bam matcher calls.)

# doing this because we use an assembly of the macaque genome that doesn't come with any associated vcf files, 
# unlike the ensembl releases. (ensemble vcfs don't match chromosome:position notation so we can't use them.)
# there's probably away to cache the results of mpileup for each bam file while making vcf file to use later
  # in bam-matcher varscan (to reduce computational time), but i didn't go that route.

# # first use samtools to make a vcf file from each bam file
# samtools mpileup -vo FGF21_vcf/sam_US-1576044.vcf -f Genomes/MacaM_Rhesus_Genome_v7.fasta FGF21_BAM/FGF21_treated/US-1576044.bam 
# bcftools call -vmO z -o <study.vcf.gz> <study.bcf>

# try doing both of the above steps at once:   
samtools mpileup -ugf Genomes/MacaM_Rhesus_Genome_v7.fasta FGF21_BAM/FGF21_treated/US-1576044.bam | bcftools call -vmO v -o FGF21_vcf/sam_US-1576044.vcf



#Now run the same thing, using many bam files.  

# list the files in the directory
bamFiles <- list.files(path = "/home/smurray/Code/FGF21_BAM/FGF21_treated", full.names = T, recursive = T) # list all files in that directory, with full paths

#index all the bam files  (there's a better way to code this. see https://www.biostars.org/p/44511/ )
samtools index FGF21_BAM/FGF21_treated/US-1577056/US-1577056.bam    
samtools index FGF21_BAM/FGF21_treated/US-1577066/US-1577066.bam   
samtools index FGF21_BAM/FGF21_treated/US-1577073/US-1577073.bam   
samtools index FGF21_BAM/FGF21_treated/US-1577086/US-1577086.bam   
samtools index FGF21_BAM/FGF21_treated/US-1577087/US-1577087.bam    
samtools index FGF21_BAM/FGF21_treated/US-1577093/US-1577093.bam   
samtools index FGF21_BAM/FGF21_treated/US-1580266/US-1580266.bam   
samtools index FGF21_BAM/FGF21_treated/US-1580277/US-1580277.bam   
samtools index FGF21_BAM/FGF21_treated/US-1580279/US-1580279.bam   
samtools index FGF21_BAM/FGF21_treated/US-1580283/US-1580283.bam   
samtools index FGF21_BAM/FGF21_treated/US-1580286/US-1580286.bam   
samtools index FGF21_BAM/FGF21_treated/US-1580304/US-1580304.bam   
samtools index FGF21_BAM/FGF21_treated/US-1580308/US-1580308.bam   
samtools index FGF21_BAM/FGF21_treated/US-1580311/US-1580311.bam   
samtools index FGF21_BAM/FGF21_treated/US-1580315/US-1580315.bam   
samtools index FGF21_BAM/FGF21_treated/US-1580319/US-1580319.bam  
samtools index FGF21_BAM/FGF21_treated/US-1580325/US-1580325.bam   
samtools index FGF21_BAM/FGF21_treated/US-1581466/US-1581466.bam  
samtools index FGF21_BAM/FGF21_treated/US-1581471/US-1581471.bam   
samtools index FGF21_BAM/FGF21_treated/US-1581472/US-1581472.bam 
samtools index FGF21_BAM/FGF21_treated/US-1581473/US-1581473.bam  
samtools index FGF21_BAM/FGF21_treated/US-1581475/US-1581475.bam 
samtools index FGF21_BAM/FGF21_treated/US-1581476/US-1581476.bam 
samtools index FGF21_BAM/FGF21_treated/US-1581477/US-1581477.bam  
samtools index FGF21_BAM/FGF21_treated/US-1581478/US-1581478.bam  
samtools index FGF21_BAM/FGF21_treated/US-1581487/US-1581487.bam  
samtools index FGF21_BAM/FGF21_treated/US-1581489/US-1581489.bam   
samtools index FGF21_BAM/FGF21_treated/US-1581492/US-1581492.bam   
samtools index FGF21_BAM/FGF21_treated/US-1581493/US-1581493.bam   
samtools index FGF21_BAM/FGF21_treated/US-1581494/US-1581494.bam   
samtools index FGF21_BAM/FGF21_treated/US-1581496/US-1581496.bam   
samtools index FGF21_BAM/FGF21_treated/US-1581497/US-1581497.bam   
samtools index FGF21_BAM/FGF21_treated/US-1581499/US-1581499.bam   
samtools index FGF21_BAM/FGF21_treated/US-1581502/US-1581502.bam   
samtools index FGF21_BAM/FGF21_treated/US-1581504/US-1581504.bam   
samtools index FGF21_BAM/FGF21_treated/US-1581505/US-1581505.bam   
samtools index FGF21_BAM/FGF21_treated/US-1581506/US-1581506.bam   
samtools index FGF21_BAM/FGF21_treated/US-1581507/US-1581507.bam   
samtools index FGF21_BAM/FGF21_treated/US-1581508/US-1581508.bam   
samtools index FGF21_BAM/FGF21_treated/US-1581509/US-1581509.bam   
samtools index FGF21_BAM/FGF21_treated/US-1581513/US-1581513.bam   
samtools index FGF21_BAM/FGF21_treated/US-1581514/US-1581514.bam   
samtools index FGF21_BAM/FGF21_treated/US-1581516/US-1581516.bam   
samtools index FGF21_BAM/FGF21_treated/US-1581517/US-1581517.bam   
samtools index FGF21_BAM/FGF21_treated/US-1581521/US-1581521.bam   
samtools index FGF21_BAM/FGF21_treated/US-1581524/US-1581524.bam   
samtools index FGF21_BAM/FGF21_treated/US-1581526/US-1581526.bam   
samtools index FGF21_BAM/FGF21_treated/US-1581535/US-1581535.bam   
samtools index FGF21_BAM/FGF21_treated/US-1581538/US-1581538.bam   

# some of these didn't transfer well- will have to re-download from server
#missed the few below- run again
samtools index FGF21_BAM/FGF21_treated/US-1581499/US-1581499.bam 
samtools index FGF21_BAM/FGF21_treated/US-1577083/US-1577083.bam 
samtools index FGF21_BAM/FGF21_treated/US-1581538/US-1581538.bam
samtools index FGF21_BAM/FGF21_treated/US-1580308/US-1580308.bam
samtools index FGF21_BAM/FGF21_treated/US-1577083/US-1577083.bam


# missed some more
samtools index FGF21_BAM/FGF21_treated/US-1581507/US-1581507.bam
samtools index FGF21_BAM/FGF21_treated/US-1581509/US-1581509.bam
samtools index FGF21_BAM/FGF21_treated/US-1581514/US-1581514.bam
samtools index FGF21_BAM/FGF21_treated/US-1581517/US-1581517.bam
samtools index FGF21_BAM/FGF21_treated/US-1581521/US-1581521.bam



# run pileup with all those files    
# **!!! ADD A FILTER FOR READ DEPTH ON PER-SAMPLE BASIS!!!
# -D outputs per-sample read depth
samtools mpileup -ugf Genomes/MacaM_Rhesus_Genome_v7.fasta \
FGF21_BAM/FGF21_treated/US-1576044.bam \
FGF21_BAM/FGF21_treated/US-1577042/US-1577042.bam \
FGF21_BAM/FGF21_treated/US-1577056/US-1577056.bam \
FGF21_BAM/FGF21_treated/US-1577066/US-1577066.bam \
FGF21_BAM/FGF21_treated/US-1577073/US-1577073.bam \
FGF21_BAM/FGF21_treated/US-1577083/US-1577083.bam \
FGF21_BAM/FGF21_treated/US-1577086/US-1577086.bam \
FGF21_BAM/FGF21_treated/US-1577087/US-1577087.bam \
FGF21_BAM/FGF21_treated/US-1577093/US-1577093.bam \
FGF21_BAM/FGF21_treated/US-1580266/US-1580266.bam \
FGF21_BAM/FGF21_treated/US-1580277/US-1580277.bam \
FGF21_BAM/FGF21_treated/US-1580279/US-1580279.bam \
FGF21_BAM/FGF21_treated/US-1580283/US-1580283.bam \
FGF21_BAM/FGF21_treated/US-1580286/US-1580286.bam \
FGF21_BAM/FGF21_treated/US-1580304/US-1580304.bam \
FGF21_BAM/FGF21_treated/US-1580308/US-1580308.bam \
FGF21_BAM/FGF21_treated/US-1580311/US-1580311.bam \
FGF21_BAM/FGF21_treated/US-1580315/US-1580315.bam \
FGF21_BAM/FGF21_treated/US-1580319/US-1580319.bam \
FGF21_BAM/FGF21_treated/US-1580325/US-1580325.bam \
FGF21_BAM/FGF21_treated/US-1581466/US-1581466.bam \
FGF21_BAM/FGF21_treated/US-1581471/US-1581471.bam \
FGF21_BAM/FGF21_treated/US-1581472/US-1581472.bam \
FGF21_BAM/FGF21_treated/US-1581473/US-1581473.bam \
FGF21_BAM/FGF21_treated/US-1581475/US-1581475.bam \
FGF21_BAM/FGF21_treated/US-1581476/US-1581476.bam \
FGF21_BAM/FGF21_treated/US-1581477/US-1581477.bam \
FGF21_BAM/FGF21_treated/US-1581478/US-1581478.bam \
FGF21_BAM/FGF21_treated/US-1581487/US-1581487.bam \
FGF21_BAM/FGF21_treated/US-1581489/US-1581489.bam \
FGF21_BAM/FGF21_treated/US-1581492/US-1581492.bam \
FGF21_BAM/FGF21_treated/US-1581493/US-1581493.bam \
FGF21_BAM/FGF21_treated/US-1581494/US-1581494.bam \
FGF21_BAM/FGF21_treated/US-1581496/US-1581496.bam \
FGF21_BAM/FGF21_treated/US-1581497/US-1581497.bam \
FGF21_BAM/FGF21_treated/US-1581499/US-1581499.bam \
FGF21_BAM/FGF21_treated/US-1581502/US-1581502.bam \
FGF21_BAM/FGF21_treated/US-1581504/US-1581504.bam \
FGF21_BAM/FGF21_treated/US-1581505/US-1581505.bam \
FGF21_BAM/FGF21_treated/US-1581506/US-1581506.bam \
FGF21_BAM/FGF21_treated/US-1581508/US-1581508.bam \
FGF21_BAM/FGF21_treated/US-1581524/US-1581524.bam \
FGF21_BAM/FGF21_treated/US-1581526/US-1581526.bam \
FGF21_BAM/FGF21_treated/US-1581535/US-1581535.bam \
FGF21_BAM/FGF21_treated/US-1581538/US-1581538.bam | bcftools call -vmO v -o FGF21_vcf/sam_mostFGF21_treated.vcf

## The above worked!!  Took 3.5 days on my 1st linux desktop without any multi-threading, etc.






##### trim vcf file using bcftools ####

# can't get VariantAnnotation or VcfR packages to work to look at VCF file in R...
## try vcf tools (perl in cmd line) or bcftools (cmd line)  **bcftools is updated version of vcftools

## common settings:  (pulled from a bioconductor discussion)
# MinDP (Minimum read depth):   5 (Indels) and 3 (SNPs)
# MaxDP (Maximum read depth):  "You have a low coverage data, so I would set it to 100. Normally it is 3 times the average coverage."
# BaseQualBias (Minimum p-value for baseQ bias):  0
# MinMQ (Minimum RMS mapping quality for SNPs):  20 or 30 (to be more stringent)
# Qual (Minimum value of QUAL field):  15 or 20
# 
# StrandBias (Minimum p-value for strand bias):  0.0001
# EndDistBias (Minimum p-value for end distance bias):  0.0001
# MapQualBias (Minimum p-value for mapQ bias):  0
# VBD (Minimum Variant Distance Bias):  0 (More relevant to RNA-seq reads)
# 
# GapWin (Window size for filtering adjacent gaps):  30 bp
# SnpGap (SNP within INT bp around a gap to be filtered):   20 bp
# 
# SNPcluster (number of snps within a region): I usually drop all the snps if there are more than 3 snps within 10 bp. 
  # I don't see how to actually implement SNPcluster in the filter...


# redo path to get to correct version of bcftools:
export PATH=/home/smurray/Code/bcftools-1.4.1/bin:$PATH 
export PATH=/home/smurray/Code/samtools-1.4.1/bin:$PATH 
cd Code

# attempt filter with bcftools
bcftools filter -e "QUAL<15" FGF21_vcf/sam_mostFGF21_copy.vcf > FGF21_vcf/sam_mostFGF21_filt1.vcf 
bcftools filter -e "DP<3" FGF21_vcf/sam_mostFGF21_filt1.vcf > FGF21_vcf/sam_mostFGF21_filt2.vcf 
# average DP in first chromosome is 153, 3x153= 450.  Let's use 500.
bcftools filter -e "DP>500" FGF21_vcf/sam_mostFGF21_filt2.vcf > FGF21_vcf/sam_mostFGF21_filt3.vcf 
bcftools filter -e "MQ<20" FGF21_vcf/sam_mostFGF21_filt3.vcf > FGF21_vcf/sam_mostFGF21_filt4.vcf 
bcftools filter -e "BQB<0" FGF21_vcf/sam_mostFGF21_filt4.vcf > FGF21_vcf/sam_mostFGF21_filt5.vcf 
bcftools filter -e "MQB<0" FGF21_vcf/sam_mostFGF21_filt5.vcf > FGF21_vcf/sam_mostFGF21_filt6.vcf 
bcftools filter -e "VDB<0" FGF21_vcf/sam_mostFGF21_filt6.vcf > FGF21_vcf/sam_mostFGF21_filt7.vcf 
bcftools filter --SnpGap 10 FGF21_vcf/sam_mostFGF21_filt7.vcf > FGF21_vcf/sam_mostFGF21_filt8.vcf   #998MB
bcftools filter -e 'TYPE="INDEL"' FGF21_vcf/sam_mostFGF21_filt8.vcf > FGF21_vcf/sam_mostFGF21_filt9.vcf   #936MB
bcftools filter -e "QUAL<50" FGF21_vcf/sam_mostFGF21_filt10.vcf > FGF21_vcf/sam_mostFGF21_filt11.vcf   #448MB
bcftools filter -e "DP<50" FGF21_vcf/sam_mostFGF21_filt11.vcf > FGF21_vcf/sam_mostFGF21_filt.DP50.vcf  #133MB


# make a second vcf, with baseline quality filters and high depth coverage filter
bcftools filter --SnpGap 10 FGF21_vcf/sam_mostFGF21_copy.vcf > FGF21_vcf/sam_mostFGF21_filt2.1.vcf   #2.2GB
bcftools filter -e 'TYPE="INDEL"' FGF21_vcf/sam_mostFGF21_filt2.1.vcf > FGF21_vcf/sam_mostFGF21_filt2.2.vcf  #2.1GB
bcftools filter -e "QUAL<20" FGF21_vcf/sam_mostFGF21_filt2.2.vcf > FGF21_vcf/sam_mostFGF21_filt2.3.vcf   # 1.1GB
# average DP in first chromosome is 153, 3x153= 450.  Let's use 500.
bcftools filter -e "DP>500" FGF21_vcf/sam_mostFGF21_filt2.3.vcf > FGF21_vcf/sam_mostFGF21_filt2.4.vcf   # 986MB
bcftools filter -e "MQ<20" FGF21_vcf/sam_mostFGF21_filt2.4.vcf > FGF21_vcf/sam_mostFGF21_filt2.5.vcf   # 8.4MB
bcftools filter -e "BQB<0" FGF21_vcf/sam_mostFGF21_filt2.5.vcf > FGF21_vcf/sam_mostFGF21_filt2.6.vcf    # 8.4MB
bcftools filter -e "MQB<0" FGF21_vcf/sam_mostFGF21_filt2.6.vcf > FGF21_vcf/sam_mostFGF21_filt2.7.vcf    # 8.4MB
bcftools filter -e "VDB<0" FGF21_vcf/sam_mostFGF21_filt2.7.vcf > FGF21_vcf/sam_mostFGF21_filt2.8.vcf    # 8.4MB
# high depth coverage filtering
bcftools filter -e "DP<50" FGF21_vcf/sam_mostFGF21_filt2.8.vcf > FGF21_vcf/sam_mostFGF21_filt2.9.vcf   #155 MB, 
bcftools filter -e "DP<100" FGF21_vcf/sam_mostFGF21_filt2.9.vcf > FGF21_vcf/sam_mostFGF21_filt.DP100.vcf   #90 MB, 
bcftools filter -e "DP<150" FGF21_vcf/sam_mostFGF21_filt.DP100.vcf > FGF21_vcf/sam_mostFGF21_filt.DP150.vcf   #62 MB, 
bcftools filter -e "DP<200" FGF21_vcf/sam_mostFGF21_filt.DP150.vcf > FGF21_vcf/sam_mostFGF21_filt.DP200.vcf   #45 MB, 
bcftools filter -e "DP>400" FGF21_vcf/sam_mostFGF21_filt.DP200.vcf > FGF21_vcf/sam_mostFGF21_filt.200DP400.vcf   #35 MB, 50k SNPs
bcftools filter -e "DP>400" FGF21_vcf/sam_mostFGF21_filt.DP100.vcf > FGF21_vcf/sam_mostFGF21_filt.100DP400.vcf   #80 MB, 117k SNPs 
bcftools filter -e "QUAL<50" FGF21_vcf/sam_mostFGF21_filt.100DP400.vcf > FGF21_vcf/sam_mostFGF21_filt.DP100.400.Q50.vcf  #70 MB, 102k SNPs
bcftools filter -e "QUAL<50" FGF21_vcf/sam_mostFGF21_filt.200DP400.vcf > FGF21_vcf/sam_mostFGF21_filt.DP200.400.Q50.vcf  #32 MB


# made a shortened VCF file of only variants on chromosome 1 (for test runs) just by trimming vcf text file in text editor.
  # file = sam_mostFGF21_filt.DP100.400.Q50.Chr01.vcf, has ~10k SNPs

# now use vcf file in bam-matcher tool.






