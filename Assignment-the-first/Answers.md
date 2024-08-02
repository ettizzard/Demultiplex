# Assignment the First

## Part 1
1. Be sure to upload your Python script. Provide a link to it here:

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read1 | 101 | +33 |
| 1294_S1_L008_R2_001.fastq.gz | index1 | 8 | +33 |
| 1294_S1_L008_R3_001.fastq.gz | index2 | 8 | +33 |
| 1294_S1_L008_R4_001.fastq.gz | read2 | 101 | +33 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.  
      
    2. For index reads, we need to be more strict on quality. If even a single base in an index sequence is incorrect, the entire read will be improperly labelled. Improper read labeling is mucn more detrimental to overall downstream analysis than a single or few bases in a biological read being incorrect. So, I would want an index quality score cutoff of at least 30, as this is 99.9% accurate and is generally sufficient for homology. However, these quality scores seem to be binned (i.e., most assigned quality scores are given A, F, or J, meaning any scores in between these are encompassed within them by binning as a means of data compression). My cutoff of 30 is encompassed by "A", so a score of 32. This is more strict than my original cutoff, but I think 32 is perfectly appropriate as well. So, any index sequence with a single base with a quality score below 32 will be sorted to the low quality/unknown index fastq file pair. For our purposes, we actually don't really need a strict cutoff for biological reads, because our downstream analysis invloves genome assembly. If a biological read is of too low quality or contains too many "N"s, it simply will not align. If I had to provide a cutoff for biological reads, however, I would just want the average quality score of the entire read to be at least 30.  
      
    3. `$ zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | awk 'NR % 4 == 2' | grep "N" | wc -l`  
        7304664
    
## Part 2
1. Define the problem
    Our original biological read fastq files contain reads undistinguished by index. We need to organize all of these reads into pairs of new fastq files based on their read number and index matching status. One pair of new fastq files needs to include all reads with matching index pairs. A second pair of fastq files needs to include index-hopped reads. The final pair of fastq files needs to include reads with unknown/low quality indices.
2. Describe output
    Assuming we use every index provided, we will have 24 pairs of fastq files containing reads with kniwn dual matched indices. File1 in the pair will be biological read1, and file 2 will be biological read2. Another pair of fastq files will contain all of the reads that exhibited index hopping of known indices, again file1 being biological read1 and file2 being biological read2. The final pair of fastq files will contain all reads with unknown or low quality corresponding indices, file1 for biological read1 and file2 for biological read2.
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
