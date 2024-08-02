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
    1.  ![Read1 Plot](https://github.com/ettizzard/Demultiplex/blob/9dffbc06dbe2ab774616887783d08fb19cc1e156/Assignment-the-first/1294_S1_L008_R1_001_distributionplot.png)  
    ![Read2 Plot](https://github.com/ettizzard/Demultiplex/blob/868d00b9ff43eeaec423377a50fecd46163aeab0/Assignment-the-first/1294_S1_L008_R2_001_distributionplot.png)  
    ![Read3 Plot](https://github.com/ettizzard/Demultiplex/blob/868d00b9ff43eeaec423377a50fecd46163aeab0/Assignment-the-first/1294_S1_L008_R3_001_distributionplot.png)  
    ![Read4 Plot](https://github.com/ettizzard/Demultiplex/blob/868d00b9ff43eeaec423377a50fecd46163aeab0/Assignment-the-first/1294_S1_L008_R4_001_distributionplot.png)  
      
    2. For index reads, we need to be more strict on quality. If even a single base in an index sequence is incorrect, the entire read will be improperly labelled. Improper read labeling is mucn more detrimental to overall downstream analysis than a single or few bases in a biological read being incorrect. So, I would want an index quality score cutoff of at least 30, as this is 99.9% accurate and is generally sufficient for homology. However, these quality scores seem to be binned (i.e., most assigned quality scores are given A, F, or J, meaning any scores in between these are encompassed within them by binning as a means of data compression). My cutoff of 30 is encompassed by "A", so a score of 32. This is more strict than my original cutoff, but I think 32 is perfectly appropriate as well. So, any index sequence with a single base with a quality score below 32 will be sorted to the low quality/unknown index fastq file pair. For our purposes, we actually don't really need a strict cutoff for biological reads, because our downstream analysis invloves genome assembly. If a biological read is of too low quality or contains too many "N"s, it simply will not align. If I had to provide a cutoff for biological reads, however, I would just want the average quality score of the entire read to be at least 30.  
      
    3. `$ zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | awk 'NR % 4 == 2' | grep "N" | wc -l`  
        7304664
    
## Part 2
1. Define the problem  
    Our original biological read fastq files contain reads undistinguished by index. We need to organize all of these reads into pairs of new fastq files based on their read number and index matching status. One pair of new fastq files needs to include all reads with matching index pairs. A second pair of fastq files needs to include index-hopped reads. The final pair of fastq files needs to include reads with unknown/low quality indices.
2. Describe output  
    Assuming we use every index provided, we will have 24 pairs of fastq files containing reads with kniwn dual matched indices. File1 in the pair will be biological read1, and file 2 will be biological read2. Another pair of fastq files will contain all of the reads that exhibited index hopping of known indices, again file1 being biological read1 and file2 being biological read2. The final pair of fastq files will contain all reads with unknown or low quality corresponding indices, file1 for biological read1 and file2 for biological read2. All said and done, we should have 52 separate fastq files. Additionally, we will output the number of read pairs with properly match indices per index pair, the number of index-hopped read pairs, and the number of read pairs with low quality/unknown indices.  
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode
```
#Initialization
    #Shebang

    #Import necessary modules

    #Create a variable to hold count of number of pairs of properly matched indices, another to hold count of index-hopped pairs, and a third variable to hold count of pairs with low quality/unknown indices. All 3 of these variables will be initially set to 0.

    #Create a set to contain our 24 index sequences.
        #To accomplish this, open the index file from Talapas and iterate over to add each index to set.
    
    #Create a dictionary to contain tuples of index pairs as keys and tuples of file names as values.
        #Format for dual-matched pairs: keys = ("index1", "index2"); values =  ("R1_dualmatched_{index1}-{index2}.fq", "R2_dualmatched_{index1}-{index2}.fq")
        #Format for index-hopped pairs: keys = ("index1", "index2"); values =  ("R1_indexhopped.fq", "R2_indexhopped.fq")
        #Format for pairs with unknown/low-quality indices: keys = ("index1", "index2"); values =  ("R1_lowqual_unknown.fq", "R2_lowqual_unknown.fq")

    #Create another dictionary to keep count of occurrences of index pairs. Dictionary keys will be a tuple of the index sequences, and values will be the number of occurrences of the index pair.
        #Before our main body of code, we can initialize the 24 dual-matched index pairs as keys in this dictionary with values set to 0. As we continue in the algorithm, we will add unmatched high-quality index pairs to this dictionary as they are encountered.

#Main Body
    #Start by opening our 4 original fastq files at once.
    
    #Iterate over all 4 fastq files 1 record at a time (4 lines).
    
    #Per record, take each record line per file and save as corresponding variable.
    
    #Take the read3 index sequence, and create its reverse compliment, updating the variable with this string.
    
    #Append the read1 and read4 header lines with the index1 and index2 (read 2 and read3 reverse compliment).

    #Now that we have our record in the proper format, we can sort the reads.
        #First check to see if indices are within known list. If not, sort reads to unknown/low quality fastq files and increment unknown/low quality counter
        #Convert quality score lines of read1 and read4. Calculate average qs of entire read for both. If Average is below cutoff, also sort to unknown/low quality fastq files and increment unknown/low quality counter.
        
        #Elif:
            #Check that indices match. If they don't, add read1 and read4 records with appended headers to index-hopped fastq file pair and incrememnt index-hopped counter.
                #Check if this index pair already exists within the index pair freq dictionary. If it does, increment the corresponding value. Otherwise, add the pair and set its value to 1.

        #Elif:
            #This statement should encompass our high quality dual-matched index reads.
            #Increment matched index counter and the corresponding value in the index pair freq dictionary.
            #Add record to corresponding index in the first index pair-filename dictionary.

#Outputs
    #Print index pair-filename dictionary entries into corresponding files.
    #Print counts.
    #Close all files left open.
```
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
```
def reverse_compliment(seq: str) -> str:
    '''Takes a nucleotide sequence and returns its reverse compliment. If sequence contains N, reverse comlimentary base will also be N.'''
    return seq_rc
Input: ATCGTAT
Expected output: ATACGAT

def print_outputs(dual_matched_count: int, index_hopped_count: int, low_qual_unknown_count: int, paired_index_dict: dict, outputfilename: str) -> str:
    '''Takes counts of index pair types as well as the dictionary containing index pairs and filename conventions, then returns corresponding str to be written to file.'''
    return fileoutput
Input: above variables and dictionary
Expected output: properly formatted tab-separated summary of results to be written to file

def print_sorted_pairs(read1_header: str, read1_seq: str, read1_qscoreline: str, read1_filename: str, read4_header: str, read4_seq: str, read4_qscoreline: str, read4_filename: str) -> str:
    '''Rewrites each line of read1 and read4 records to sorted fastq file pairs based on index.
Input: above variables
Expected output: string of lines formatted in the fastq file convention
```