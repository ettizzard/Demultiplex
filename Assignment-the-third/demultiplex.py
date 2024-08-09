#!/usr/bin/env python

#Importing necessary modules
import bioinfo
import matplotlib.pyplot as plt
import gzip
import argparse
import itertools

#Initializing argparse
def get_args():
    parser = argparse.ArgumentParser(description="Setting up fastq and read length arguments")
    parser.add_argument("-r1", help="input r1 fastq file", type=str)
    parser.add_argument("-r2", help="input r2 fastq file", type=str)
    parser.add_argument("-r3", help="input r3 fastq file", type=str)
    parser.add_argument("-r4", help="input r4 fastq file", type=str)
    parser.add_argument("-i", help="known index file", type=str)
    return parser.parse_args()
args = get_args()

#Starting counts for index pair types
index_hopped_pairs_count = 0
lowqual_unknown_pairs_count = 0

#Creating and populating set to contain known index sequences
indices = set()
with open(args.i, "r") as index_textfile:
    for line_num, line in enumerate(index_textfile):
        line = line.strip()
        line = line.split("\t")
        if line_num > 0:
            indices.add(line[4])


#Initializing filename and pair count dictionaries for dual matched indices.
dualmatched_filename_dict = {} # {keys: (index1, index2); values: ({index1}-{index2}_R1.fq, {index1}-{index2}_R2.fq)}
dualmatched_pair_count = {} # {keys: (index1, index2); values: count of index pair occurences}


#For loop iterating over list of indices (index1) to generate their reverse compliments, index2 in the case of dual-matching. Then populating dictionaries
#to hold indices, filenames, and counts of index pair occurrences.
for index1 in indices:
    index2 = index1
    dualmatched_filename_dict[(index1, index2)] = (open(str(index1)+"-"+str(index2)+"_R1.fq", "w"), open(str(index1)+"-"+str(index2)+"_R2.fq", "w"))
    dualmatched_pair_count[(index1, index2)] = 0

#Opening index hopped and unknown/lowquality read fq files for writing
r1_index_hopped_fq = open("index_hopped_R1.fq", "w")
r2_index_hopped_fq = open("index_hopped_R2.fq", "w")
r1_unknown_lowqual_fq = open("unknown_lowqual_R1.fq", "w")
r2_unknown_lowqual_fq = open("unknown_lowqual_R2.fq", "w")

#Setting a record count for later stats
total_record_count = 0

#Now starting to iterate record by record in all 4 fastq files.
with gzip.open(args.r1, "rt") as read1, gzip.open(args.r2, "rt") as read2, gzip.open(args.r3, "rt") as read3, gzip.open(args.r4, "rt") as read4:
    while True:
        #Iterating 4 lines at a time per iteration
        r1_header = read1.readline().strip()
        r1_sequence = read1.readline().strip()
        r1_plus = read1.readline().strip()
        r1_qscore = read1.readline().strip()

        #Breaking loop when reaching EOF
        if r1_header == "":
            break
        
        total_record_count += 1
        r2_header = read2.readline().strip()
        r2_sequence = read2.readline().strip()
        r2_plus = read2.readline().strip()
        r2_qscore = read2.readline().strip()

        r3_header = read3.readline().strip()
        r3_sequence = read3.readline().strip()
        #Creating the reverse compliment of R3 to verify it is the same as R2 for dual matching
        r3_sequence_rc = bioinfo.reverse_compliment(r3_sequence)
        r3_plus = read3.readline().strip()
        r3_qscore = read3.readline().strip()

        r4_header = read4.readline().strip()
        r4_sequence = read4.readline().strip()
        r4_plus = read4.readline().strip()
        r4_qscore = read4.readline().strip()

        #appending header lines with index sequences
        indexed_r1_header = r1_header + " " + r2_sequence + "-" + r3_sequence_rc
        indexed_r4_header = r4_header + " " + r2_sequence + "-" + r3_sequence_rc

            
        #Sorting known reads
        if r2_sequence in indices and r3_sequence_rc in indices:
            
            #Sorting low-quality reads
            both_index_qscores = r2_qscore + r3_qscore
            read_written_flag = False
            for base in both_index_qscores:
                if bioinfo.convert_phred(base) < 32:
                    r1_unknown_lowqual_fq.write(indexed_r1_header + "\n" + r1_sequence + "\n" + "+\n" + r1_qscore + "\n")
                    r2_unknown_lowqual_fq.write(indexed_r4_header + "\n" + r4_sequence + "\n" + "+\n" + r4_qscore + "\n")
                    lowqual_unknown_pairs_count += 1
                    read_written_flag = True
                    break
            
            if read_written_flag:
                continue
            
            #Sorting dual matched reads
            elif r2_sequence == r3_sequence_rc:
                
                #Setting file handles based on indices
                r1_dualmatched_fq = dualmatched_filename_dict[(r2_sequence, r2_sequence)][0]
                r2_dualmatched_fq = dualmatched_filename_dict[(r2_sequence, r2_sequence)][1]
            
                r1_dualmatched_fq.write(indexed_r1_header + "\n" + r1_sequence + "\n" + "+\n" + r1_qscore + "\n")
                r2_dualmatched_fq.write(indexed_r4_header + "\n" + r4_sequence + "\n" + "+\n" + r4_qscore + "\n")
                dualmatched_pair_count[(r2_sequence, r3_sequence_rc)] += 1
                
            
            #Sorting index-hopped reads
            elif r2_sequence != r3_sequence_rc:
                r1_index_hopped_fq.write(indexed_r1_header + "\n" + r1_sequence + "\n" + "+\n" + r1_qscore + "\n")
                r2_index_hopped_fq.write(indexed_r4_header + "\n" + r4_sequence + "\n" + "+\n" + r4_qscore + "\n")
                index_hopped_pairs_count += 1
            

        #Sorting unknown reads
        else:
            r1_unknown_lowqual_fq.write(indexed_r1_header + "\n" + r1_sequence + "\n" + "+\n" + r1_qscore + "\n")
            r2_unknown_lowqual_fq.write(indexed_r4_header + "\n" + r4_sequence + "\n" + "+\n" + r4_qscore + "\n")
            lowqual_unknown_pairs_count += 1
  

#Output print statements
for key in dualmatched_pair_count:
    print(f"{key[0]}-{key[1]}\t{dualmatched_pair_count[key]}\t{dualmatched_pair_count[key]/total_record_count*100:.2f}%")

print(f"Number of Index-Hopped Pairs: {index_hopped_pairs_count}\nPercentage of Read Pairs Exhibiting Index-Hopping: {index_hopped_pairs_count/total_record_count*100:.2f}%")
print(f"Number of Unknown or Low Quality Pairs: {lowqual_unknown_pairs_count}\nPercentage of Reads Containing Unknown or Low Quality Indices: {lowqual_unknown_pairs_count/total_record_count*100:.2f}%")


#Closing manually opened files
for key in dualmatched_filename_dict:
    dualmatched_filename_dict[key][0].close()
    dualmatched_filename_dict[key][1].close()


r1_index_hopped_fq.close()
r2_index_hopped_fq.close()
r1_unknown_lowqual_fq.close()
r2_unknown_lowqual_fq.close()