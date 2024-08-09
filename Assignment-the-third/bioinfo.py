#!/usr/bin/env python

# Author: Evan Tizzard tizzard@uoregon.edu

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "0.4"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = "ATCGNatcgn"
RNA_bases = "AUCGNaucgn"

def convert_phred(letter: str) -> int:
    '''Converts a single ASCII character into its corresponding Phred+33 score.'''
    return ord(letter)-33

def qual_score(phred_score: str) -> float:
    '''Calculates the average quality score of an entire Phred string.'''
    qual_sum = 0
    for character in (phred_score):
        qual_sum += convert_phred(character)
    return qual_sum/len(phred_score)

def validate_base_seq(seq:str, RNAflag:bool=False) -> bool:
    '''This function takes a sequence string and flag. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    DNA = set(DNA_bases)
    RNA = set(RNA_bases)
    return set(seq)<=(RNA if RNAflag else DNA)

def gc_content(seq: str) -> float:
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    gc_count = 0
    g_or_c = "GgCc"
    for character in seq:
        if character in g_or_c:
            gc_count += 1
    GC_percent = gc_count/len(seq)
    assert validate_base_seq(seq), "not a DNA string"
    return GC_percent

def calc_median(lst: list) -> float:
    '''Sorts then returns the median of a one-dimensional list.'''
    lst.sort()
    list_length = int(len(lst))
    
    if list_length % 2 == 0:
        median_position1 = int(list_length/2)
        median_position2 = int(median_position1 - 1)
        median = (lst[median_position1] + lst[median_position2]) / 2

    if list_length % 2 != 0:
        median_position = list_length//2
        median = lst[median_position]
    return median

def oneline_fasta(input_fasta: str) -> bool:
    '''Function takes an input FASTA file with multiple sequence lines per record and generates a FASTA file with one sequence
    line per one header line with the naming convention "oneseqlineperrecord_<inputfilename>; returns True if half the number of 
    output file lines is equal to the number of records in the input file."'''

    #initializing first line flag and line count variables
    first_line = True
    record_count = 0
    output_lines = 0

    #opening output and input files for writing/reading
    with open("oneseqlineperrecord_"+input_fasta, "w") as output_fasta, open(input_fasta, "r") as fasta_file:

        # for every line in input file, strip newlines
        for line in fasta_file:
            line = line.strip()
            
            #when first line is encountered, write the first line to output
            #file and set first line flag to false, adding newline to end
            if first_line:
                output_fasta.write(line+"\n")
                first_line = False
                record_count += 1
            
            #else if any line in file begins with ">" (is a header line),
            #print line to file and add leading and tailing newlines
            elif line.startswith(">"):
                record_count += 1
                output_fasta.write("\n"+line+"\n")
                
            #finally when sequence lines are encountered, print with no newlines
            #to force concatenation until another header line is encountered
            else:
                output_fasta.write(line)
    
    with open("oneseqlineperrecord_"+input_fasta, "r") as output_fasta:
        for line in output_fasta:
            output_lines += 1
        output_records = output_lines / 2
    return output_records == record_count

def reverse_compliment(seq: str) -> str:
    '''Takes a nucleotide sequence and returns its reverse compliment. If sequence contains N, reverse complimentary base will also be N.'''
    seq_rc = ""
    for character in seq[::-1]:
        if character == "A":
            seq_rc += "T"
        elif character == "T":
            seq_rc += "A"
        elif character == "G":
            seq_rc += "C"
        elif character == "C":
            seq_rc += "G"
        elif character == "N":
            seq_rc += "N"
    return seq_rc



if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    # These tests are run when you execute this file directly (instead of importing it)
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("convert_phred function working properly")

    assert qual_score("A") == 32.0, "wrong average phred score for 'A'"
    assert qual_score("AC") == 33.0, "wrong average phred score for 'AC'"
    assert qual_score("@@##") == 16.5, "wrong average phred score for '@@##'"
    assert qual_score("EEEEAAA!") == 30.0, "wrong average phred score for 'EEEEAAA!'"
    assert qual_score("$") == 3.0, "wrong average phred score for '$'"
    print("qual_score function working properly")
    
    assert validate_base_seq("AATAGAT"), "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True), "Validate base seq does not work on RNA"
    assert validate_base_seq("R is the best!")==False, "Not a DNA string"
    assert validate_base_seq("aatagat"), "Validate base seq does not work on lowercase DNA"
    assert validate_base_seq("aauagau", True), "Validate base seq does not work on lowercase RNA"
    assert validate_base_seq("TTTTtttttTTT")
    print("validate_base_seq function working properly")

    assert gc_content("GCGCGC") == 1, "wrong GC content for string entirely composed of G and C"
    assert gc_content("AATTATA") == 0, "wrong GC content for string composed of no Gs or Cs"
    assert gc_content("GCATCGAT") == 0.5, "wrong GC content for DNA string composed of half Gs or Cs"
    print("gc_content function working properly")
    
    assert calc_median([1,2,3]) == 2, "wrong median for odd-numbered list"
    assert calc_median([5,6,7,8]) == 6.5, "wrong median for even-numbered list"
    assert calc_median([1,1,1,1,1,1,1,1,100]) == 1, "wrong median for heavily duplicated value list"
    assert calc_median([7]) == 7, "wrong median for single value"
    assert calc_median([50,100]) == 75, "wrong median for list of 2 values"
    print("calc_median function working properly")

    assert oneline_fasta("Unit_test_01.fa") == True, "output file does not contain one-line sequences"
    print("oneline_fasta function working properly")