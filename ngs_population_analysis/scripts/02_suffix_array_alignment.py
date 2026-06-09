#------ Assignment 02 Suffix Array ------

#libraries needed
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser

'''Functions for code'''

def ImportFastaSeq(file_in):
    with open (file_in) as file:                        #import file
        for record in SimpleFastaParser(file):
            seq = record[1]                             
##            print(seq[0:5])
    return seq

def suffix_array(input_str):                            #sorts genome into suffix array of indices
    SA = sorted(range(len(input_str)), key=lambda i: input_str[i:])
##    print(SA[0:5])
    return SA

def findbisection(pattern, text, suf_arr):
    low,lo,high,hi = 0,0,len(text),len(text)
    bounds_pair = []
    
    while (low < high):                                 #finds lower bound where first read occurs
        midpoint = int((low+high)/2)
        if text[suf_arr[midpoint]:] < pattern: 
            low = midpoint + 1
        else:
            high = midpoint
##    print(low)        
    bounds_pair.append(low)
              
    while (lo < hi):                                    #finds upperbound where read occurence ends
        middle = int((lo+hi)/2)
        if text[suf_arr[middle]:suf_arr[middle]+len(pattern)] <= pattern:
            lo = middle + 1
        else:
            hi = middle

##    print(lo)       
    bounds_pair.append(lo)
    
    return bounds_pair

def output_info(bounds_list, Suf_ar, Genome, Read):     #print location of reads and sequences
    i = bounds_list[0]
    with open("output.txt", "w") as file:
        text = f"Read:{Read}\nSuffix Array Index:{bounds_list[0]} to {bounds_list[1]}\nRead Locations: {SA[read_bounds[0]:read_bounds[1]]}\n"
        file.write(text)
        while i < bounds_list[1]:
            text = f"Location:{Suf_ar[i]}\nSequence:{Genome[Suf_ar[i]:]}\n"
            file.write(text)
            i +=1
    
'''Function Calls'''


##genome = "ATCGTCAGTACGATGCTGGGATACTTAGATAAGCAATTGCTCGAT"
##read = "AG"


genome = ImportFastaSeq("Assignment2_refgenome.fasta")
read = ImportFastaSeq("Assignment2_read.fasta")
SA = suffix_array(genome)
read_bounds = findbisection(read,genome,SA)
output_info(read_bounds, SA, genome, read)



##print(read_bounds)
##
##i = read_bounds[0]
##print(SA[read_bounds[0]:read_bounds[1]])
##while i < read_bounds[1]:
##    print(SA[i], genome[SA[i]:])
##    i +=1 



