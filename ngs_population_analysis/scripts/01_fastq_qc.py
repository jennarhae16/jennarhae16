#------ Assignment 01 FASTQ ------ 
#libraries needed
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO


"Functions for code"
def ImportFastq(file_in):
    with open (file_in) as file:        #import file
        quality_scores=[]               #empty list for quality scores
        for record in SeqIO.parse(file, "fastq"):
            quality = record.letter_annotations['phred_quality'] #allows ASCII code to be converted to Q-Score list
            quality_scores.append(quality)

    return quality_scores

def mean_of_Q_scores(list):
    arr = np.asarray(list)      #convert nested list to multidimentional array
    mean_array = np.mean(arr, axis = 0) #take means of each column
    
    return mean_array



fastq_list = ImportFastq("assignment1.fastq")
mean_array = mean_of_Q_scores(fastq_list)

x = range(len(mean_array))
y = mean_array
plt.title("Avg Q Score by position")
plt.xlabel("Position")
plt.ylabel("Average Q score")
plt.plot(x,y)
plt.show()



