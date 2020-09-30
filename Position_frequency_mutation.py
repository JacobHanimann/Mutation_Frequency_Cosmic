FoundationOne_Gene_Set = True

#Calculation of frequencies of selected mutations from cosmic for the tumour board online tool

# improvements:
#downloading data vith API

import csv
import gzip
import xlrd
import time
import pickle
import pprint
import math
from collections import OrderedDict
import copy
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from sklearn.ensemble import IsolationForest
import re
import pickle
import pandas as pd
import argparse



# Read Comic Data and Drop Duplicates
def cosmic_data_freq(filename):
    cosmic_data = pd.read_csv(filename, sep='\t')
    cosmic_mut_file_wo_dups = cosmic_data.drop_duplicates(
        subset=["Gene name", "Accession Number", "Gene CDS length", "Primary site", "Mutation CDS", "Mutation AA",
                "Mutation genome position", "Age", 'Sample name'])
    return cosmic_mut_file_wo_dups



# start time count
start_time = time.time()

if FoundationOne_Gene_Set == True:
    print('Selecting genes from the gene-panel of FoundationOnex')
# creating a list of the 324 genes from the gene panel

# create nested dictionary for gene-set panel (too many apparently) with value= None:
    gene_set_foundation = {}
    with open("genepanel.txt") as fh: #current version: 336 units, probably because of alternative gene names
        for line in fh:
            gene_set_foundation = dict.fromkeys(line.split(), {'count': 0})

print("Removing duplicates from the dataset...") #takes approximately 30 min
cosmic_without_duplicates = cosmic_data_freq('CosmicMutantExport.tsv.gz')  #remove duplicates
cosmic_without_duplicates.to_csv(r'Cosmic_wo_duplicated_samples', index=False, sep='\t')  #converting panda file to csv, maybe better zip it?
# extract and select mutations from the cosmic dataset

print('Extracting, selecting and counting COSMIC data...')
# match each COSMIC mutation to NCI numbers by matching the tissues using nci_to_tissues
# write matches in dictionary nci_to_genes
gene_mut_count = {}
i = 0  # to visualize progress while computing

#transcript counter
transcript_counter= 0
with open("Cosmic_wo_duplicated_samples") as cosmic_db:
    for line in cosmic_db:  # each line goes through the loop
        data = str(line).split('\t')
        gene = data[0].replace('b\'', '')  # gene name, remove unicode, how does it work?
        if "_" in gene:
            continue
        if FoundationOne_Gene_Set==True: #selecting for selected gene-set from FoundationOne
            if gene not in gene_set_foundation:
                continue
        transcript= data[1]
        snp = data[27]  # SNP info
        if snp == 'y':  # skip SNPs
            continue
        # Counting mutations in nested gene_set dictionary
        mutation = data[20]  # AA mutation
        if mutation == '':  # filter unknown mutations which are not the same (position not even known)
            continue
        if mutation == 'p.?':  # filter unknown mutation which are not the same (position mostly known), could be deleted probably because of regular expression selection below
            continue
        # deleting the "p." string before mutation name with regular expression
        pattern = '[^p\.]+'
        mutation = re.findall(pattern, mutation)[0]
        #saving just the position of the mutation, deleting specific change in the amino acid sequence
        pattern = '[A-Z]\d+'
        match = re.findall(pattern, mutation)
        if match == []:
            continue
        else:
            mutation= match[0]
        #maybe count how many unknown mutations there are...

    #counting number of mutation for each gene and each transcript resulting in a 3x nested dictionary: gene_mut_count
    #Syntax is: key1= gene, key2= transcript,count, key3=mutation, value3=count of specific mutation
        if gene in gene_mut_count.keys(): #gene seen before
            gene_mut_count[gene]['count'] += 1
            if transcript in gene_mut_count[gene]: #transcript seen before
                if mutation in gene_mut_count[gene][transcript]: #mutation seen before
                    gene_mut_count[gene][transcript][mutation] +=1
                else: #mutation never seen before
                    gene_mut_count[gene][transcript][mutation]=1
            else: #transcript never seen before
                gene_mut_count[gene][transcript]={mutation:1}
                transcript_counter +=1
        else: #gene never seen before
            gene_mut_count[gene]={'count':1}
            gene_mut_count[gene][transcript]={mutation:1}

    # visualize progress
        i += 1
        if (i % 100000 == 0):
            print('processed', round(i / 100000, 0), '100k lines')

print('It occured '+str(transcript_counter)+' times where another transcript was identified')


print("Calculating frequency of each mutation type...")

gene_mut_freq = copy.deepcopy(gene_mut_count) #copy dictionary as independent
#frequency is calculated with total gene count
#iterating through 3x nested dictionary to calculate the relative frequency of each mutation and saving results in new dict.
#Syntax is: key1= gene, key2= transcript,count, key3=mutation, value3=relative frequency of each mutation
for gene, transcript_list in gene_mut_count.items():
    for transcript,mutation_list in transcript_list.items()[1:]: #starting at second because value2 of count is an integer
        for mutation,count in mutation_list.items():
            gene_mut_freq[gene][transcript][mutation]= round(float(count)/ (gene_mut_count[gene]['count']),3)
    gene_mut_freq[gene]['count'] = gene_mut_count[gene]['count']  # changing the count back to true value, delete?


print('Reordering the mutations by descending frequency...') #creating new sorted dictionary
ordered_gene_mut_freq = {}
#Syntax is: key1= gene, key2= transcript dict, key3= mutation_dict, value3= relativ mutation frequencies
for gene in gene_mut_freq:
    ordered_gene_mut_freq[gene] = {} #dict must be created before inserting an object to it
    for transcript,mutation_dict in gene_mut_freq[gene].items()[1:]: #starting at second index because value2 of count is an integer
        ordered = OrderedDict(sorted(mutation_dict.items(), key=lambda i: i[1], reverse=True))
        ordered_gene_mut_freq[gene][transcript] = ordered

#reordering gene count dictionary
ordered_gene_mut_count = {}
ordered = 0
# Syntax is: key1= gene, key2= transcript dict, key3= mutation_dict, value3= relativ mutation frequencies
for gene in gene_mut_count:
    ordered_gene_mut_count[gene] = {}  # dict must be created before inserting an object to it
    for transcript, mutation_dict in gene_mut_count[gene].items()[1:]:  # starting at second index because value2 of count is an integer
        ordered = OrderedDict(sorted(mutation_dict.items(), key=lambda i: i[1], reverse=True))
        ordered_gene_mut_count[gene][transcript] = ordered


print('Assign a rank to each mutation...')
#iterating through list and assigning the ranking so that value2 of dict is a list of freq and rank number
#Syntax is: key1= gene, key2= transcript_dict, key3= mutation_dict, "rank number" value3 = [relative frequency, rank number,relativ rank]

for gene, transcript_dict in ordered_gene_mut_freq.items():
    s= 0 #total rank count
    for transcript, mutation_dict in transcript_dict.items():
        freq_value_before = 0 #to compare the frequency to the last frequency to decide if the rank should be the same or higher
        i=1
        for mutation in mutation_dict.keys():
            ordered_gene_mut_freq[gene][transcript][mutation]= [ordered_gene_mut_freq[gene][transcript][mutation]]
            if ordered_gene_mut_freq[gene][transcript][mutation][0] < freq_value_before: #compare to previous frequency
                i +=1
                ordered_gene_mut_freq[gene][transcript][mutation].append(i) #append rank
            else:
                ordered_gene_mut_freq[gene][transcript][mutation].append(i) #the frequency is the same as the previous one
            freq_value_before = ordered_gene_mut_freq[gene][transcript][mutation][0]
            s +=1
        rank_count= s #save total rank number for a gene


#calculate and append relative rank number
        for mutation in mutation_dict.keys():
            ordered_gene_mut_freq[gene][transcript][mutation].append(round(ordered_gene_mut_freq[gene][transcript][mutation][1] / float(rank_count),3))
            ordered_gene_mut_freq[gene][transcript]['rank count']= rank_count #append total rank number to the transcript


#Saving lists to text file to be used in other script
#Syntax is: key1= gene, key2= transcript_dict, key3= mutation_dict, "rank number" value3 = [relative frequency, rank number,relativ rank]

with open("position_gene_mut_count.txt", "wb") as fp:   #Pickling
    pickle.dump(gene_mut_count, fp)

with open("position_gene_mut_freq.txt", "wb") as fp:   #Pickling
    pickle.dump(gene_mut_freq, fp)

with open("position_ordered_gene_mut_freq.txt", "wb") as fp:   #Pickling
    pickle.dump(OrderedDict(ordered_gene_mut_freq), fp)

with open("position_ordered_gene_mut_count.txt", "wb") as fp:   #Pickling
    pickle.dump(OrderedDict(ordered_gene_mut_count), fp)


print('Writing results to .tsv file...')
# write results in csv file
i = 0
with open('position_mutation_frequencies.tsv', 'w') as output_file:
    tsv_writer = csv.writer(output_file, delimiter='\t')
    tsv_writer.writerow(['mutation', 'gene', 'transcript', 'mutation_count', 'Gene_count', 'relative_frequency', 'rank', 'rank_score', 'total_ranks', 'relative_position'])
# iterate through the gene_mut_count dictionary while also extracting from the gene_count_freq dictionary
    for gene, transcript_dict in ordered_gene_mut_count.items():
        gene_count = gene_mut_count[gene]['count'] #Gene count
        for transcript, mutation_dict in transcript_dict.items():
            for mutation,count in mutation_dict.items():
                relative_frequency = ordered_gene_mut_freq[gene][transcript][mutation][0]
                rank = ordered_gene_mut_freq[gene][transcript][mutation][1]
                total_ranks = ordered_gene_mut_freq[gene][transcript]['rank count']
                relative_position = round(1-ordered_gene_mut_freq[gene][transcript][mutation][2],3)
                rank_score = round((relative_frequency+rank),3)
                #insert variables per row
                tsv_writer.writerow([mutation, gene, transcript, count, gene_count,relative_frequency,rank,rank_score, total_ranks,relative_position])
        i += 1
        if (i % 10000 == 0):
            print('processed', round(i / 10000, 0), '10k lines')

#count how many mutations were identified
file = open("position_mutation_frequencies.tsv", "r")
number_of_lines = 0
for line in file:
  line = line.strip("\n")

  number_of_lines += 1
file.close()

print(str(number_of_lines) + " mutations were identified.")

print("--- %s seconds ---" % (time.time() - start_time))


print('Classifying, Plotting Data and Saving...')
#
#plt.style.use('ggplot')
#z = 0 # to visualize how many plots were already generated
#for gene, transcript_dict in ordered_gene_mut_freq.items(): #iterating through the list to extract values
#    for transcript in transcript_dict.keys():
#        mutation = ordered_gene_mut_freq[gene][transcript].keys()[:-1]
#        frequency = []
#        frequency_tree =[] #for trainingdata set, has to be in the format [[1],[2],[2]] for example
#        count = [] #extracting the discrete counting values i suppose
#        for i in range(0, (len(ordered_gene_mut_freq[gene][transcript].values()) - 1)):
#            frequency.append(ordered_gene_mut_freq[gene][transcript].values()[i][0]) #frequency extraction
#            frequency_tree.append([ordered_gene_mut_freq[gene][transcript].values()[i][0]]) #data for isolation forest
#            total_count = len(ordered_gene_mut_count[gene][transcript].keys()) #for count label per mutation
#            steps = int(0.05 *total_count+1) #steps which should be displayed as labels
#            if (i % steps == 0):
#                count.append(ordered_gene_mut_count[gene][transcript].values()[i]) #get count of mutation
#
#        remove_value= min(frequency_tree) #removing unwanted data for trainingdata set
#        number = 0
#        for i in frequency_tree:
#            if i == remove_value:
#                number +=1
#        if (len(frequency_tree) -number) > 3: #make sure input for isolation forest is bigger than three
#            frequency_tree = frequency_tree[:-number] #is short for [0:len(l)-1]
#
#    #Create the isolation forest model and train and test.
#    clf = IsolationForest(random_state=0,bootstrap=False).fit(frequency_tree) #have to be moved to the left when there is more than one transcript
#    outlier_score= clf.decision_function(frequency_tree) #there is also predict as a method (value 1 or -1) output
#
#
#    #label the results of the classification
#    x_pos = [i for i, _ in enumerate(mutation)]
#    col = [] #coloring data
#    #if there is no negative value, blau und green
#    outcome = sum(1 for number in outlier_score if number < 0) #checking if there is a negative value #number problem?
#    #check outcome of outcome by printing the array and the outcome lol
#
#    #coloring the data according to the outliner score generated with following criteria
#    if outcome > 0: #wenn es negativi values het
#     for val in outlier_score:
#        if val < 0.666 *min(outlier_score):
#             col.append('red')
#        elif val > 0.666 *max(outlier_score):
#            col.append('blue')
#        else:
#            col.append('green')
#    else: # if there are negative values, respectively when there are no outliners
#        for val in outlier_score:
#            if val < 0.666 * max(outlier_score):
#                col.append('blue')
#            else:
#                col.append('green')
#
#    appendix = len(frequency) - len(frequency_tree)  #counting how many datapoints were not used for training
#
#    for i in range(0, appendix): #coloring data that wasnt used for training model
#        col.append('purple')
#
##plotting data
#    plt.bar(x_pos, frequency, color=col)
#    plt.xlabel("Mutations total: " + str(len(x_pos)))
#    plt.ylabel("Frequency")
#    plt.title("Classification for " + str(gene) + ' mutations')
#
#    # label with count every 10th mutation
#    x_pos_anot = []
#    for i in range(0, steps * len(count), steps):
#        x_pos_anot.append(i)
#    plt.xticks(x_pos_anot, count) #labeling
#
#    plt.show()
#    # plt.savefig("Classification_" + str(gene) + "_mutations.png")
#    print(z)
#    z += 1 #counting how many plots were already made
#
#print("--- %s seconds ---" % (time.time() - start_time))