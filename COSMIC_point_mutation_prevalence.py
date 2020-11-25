# Author: Jacob Haninmann
# Date: 20.11.20
# Python 2.7

Gene_Set = raw_input('Do you want to analyse a specific gene-set? (True/False):')


import csv
import time
from collections import OrderedDict
import copy
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from sklearn.ensemble import IsolationForest
import re
import pickle
import pandas as pd
from csv import writer
from csv import reader
import random
import statistics


#Read Cosmic Data and Drop Duplicates
def cosmic_data_freq(filename):
    print('Files bigger than 1 GB could take more than 30 minutes, please wait...')
    cosmic_data = pd.read_csv(filename, sep='\t')
    cosmic_mut_file_wo_dups = cosmic_data.drop_duplicates(
        subset=["Gene name", "Accession Number", "Gene CDS length", "Primary site", "Mutation CDS", "Mutation AA",
                "Mutation genome position", "Age", 'Sample name'])
    return cosmic_mut_file_wo_dups

# timer
start_time = time.time()

if Gene_Set == "True":
# create nested dictionary for gene-set panel with value= None:
    gene_set_dict = {}
    with open(raw_input("Please enter path of gene-set file: ")) as fh:
        for line in fh:
            gene_set_dict = dict.fromkeys(line.split(), {'count': 0})


print("Removing duplicates from the dataset...") #takes approximately 30 min for 1 GB file
cosmic_without_duplicates = cosmic_data_freq(raw_input('Please enter path of COSMIC Mutant file: '))
cosmic_without_duplicates.to_csv(r'Cosmic_wo_duplicated_samples', index=False, sep='\t')  #converting pandas file to csv


print('Extracting, selecting and counting COSMIC data...')
gene_mut_count = {}
i = 0  #visualize progress


#transcript counter
transcript_counter= 0
multiple_transcripts=[]
with open("Cosmic_wo_duplicated_samples") as cosmic_db:
    for line in cosmic_db:
        data = str(line).split('\t')
        gene = data[0].replace('b\'', '')
        if "_" in gene:
            continue
        if Gene_Set=="True": #selecting for gene-set
            if gene not in gene_set_dict:
                continue
        transcript= data[1]
        snp = data[27]  # SNP info
        if snp == 'y':  # skip SNPs
            continue
        # Counting mutations in nested gene_set dictionary
        mutation = data[20]  # AA mutation
        if mutation == '':  # filter unknown mutations
            continue
        if mutation == 'p.?': #filter unknown mutations
            continue
        # deleting the "p." string before mutation name
        pattern = '[^p\.]+'
        mutation = re.findall(pattern, mutation)[0]
        #saving just the position of the mutation, deleting specific change in the amino acid sequence
        pattern = '[A-Z]\d+'
        match = re.findall(pattern, mutation)
        if match == []:
            continue
        else:
            mutation= match[0]


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
                print('new isoform of: '+gene)
                print(transcript)
                multiple_transcripts.append(gene)
        else: #gene never seen before
            gene_mut_count[gene]={'count':1}
            gene_mut_count[gene][transcript]={mutation:1}

    # visualize progress
        i += 1
        if (i % 100000 == 0):
            print('processed', round(i / 100000, 0), '100k lines')

print('It occured '+str(transcript_counter)+' times where another transcript was identified')


print("Calculating frequency of each mutation type...(per gene)")

gene_mut_freq = copy.deepcopy(gene_mut_count) #copy dictionary as independent dict
#frequency is calculated with total gene count, not per transcript
#iterating through 3x nested dictionary to calculate the relative frequency of each mutation and saving results in new dict.
#Syntax is: key1= gene, key2= transcript,count, key3=mutation, value3=relative frequency of each mutation
for gene, transcript_list in gene_mut_count.items():
    for transcript,mutation_list in transcript_list.items()[1:]: #starting at second because value2 of count is an integer
        for mutation,count in mutation_list.items():
            gene_mut_freq[gene][transcript][mutation]= round(float(count)/ (gene_mut_count[gene]['count']),3)
    gene_mut_freq[gene]['count'] = gene_mut_count[gene]['count']


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


print('Assign a rank to each mutation...(rank is per transcript, not per gene)')
#Rank is per transcript and not per gene
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


# write results in csv file
i = 0
name= raw_input('Please name output file (format is "name"_mutation_frequency.csv): ')
print('Writing results to .tsv file...')
with open(name+'_mutation_frequencies.tsv', 'w') as output_file:
    tsv_writer = csv.writer(output_file, delimiter='\t')
    tsv_writer.writerow(['mutation', 'gene', 'transcript', 'mutation_count', 'Gene_count', 'relative_frequency', 'rank', 'rank_score', 'total_ranks', 'relative_position'])
#iterate through the gene_mut_count dictionary while also extracting from the gene_count_freq dictionary
    for gene, transcript_dict in ordered_gene_mut_count.items():
        gene_count = gene_mut_count[gene]['count'] #Gene count
        for transcript, mutation_dict in transcript_dict.items():
            for mutation,count in mutation_dict.items():
                relative_frequency = ordered_gene_mut_freq[gene][transcript][mutation][0]
                rank = ordered_gene_mut_freq[gene][transcript][mutation][1]
                total_ranks = ordered_gene_mut_freq[gene][transcript]['rank count']
                relative_position = round(1-ordered_gene_mut_freq[gene][transcript][mutation][2],3)
                rank_score = round((relative_frequency+rank),3)
                tsv_writer.writerow([mutation, gene, transcript, count, gene_count,relative_frequency,rank,rank_score, total_ranks,relative_position])
        i += 1
        if (i % 10000 == 0):
            print('processed', round(i / 10000, 0), '10k lines')

#count how many different mutations were identified
file = open(name+"_mutation_frequencies.tsv", "r")
number_of_lines = 0
for line in file:
  line = line.strip("\n")

  number_of_lines += 1
file.close()

print(str(number_of_lines) + " mutations were identified.")


print('Cashing data in textfiles...')

with open("position_gene_mut_count.txt", "wb") as fp:   #Pickling
    pickle.dump(gene_mut_count, fp)

with open("position_gene_mut_freq.txt", "wb") as fp:   #Pickling
    pickle.dump(gene_mut_freq, fp)

with open("position_ordered_gene_mut_freq.txt", "wb") as fp:   #Pickling
    pickle.dump(OrderedDict(ordered_gene_mut_freq), fp)

with open("position_ordered_gene_mut_count.txt", "wb") as fp:   #Pickling
    pickle.dump(OrderedDict(ordered_gene_mut_count), fp)

#--------------------------------------------------------------------------------------------------

print("Load data from textfiles...")

with open("position_gene_mut_count.txt", "rb") as fp:   #Pickling
    gene_mut_count = pickle.load(fp)

with open("position_gene_mut_freq.txt", "rb") as fp:   #Pickling
    gene_mut_freq = pickle.load(fp)

with open("position_ordered_gene_mut_freq.txt", "rb") as fp:   #Pickling
    ordered_gene_mut_freq = pickle.load(fp)

with open("position_ordered_gene_mut_count.txt", "rb") as fp:   #Pickling
    ordered_gene_mut_count = pickle.load(fp)


#Classification and visualisation of the data

#Selecting mode: 'Plot' or 'Save'
mode = raw_input("Do you want to plot the data or save the classification as a csv file? (Plot/Save) ")

if mode == 'Plot':
    red = 'red'
    green = 'green'
    blue = 'blue'
else: #to save Classification in the file
    red = 1
    green = 2
    blue = 3


#ignore warnings of scitkit isolation forest function
import warnings
warnings.filterwarnings("ignore")


# timer
start_time = time.time()


print('Classifying, Plotting Data and Saving...')

hotspot_classification = []

plt.style.use('ggplot')
random = (random.randint(1,20)) #to randomly plot genes to avoid biased overfitting
iteration = 1
for gene, transcript_dict in ordered_gene_mut_count.items(): #iterating through the list to extract values
    print('Mutation entries for '+gene+': '+str(gene_mut_count[gene]['count']))
    for transcript in transcript_dict.keys():
        mutation = ordered_gene_mut_freq[gene][transcript].keys()[:-1]
        frequency = []
        col = []  # list to store classification for plot
        frequency_tree =[] #for trainingdata set, has to be in the format [[1],[2],[2]] for example
        count = [] #extracting the discrete count of each mutation, will be used for labeling the plot
        wholecount = [] #saving the discrete count of each mutation, will be used for classification
        for i in range(0, (len(ordered_gene_mut_freq[gene][transcript].values()) - 1)):
            frequency.append(ordered_gene_mut_freq[gene][transcript].values()[i][0]) #frequency extraction
            frequency_tree.append([ordered_gene_mut_freq[gene][transcript].values()[i][0]]) #data for isolation forest
            total_count = len(ordered_gene_mut_count[gene][transcript].keys()) #for count label per mutation
            wholecount.append(ordered_gene_mut_count[gene][transcript].values()[i]) #save count of each mutation
            steps = int(0.05 *total_count+1) #steps which should be displayed as labels
            if (i % steps == 0):
                count.append(ordered_gene_mut_count[gene][transcript].values()[i]) #get count of mutation

        remove_value= min(frequency_tree) #removing unwanted data for trainingdata set
        number = 0
        for i in frequency_tree:
            if i == remove_value:
                number +=1
        if (len(frequency_tree) -number) > 3: #make sure input for isolation forest is bigger than three
            frequency_tree = frequency_tree[:-number] #is short for [0:len(l)-1]
        strong_outlier = False
        try:
            if frequency_tree[0][0] > (4*float(frequency_tree[1][0])): #if highest frequency is 4 times higher as the second one, remove from trainingdata
                frequency_tree= frequency_tree[1:]
                col= [red] #classify highest value as Hot
                strong_outlier = True
        except:
            pass

    #Create the isolation forest model and train and test.
        clf = IsolationForest(random_state=0,bootstrap=False).fit(frequency_tree)
        outlier_score= clf.decision_function(frequency_tree) #there is also predict as a method (value 1 or -1) output
        outlier_classification = clf.predict(frequency_tree) # -1 is outlier and 1 is inliner


        #label the results of the classification
        x_pos = [i for i, _ in enumerate(mutation)]
        #if there is no negative value, blau und green
        outcome_score = sum(1 for number in outlier_score if number < 0) #checking if there is a negative value
        outcome_classification = sum(1 for number in outlier_classification if number < 0) #checking if classification thinks the same because sometimes a outliner-score accidentally swaps to a negative number


        #coloring the data according to the outliner score generated with following criteria
        if outcome_score > 0 and outcome_classification > 0: #make sure that both predictions have the same amount of outliners
            b = -1  # keeping track of indices
            switch_green = 0  #isolation forest sometimes classifies data clearly wrong (issue of the method apparently): switch will be activated so that after the first classification of green, only green and blue can follow to the end
            switch_blue = 0  # isolation forest sometimes classifies data clearly wrong: switch will be activated so that after the first classification of blue, only blue can follow to the end
            for val in outlier_score:
               b += 1 #val score is corresponding for wholecount[b]
               if outlier_score[b] <0 and outlier_classification[b] <0 and frequency[b] >= (10*statistics.median(frequency)): #double check that both prediction calculated the same score for a mutation
                    if val <= 0.85 *min(outlier_score) and switch_green==0 and switch_blue==0: #to be a hotspot, the discrete count has to be 20 or more and be in the top 33% of outliners
                         col.append(red)
                    elif val <= 0.3 *min(outlier_score) and frequency[b] >= (5*statistics.median(frequency)) and switch_blue==0: #threshold for warm mutations is set to at least 10 counts
                         col.append(green)
                         switch_green=1
                    else: #if score is not in the top 66% outliners, it is a cold mutation
                         col.append(blue)
                         switch_blue =1  #after the first classifcation as blue, only blue can follow because ranking is descendend
               else:
                   col.append(blue)  #if results are not consistent (predict and decision_function), all mutations are unsignificanct, thus cold mutations
                   switch_blue=1


        else: # if there are negative values, respectively when there are no outliners
            for val in outlier_score:
                    col.append(blue)

        if strong_outlier == False:
            appendix = len(frequency) - len(frequency_tree)  #counting how many datapoints were not used for training
        else: appendix = len(frequency) - (len(frequency_tree)+1)

        for i in range(0, appendix): #coloring data that was not used for training the model
            col.append(blue)


        #make list for classifications to later add as a column to the frequency file

        hotspot_classification.extend(col)
        if mode== "Plot":
           if iteration % random==0:  # to avoid bias of looking at the same plots over and over again
               plt.bar(x_pos, frequency, color=col)
               plt.xlabel("Mutations total: " + str(len(x_pos)))
               plt.ylabel("Frequency")
               plt.title("Classification for " + str(gene) + ' mutations')

               #label with count every 10th mutation
               x_pos_anot = []
               for i in range(0, steps * len(count), steps):
                   x_pos_anot.append(i)
               plt.xticks(x_pos_anot, count) #labeling
               countred = sum(map(lambda x: x == 'red', col)) #print the mutation names, get index
               countgreen = sum(map(lambda x: x == 'green', col)) #get index
               print('red'+str(mutation[0:countred]))
               print('green' + str(mutation[countred:(countred+countgreen)]))
               plt.show()
               iteration +=1
           else:
               iteration +=1


if mode =='Save':

    print('Add new column <Classification> to existing frequency file...')

    def add_column_in_csv(input_file, output_file, transform_row):
        """ Append a column in existing csv using csv.reader / csv.writer classes"""
        # Open the input_file in read mode and output_file in write mode
        with open(input_file, 'r') as read_obj, \
                open(output_file, 'w') as write_obj:
            # Create a csv.reader object from the input file object
            csv_reader = reader(read_obj, delimiter='\t')
            # Create a csv.writer object from the output file object
            csv_writer = writer(write_obj, delimiter='\t')
            # Read each row of the input csv file as list
            for row in csv_reader:
                # Pass the list / row in the transform function to add column text for this row
                transform_row(row, csv_reader.line_num)
                # Write the updated row / list to the output file
                csv_writer.writerow(row)


    #add column to frequency file
    header_of_new_col = 'Classification'
    add_column_in_csv(name+'_mutation_frequencies.tsv', name+'_frequencies_namitag_Cdx_w_classification.csv', lambda row, line_num: row.append(header_of_new_col)
                                                                if line_num == 1 else row.append(hotspot_classification[line_num - 2]))
    print('Finished')

print('It is recommended to only keep one transcript per gene for consistency...')
print('Genes with multiple isoforms are the following:')
print(multiple_transcripts)