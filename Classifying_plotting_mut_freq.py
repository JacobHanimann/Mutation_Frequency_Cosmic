#Selecting mode: 'Plot' or 'Save'
mode = 'Plot'

if mode == 'Plot':
    red = 'red'
    green = 'green'
    blue = 'blue'
else: #to save Classification in the file
    red = 1
    green = 2
    blue = 3


#Calculation of frequencies of selected mutations from cosmic for the tumour board online tool

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
from csv import writer
from csv import reader
import random
import statistics



#ignore warnings because scitkit isolation forest function
import warnings
warnings.filterwarnings("ignore")


# start time count
start_time = time.time()


print('Reading dictionaries from first script into this one...')
#Saving lists to text file to be used in other script

with open("position_gene_mut_count.txt", "rb") as fp:   #Pickling
    gene_mut_count = pickle.load(fp)

with open("position_gene_mut_freq.txt", "rb") as fp:   #Pickling
    gene_mut_freq = pickle.load(fp)

with open("position_ordered_gene_mut_freq.txt", "rb") as fp:   #Pickling
    ordered_gene_mut_freq = pickle.load(fp)

with open("position_ordered_gene_mut_count.txt", "rb") as fp:   #Pickling
    ordered_gene_mut_count = pickle.load(fp)


print('Classifying, Plotting Data and Saving...')

hotspot_classification = []

plt.style.use('ggplot')
z = 0 # to visualize how many plots were already generated
random = (random.randint(1,20))
print(random)
for gene, transcript_dict in ordered_gene_mut_count.items(): #iterating through the list to extract values
    for transcript in transcript_dict.keys():
        mutation = ordered_gene_mut_freq[gene][transcript].keys()[:-1]
        frequency = []
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
            col = [] #list to store classification for plot
        strong_outlier = False
        if frequency_tree[0][0] > (4*float(frequency_tree[1][0])): #if highest frequency is 4 times higher as the second one, remove from trainingdata
            frequency_tree= frequency_tree[1:]
            col= [red] #classify highest value as Hot
            strong_outlier = True

    #Create the isolation forest model and train and test.
    clf = IsolationForest(random_state=0,bootstrap=False).fit(frequency_tree) #have to be moved to the left when there is more than one transcript
    outlier_score= clf.decision_function(frequency_tree) #there is also predict as a method (value 1 or -1) output
    outlier_classification = clf.predict(frequency_tree) # -1 is outlier and 1 is inliner


    #label the results of the classification
    x_pos = [i for i, _ in enumerate(mutation)]
    #if there is no negative value, blau und green
    outcome_score = sum(1 for number in outlier_score if number < 0) #checking if there is a negative value
    outcome_classification = sum(1 for number in outlier_classification if number < 0) #checking if classification thinks the same because sometimes outliner-score accidentally swaps to negative numbers


    #coloring the data according to the outliner score generated with following criteria
    if outcome_score > 0 and outcome_classification > 0: #make sure that both predictions have the same amount of outliners
        b = -1  # keeping track of indices
        switch_green = 0  #isolation forest sometimes classifies data clearly wrong (problem of the method I guess): switch will be activated so that after the first classification of green, only green and blue can follow to the end
        switch_blue = 0  # isolation forest sometimes classifies data clearly wrong: switch will be activated so that after the first classification of blue, only blue can follow to the end
        for val in outlier_score:
           b += 1 #val score is corresponding for wholecount[b]
           if outlier_score[b] <0 and outlier_classification[b] <0 and frequency[b] >= (10*statistics.median(frequency)): #double check that both prediction calculated the same score for a mutation
                if val <= 0.85 *min(outlier_score) and switch_green==0 and switch_blue==0: #to be a hotspot, discrete count has  be 20 or more and be top 33% of outliner
                     col.append(red)
                elif val <= 0.3 *min(outlier_score) and frequency[b] >= (5*statistics.median(frequency)) and switch_blue==0: #threshold for warm mutations is set to at least 10 counts, 0.6 isches gsi
                     col.append(green)
                     switch_green=1
                else: #if score is not in the top 66% outliners, it is a cold mutation
                     col.append(blue)
                     switch_blue =1  #after the first classifcation as blue, only blue can follow because ranking is descendend
           else:
               col.append(blue)  #if results are not consistent (predict and decision_function), all mutations are unsignificanct thus cold mutations
               switch_blue=1

    # maybe no es else da
    else: # if there are negative values, respectively when there are no outliners
        for val in outlier_score:
                col.append(blue)

    if strong_outlier == False:
        appendix = len(frequency) - len(frequency_tree)  #counting how many datapoints were not used for training
    else: appendix = len(frequency) - (len(frequency_tree)+1)

    for i in range(0, appendix): #coloring data that wasnt used for training model
        col.append(blue)

#make list for classifications to later add as a column to the frequency file

    hotspot_classification.extend(col)
    z += 1 #counting how many plots were already made
    print(z)

    if mode== "Plot":
       if gene =='BRAF' or gene == 'KRAS' or gene == 'EGFR' or gene == 'CTNNB1' or gene == 'ERBB2' or gene == 'FGFR2' or gene == 'IDH1' or gene == 'MET' or gene == 'MTOR' or gene == 'ROS1':  # to avoid bias of looking at the same plots over and over again
           plt.bar(x_pos, frequency, color=col)
           plt.xlabel("Mutations total: " + str(len(x_pos)))
           plt.ylabel("Frequency")
           plt.title("Classification for " + str(gene) + ' mutations')

           # label with count every 10th mutation
           x_pos_anot = []
           for i in range(0, steps * len(count), steps):
               x_pos_anot.append(i)
           plt.xticks(x_pos_anot, count) #labeling
           countred = sum(map(lambda x: x == 'red', col)) #print the mutation names, get index
           countgreen = sum(map(lambda x: x == 'green', col)) #get index
           print('red'+str(mutation[0:countred]))
           print(countgreen)
           print('green' + str(mutation[countred:(countred+countgreen)]))
           plt.show()
           #plt.savefig("Classification_" + str(gene) + "_mutations.png")


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
    add_column_in_csv('position_mutation_frequencies.tsv', 'position_frequencies_Cdx_final.csv', lambda row, line_num: row.append(header_of_new_col)
                                                                if line_num == 1 else row.append(hotspot_classification[line_num - 2]))
    print('Finished')