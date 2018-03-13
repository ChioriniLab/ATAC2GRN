import numpy as np
from sklearn import metrics as skm
import argparse
import random
import math
import os
import csv
import matplotlib.pyplot as plt

random.seed(10000)

def read_tfpn_bedfile(file):
    
    # This is a function that reads a .bed file
    # of the open chromatin intervals and relays
    # that to a list object.
    
    # Here we open the .bed file.
    tfpn_bedfile = open(file,'rU')
    
    # Here we initialize the returned list.
    tfpn_list = []
    
    # Here we iterate through it, pulling each tab-
    # delimited line as a list into our list.
    for line in csv.reader(tfpn_bedfile, delimiter='\n'):
        for coordinate in csv.reader(line, delimiter='\t'):
            tfpn_list.append(coordinate)
    # Here we check and make sure that there are
    # values to work with here.
    if len(tfpn_list) == 0:
        raise IOError('The file at '+file+' could not be read.')
    # Closing files is hygienic.
    tfpn_bedfile.close()

    return tfpn_list


def find_file_length(file):
    
    # This is a function that reads a file, split it
    # by line and returns the length of the file.
    
    # Here we open the file.
    file_to_check = open(file,'rU')
    
    # Here we find the file length.
    length = len(file_to_check.readlines())

    # Closing files is hygienic.
    file_to_check.close()

    # We return the file length.
    return length


def generate_lognormal():

    # This is a very simple function that generates
    # a number in a lognormal distribution that just
    # has to be above a certain limit.

    # Here we initialize our returned list.
    lognorm_value = np.random.lognormal()
    
    # Here we ensure the value results in a p-value
    # that is high enough to pass the limit. If not,
    # we reroll the number.
    while (math.pow(10.0,-lognorm_value) < math.pow(10.0,-5)):
        lognorm_value = np.random.lognormal()

    # We return the value.
    return math.pow(10.0,-lognorm_value)


# This function generates the values for true and false
# positives and true and false negatives.
def produce_pred_values(directory):

    # Here we find the number of PWMs in each category.
    tp_len = find_file_length(directory+'/true_positives.bed')
    fp_len = find_file_length(directory+'/false_positives.bed')
    fn_len = find_file_length(directory+'/false_negatives.bed')
    tn_len = find_file_length(directory+'/true_negatives.bed')

    # We have to hold onto the p-values given by the footprint
    # algorithm for our positive assertions.
    tp_list = read_tfpn_bedfile(directory+'/true_positives.bed')
    fp_list = read_tfpn_bedfile(directory+'/false_positives.bed')

    # We initialize vectors to return.
    tp_values = np.zeros(tp_len)
    fp_values = np.zeros(fp_len)
    fn_values = np.zeros(fn_len)
    tn_values = np.zeros(tn_len)

    # We pass through the categories, altering the values to
    # match either the score from the footprinting algorithm
    # or a random lognormal p-value.
    for i in range(len(tp_values)):
        tp_values[i] = 1.0 - math.pow(10.0,-1*abs(float(tp_list[i][14])))
    for i in range(len(fp_values)):
        fp_values[i] = 1.0 - math.pow(10.0,-1*abs(float(fp_list[i][14])))

    for i in range(len(fn_values)):
        fn_values[i] = generate_lognormal()
    for i in range(len(tn_values)):
        tn_values[i] = generate_lognormal()

    # We return the vectors for each category.
    return tp_values,fp_values,fn_values,tn_values


# Here we turn the categories into prediction vs. labels
# vectors.
def craft_predlabel_vectors(tpv,fpv,fnv,tnv):

    # We initialize our starting arrays.
    predictions = np.zeros(len(tpv)+len(tnv)+len(fpv)+len(fnv))
    labels = np.zeros(len(tpv)+len(tnv)+len(fpv)+len(fnv))

    # We run through the list of categories,
    # assigning values based on the passed p-values.
    for j in range(len(tpv)):
        predictions[j] = tpv[j]
        labels[j] = 1
    for j in range(len(fpv)):
        predictions[j+len(tpv)] = fpv[j]
        labels[j+len(tpv)] = 0
    for j in range(len(fnv)):
        predictions[j+len(tpv)+len(fpv)] = fnv[j]
        labels[j+len(tpv)+len(fpv)] = 1
    for j in range(len(tnv)):
        predictions[j+len(tpv)+len(fpv)+len(fnv)] = tnv[j]
        labels[j+len(tpv)+len(fpv)+len(fnv)] = 0

    # We return the predictions and labels.
    return predictions, labels


# The final task is to write the AUC out to a file.
def iterate_through_chip():
    
    mean_tpr = []
    mean_fpr = np.linspace(0,1,100)
    j = 0
    
    for chip_target in os.listdir('ChIP-validation'):
        tpv,fpv,fnv,tnv = produce_pred_values('ChIP-validation/'
                                              +chip_target)
        if len(tpv)+len(fpv)+len(fnv) > 0:
            preds, labels = craft_predlabel_vectors(tpv,fpv,fnv,tnv)
            fpr,tpr,thresholds = skm.roc_curve(labels,preds)
            new_mean = np.interp(mean_fpr,fpr,tpr).tolist()
            if j == 0: mean_tpr += new_mean
            else:
                for point in range(len(mean_tpr)):
                    mean_tpr[point] += new_mean[point]
            mean_tpr[0] = 0.0
            roc_auc = skm.auc(fpr, tpr)
            print chip_target,"AUC:",roc_auc
            j += 1

    for point in range(len(mean_tpr)):
        mean_tpr[point] /= j
    mean_tpr[-1] = 1.0
    mean_auc = skm.auc(mean_fpr,mean_tpr)
    print "Mean AUC:",mean_auc


def __main__():

    # This function governs all other functions and runs
    # only when the program is run directly.

    # Here we read the .bed file input.
    iterate_through_chip()


# The condition that the program is run directly.
# This makes unit testing easy.
if __name__ == '__main__':
    __main__()