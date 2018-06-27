import matplotlib
matplotlib.use('Agg')

import numpy as np
from sklearn import metrics as skm
import random
import math

import argparse

import os
import csv
from glob import glob

import matplotlib.pyplot as plt
import seaborn as sns

# It's important to set a seed so that results are
# reproducible.
random.seed(10000)

# Using the files that are in JUBAL_figures, the
# folders can be linked to specific colors and
# labels which make plotting easier later.
def make_color_dictionary(c):

    filepath = '/data/pranzatellitj/JUBAL_figures/i-'+c+'.txt'
    readfile = open(filepath).read()
    ideals = []; names = {}; colors  = {}
    for line in readfile.split('\n'):
        tabs = line.split('\t')
        ideals.append(tabs[0])
        names[tabs[0]] = tabs[1]
        colors[tabs[0]] = tabs[2]
    return ideals,names,colors

# This function reads the ensembl gene space
# and translates to a the slightly more human-
# readable alphabet soup of old school bio.
def read_ensg_name():

    # Initialize the dictionary.
    ensg_name = {}
    # Here's where the file is.
    filepath = '/data/ChioriniCompCor/experiments/2017-09-01_Metamachine_Footprinting_Pipelines/genomes/hg19/name_to_ensg.txt'
    # We read the file into memory.
    openfile = open(filepath).read()
    # The file is split on lines.
    lines = openfile.split('\n')
    # At each line, the name is the value
    # and the key is the ENSG.
    for line in lines[:-1]:
        tab = line.split('\t')
        ensg_name[tab[1]] = tab[0]
    # The dictionary is returned.
    return ensg_name


def read_tfpn_bedfile(file):
    
    # This is a function that reads a .bed file
    # of the true/false positives/negatives and
    # relays that to a list object.
    
    # Here we open the .bed file.
    tfpn_bedfile = open(file,'rU')
    
    # Here we initialize the returned list.
    tfpn_list = []
    
    # Here we iterate through it, pulling each tab-
    # delimited line as a list into our list.
    for line in csv.reader(tfpn_bedfile, delimiter='\n'):
        for coordinate in csv.reader(line, delimiter='\t'):
            tfpn_list.append(coordinate)

    # Closing files is hygienic.
    tfpn_bedfile.close()

    return tfpn_list


def find_file_length(file):
    
    # This is a function that takes a .bed file
    # and returns the length of it.
    
    # Here we open the .bed file.
    file_to_check = open(file,'rU')
    
    # Here we determine the length.
    length = len(file_to_check.readlines())

    # Closing files is hygienic.
    file_to_check.close()
    
    return length


def generate_lognormal(limit):

    # This is a function that takes a limit
    # which, in this case, is the p-value
    # cutoff of Wellington, and generates a
    # value on a lognormal distribution. An
    # alternative is to use an exponential.

    # Here we initialize the value.
    lognorm_value = np.random.lognormal()
    
    # Here we iterate until we get a value
    # that falls within the cutoff.
    while (math.pow(10.0,-lognorm_value) < math.pow(10.0,limit)):
        lognorm_value = np.random.lognormal()

    return math.pow(10.0,-lognorm_value)


def produce_pred_values(directory, limit):

    # This is a function that returns prediction
    # values for a directory, hampered by a limit.
    
    # Here we find the lengths of each kind of
    # result (true positive,..., false negative).
    tp_len = find_file_length(directory+'/true_positives.bed')
    fp_len = find_file_length(directory+'/false_positives.bed')
    fn_len = find_file_length(directory+'/false_negatives.bed')
    tn_len = find_file_length(directory+'/true_negatives.bed')

    # The positives have a value, which is the
    # p-value of the predictive footprint.
    tp_list = read_tfpn_bedfile(directory+'/true_positives.bed')
    fp_list = read_tfpn_bedfile(directory+'/false_positives.bed')

    # Here the list of values is initialized.
    tp_values = np.zeros(tp_len)
    fp_values = np.zeros(fp_len)
    fn_values = np.zeros(fn_len)
    tn_values = np.zeros(tn_len)

    # Here the positives are assigned their
    # p-values.
    for i in range(len(tp_values)):
        tp_values[i] = 1.0 - math.pow(10.0,-1*abs(float(tp_list[i][14])))
    for i in range(len(fp_values)):
        fp_values[i] = 1.0 - math.pow(10.0,-1*abs(float(fp_list[i][14])))

    # Here the negatives are given a
    # value generated along a distribution.
    for i in range(len(fn_values)):
        fn_values[i] = generate_lognormal(limit)
    for i in range(len(tn_values)):
        tn_values[i] = generate_lognormal(limit)

    return tp_values,fp_values,fn_values,tn_values


def craft_predlabel_vectors(tpv,fpv,fnv,tnv):
    
    # This is a function that generates two long
    # lists, one of predictions and one of labels.

    # Here the predictions and labels are initialized.
    predictions = np.zeros(len(tpv)+len(tnv)+len(fpv)+len(fnv))
    labels = np.zeros(len(tpv)+len(tnv)+len(fpv)+len(fnv))

    # In each of these for loops, a prediction is
    # added to the master list of predictions, and
    # a label is added to the master list of labels.
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

    return predictions, labels

# This function plots the AUC for all TFs for each
# pipeline. The X axis is each TF and the Y axis is
# the AUC.
def Barozzi_plot(experiment,color,name,axis,ensg_name):

    filepath = ('/data/pranzatellitj/experiments/2016-10-17-JUBAL_Cotillion_DNase_Grid_Search/'+experiment+'/AUC.txt')
    # The file is read into memory.
    openfile = open(filepath).read()
    # Initialize the x and y vectors.
    x = []; y = []
    # The file is split by line and for each TF...
    for line in openfile.split('\n')[:-2]:
        # Each line is TF to AUC value.
        tabs = line.split(' AUC: ')
        # The x and y vectors are populated.
        x.append(ensg_name[tabs[0]]); y.append(float(tabs[1]))
    # The vectors are sorted by X axis. Easier to plot.
    x,y = (list(k) for k in zip(*sorted(zip(x,y))))
    # There are x ticks and they're 0-something.
    ticks = range(len(x))
    # We plot the y value on the x axis.
    axis.plot(ticks,y,color=color,label=name)
    # And, we get to name the ticks by TF.
    plt.xticks(ticks,x,rotation=45)

# This plot is of the overlap as the allowed overlap
# scales.
def span_plot(experiment,color,name,axis):

    filepath = ('/data/pranzatellitj/experiments/2016-10-17-JUBAL_Cotillion_DNase_Grid_Search/'+experiment+'/overlap.txt')
    # The file is read into memory.
    openfile = open(filepath).read()
    # Vectors are initialized.
    x = []; y = []
    # Lines are iterated through.
    for line in openfile.split('\n')[:-1]:
        # Vectors are populated.
        spaces = line.split(' ')
        x.append(int(spaces[0])); y.append(float(spaces[1]))
    # This is simply plotted as a line.
    axis.plot(x,y,color=color,label=name)


def iterate_through_chip(experiment,color,name,axis):
    
    mean_tpr = []
    mean_fpr = np.linspace(0,1,100)
    limit = -5
    i = 0
    for chip_target in os.listdir('/data/pranzatellitj/experiments/2016-10-17-JUBAL_Cotillion_DNase_Grid_Search/'+experiment+'/ChIP-validation'):
        tpv,fpv,fnv,tnv = produce_pred_values('/data/pranzatellitj/experiments/2016-10-17-JUBAL_Cotillion_DNase_Grid_Search/'+experiment+'/ChIP-validation/'+chip_target, limit)
        if len(tpv)+len(fpv)+len(fnv) > 0:
            preds, labels = craft_predlabel_vectors(tpv,fpv,fnv,tnv)
            fpr,tpr,thresholds = skm.roc_curve(labels,preds)
            new_mean = np.interp(mean_fpr,fpr,tpr).tolist()
            if i == 0: mean_tpr += new_mean
            else:
                for point in range(len(mean_tpr)):
                    mean_tpr[point] += new_mean[point]
            mean_tpr[0] = 0.0
            roc_auc = skm.auc(fpr, tpr)
            i += 1

    for point in range(len(mean_tpr)):
        mean_tpr[point] /= i
    mean_tpr[-1] = 1.0
    mean_auc = skm.auc(mean_fpr,mean_tpr)
    plt.plot(mean_fpr, mean_tpr, color=color, label=name)
    plt.plot([0,1],[0,1],'--',color=(.7,.7,.7))
    plt.xlim([-0.05,1.05])
    plt.ylim([-0.05,1.05])


def __main__():

    # This function governs all other functions and runs
    # only when the program is run directly.

    # This function differs from the standard functional
    # module style of coding because it makes plotting and
    # editing easier.

    ensg_name = read_ensg_name()
    
    ideals,names,colors = make_color_dictionary('d')

    plt.clf(); plt.cla()
    f, ax = plt.subplots(1,1)
    for ideal in ideals:
        iterate_through_chip(ideal,colors[ideal],names[ideal],ax)
        ax.legend(loc='center left',bbox_to_anchor=(1,0.5),fancybox=True)
        plt.xlabel('Sensitivity')
        plt.ylabel('Specificity')
        plt.title('ROC Curves for ChIP Validation')
        plt.savefig('/data/pranzatellitj/JUBAL_figures/cotillion_fig1_DNase.pdf',format='pdf',dpi=600,bbox_inches='tight')
    plt.clf(); plt.cla()
    f, ax = plt.subplots(1,1)
    for ideal in ideals:
        Barozzi_plot(ideal,colors[ideal],names[ideal],ax,ensg_name)
        ax.legend(loc='center left',bbox_to_anchor=(1,0.5),fancybox=True)
        plt.xlabel('Transcription Factors')
        plt.ylabel('Area Under Curve')
        plt.title('Per-TF AUC Values')
        plt.savefig('/data/pranzatellitj/JUBAL_figures/cotillion_fig2_DNase.pdf',format='pdf',dpi=600,bbox_inches='tight')
    plt.clf(); plt.cla()
    f, ax = plt.subplots(1,1)
    for ideal in ideals:
        span_plot(ideal,colors[ideal],names[ideal],ax)
        ax.legend(loc='center left',bbox_to_anchor=(1,0.5),fancybox=True)
        plt.xlabel('Percentage Overlap Allowed')
        plt.ylabel('Percentage of Footprints Overlapping')
        plt.title('Footprint Overlap Between Replicates')
        plt.savefig('/data/pranzatellitj/JUBAL_figures/cotillion_fig3_DNase.pdf',format='pdf',dpi=600,bbox_inches='tight')


    ideals,names,colors = make_color_dictionary('a')

    plt.clf(); plt.cla()
    f, ax = plt.subplots(1,1)
    for ideal in ideals:
        iterate_through_chip(ideal,colors[ideal],names[ideal],ax)
        ax.legend(loc='center left',bbox_to_anchor=(1,0.5),fancybox=True)
        plt.xlabel('Sensitivity')
        plt.ylabel('Specificity')
        plt.title('ROC Curves for ChIP Validation')
        plt.savefig('/data/pranzatellitj/JUBAL_figures/cotillion_fig1_ATAC.pdf',format='pdf',dpi=600,bbox_inches='tight')
    plt.clf(); plt.cla()
    f, ax = plt.subplots(1,1)
    for ideal in ideals:
        Barozzi_plot(ideal,colors[ideal],names[ideal],ax,ensg_name)
        ax.legend(loc='center left',bbox_to_anchor=(1,0.5),fancybox=True)
        plt.xlabel('Transcription Factors')
        plt.ylabel('Area Under Curve')
        plt.title('Per-TF AUC Values')
        plt.savefig('/data/pranzatellitj/JUBAL_figures/cotillion_fig2_ATAC.pdf',format='pdf',dpi=600,bbox_inches='tight')
    plt.clf(); plt.cla()
    f, ax = plt.subplots(1,1)
    for ideal in ideals:
        span_plot(ideal,colors[ideal],names[ideal],ax)
        ax.legend(loc='center left',bbox_to_anchor=(1,0.5),fancybox=True)
        plt.xlabel('Percentage Overlap Allowed')
        plt.ylabel('Percentage of Footprints Overlapping')
        plt.title('Footprint Overlap Between Replicates')
        plt.savefig('/data/pranzatellitj/JUBAL_figures/cotillion_fig3_ATAC.pdf',format='pdf',dpi=600,bbox_inches='tight')


# The condition that the program is run directly.
# This makes unit testing easy.
if __name__ == '__main__':
    __main__()