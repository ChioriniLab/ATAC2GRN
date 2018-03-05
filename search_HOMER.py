# Thomas Pranzatelli
# This code reads the results of running the pipelines with different
# HOMER peak sizes and measuring the number of footprints produced
# and the biological reproducibility of the footprints.

# This odd block is where we plot on the cluster. We can also plot
# in 3D.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# This is the main import block for most of our data analysis tools.
import csv
import os
from glob import glob
import numpy as numpy
import pandas as pd
import seaborn as sns

# We can perform systems testing as well.
import unittest


# This is a function that reads in a file and returns
# a list as though it were tab-delimited (which it
# hopefully is).
def read_file(file):
    # The file must first be opened to RAM.
    openfile = open(file)
    # We initialize the return list.
    return_list = []
    # First we split the file on its lines...
    for line in csv.reader(openfile,delimiter='\n'):
        # Then we split the lines by their tabs...
        for tab in csv.reader(line,delimiter='\t'):
            # Finally we form a list.
            return_list.append(tab)
    # If the file is empty, we complain.
    if len(return_list) == 0:
        raise IOError('The file '+file+' does not have multiple lines to read.')
    # We return this final list.
    return return_list

# This is just a function that finds all folders that match
# the sequence passed to glob and returns the list.
def read_folders(experiment):
    # So this will be either homer_DNase/ or homer_ATAC/ folders.
    folders = glob('/data/pranzatellitj/check_HOMER_stitching/homer_'+experiment+'/*')
    # We make sure that multiple folders are returned.
    if len(folders) == 0:
        raise IOError('The folders for '+experiment+' do not exist in this directory.')
    # We return this list of folders.
    return folders

# This function gets the AUC values from a folder and adds
# them to a dictionary.
def get_auc(folder,dictionary):
    # We read the file and get the tab-delimited list.
    auclist = read_file(folder+'/AUC.txt')
    # Making sure it isn't empty...
    if auclist != []:
        # We collect the AUC.
        auc = auclist[-1][0].split(': ')[1]
        # We then update the dictionary.
        dictionary['Mean AUC'] = float(auc)
    else:
        # We need to notify the code-runner that something
        # is awry.
        raise RuntimeError('The AUC.txt file at '+folder+' is empty.')
    # We finally return the updated dictionary.
    return dictionary

# We need to generate each series to form the DataFrame.
# The series are populated with information found in the
# folder for each argument.
def generate_series(folder,experiment):
    try:
        # We initialize the series.
        dict_to_Series = {}
        # We pass either DNase or ATAC to hold onto
        # that information.
        dict_to_Series['Footprinting'] = experiment
        # The name of the folder is the size of the homer
        # peaks that we're passing to HOMER. We need to
        # populate the dictionary with this info as well.
        size = int(folder.split('/')[-1])
        dict_to_Series['Peak Size'] = size
        # We finally tell the dictionary to add AUC 
        # information for the y axis.
        dict_to_Series = get_auc(folder,dict_to_Series)
    except:
        # If the file is empty, we complain.
        raise RuntimeError('Series could not be generated from '+folder+' .')

# The DataFrame is produced from Series objects.
def generate_DF():
    # We initialize the DataFrame dictionary and the index.
    seriesdict= {}; index = 0
    # We do this for both ATAC and DNase...
    for experiment in ['ATAC','DNase']:
        # We read the folder list into memory...
        folders = read_folders(experiment)
        # We iterate through the folders...
        for folder in folders:
            index += 1
            # Generating series of data...
            series = generate_series(folder,experiment)
            # And appending them to the dictionary.
            appendable = pd.Series(dict_to_Series)
            seriesdict[index] = series
    # Final step is to convert to DataFrame and return.
    DF = pd.DataFrame(seriesdict).transpose().fillna(value=0)
    return DF

# It's important to plot our results!
def plot(DF):
    # We'll use seaborn's pointplot, where X is the peak size
    # from the folder name and Y is the mean AUC. We'll
    # separate the lines by DNase/ATAC.
    sns.pointplot(x="Peak Size",y="Mean AUC",hue="Footprinting",palette="muted",data=DF)
    # And, of course, we need to save the figure we produce.
    plt.savefig('/data/pranzatellitj/JUBAL_figures/cotil_fig_HomerAUC.pdf',format='pdf',dpi=600)

def __main__():
    DF = generate_DF():
    plot(DF)

if __name__ == '__main__':
    __main__()
    unittest.main()