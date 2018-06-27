# Thomas Pranzatelli
# This code reads in the results of running the pipelines with different
# read depths and measuring the number of footprints produced and the
# biological reproducibility of the footprints produced.

# This is a weird block that allows us to plot from the cluster.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Here is where we import our modules of interest.
import csv
from os import system
from glob import glob
import numpy as numpy
import pandas as pd
import seaborn as sns

# We can also perform unit and systems testing.
import unittest


# This is just a function that finds all folders that match
# the sequence passed to glob and returns the list.
def read_folders(experiment):
    # So this will be either experiments_DNase/ or experiments_ATAC/ folders.
	folders = glob('/data/pranzatellitj/check_HOMER_stitching/experiments_'+experiment+'/*')
	# We make sure that multiple folders are returned.
    if len(folders) == 0:
        raise IOError('No folders were discovered for '+experiment+' .')
    # We return this list of folders.
	return folders

# This is the generic file-reading function. Files should be
# in a tab-delimited format.
def read_file(file):
	# The first step is to open the file to RAM.
	openfile = open(file)
	# We initialize the list object.
	return_list = []
	# We iterate through the file, split by line.
	for line in csv.reader(openfile,delimiter='\n'):
		# The lines are split by tab.
		for tab in csv.reader(line,delimiter='\t'):
			# The tab-split lines are added to the return list.
			return_list.append(tab)
	# If the file is empty, we complain.
    if len(return_list) == 0:
        raise IOError('The file '+file+' does not have multiple lines to read.')
	# Of course, that list is returned.
	return return_list

# This function gets the AUC values for each TF from the
# AUC.txt file. The one we care about is the Mean AUC.
def get_AUC(folder,dictionary):
	# We read in the file to our tab-split standard.
	auclist = read_file(folder+'/AUC.txt')
	# Assuming it isn't empty...
	if auclist != []:
		# We iterate through it...
		for auc in auclist:
			# We populate the dictionary with each TF and its
			# associated AUC.
			key = auc[0].split(': ')[0]
			value = float(auc[0].split(': ')[1])
			dictionary[key] = value
	else:
		# This is a problem that requires an error.
		raise RuntimeError('The AUC.txt file at '+folder+' is empty.')
	# We return the dictionary.
	return dictionary

# Determine the number of reads and the round of a project from
# the folder name. As an example, if a folder in named
# 160000000-3, it's the third random sampling of sixteen
# million reads.
def get_reads_rounds(folder):
	# The reads are the first part of the folder name and are
	# divided by a million for a number in millions.
	reads = int(folder.split('/')[-1].split('-')[0])/1000000
	# The round is the second component of the folder name.
	round = int(folder.split('/')[-1].split('-')[1])
	# These values are returned.
	return reads,round

# This function gets the number of footprints produced at each
# read depth.
def get_footprint_num(folder,dictionary):
	# The file is read to tab-split format.
	fplist = read_file(folder+'/footprints.bed')
	# Assuming the file isn't empty...
	if auclist != []:
		# The number of footprints is added to the dictionary.
		dictionary['Number of Footprints'] = len(fplist) - 1
	else:
		# We need to inform that the footprints file doesn't
		# exist.
		raise RuntimeError('No footprints file found in '+folder+'.')
	# Finally, the dictionary is returned.
	return dictionary

# This function generates the series that will populate
# the dataframe.
def generate_Series(experiment,reads,round):
	# The dictionary is initialized.
	dict_to_Series = {}
	# The first entries are filled based on given information.
	dict_to_Series['Footprinting'] = experiment
	dict_to_Series['Reads'] = reads
	dict_to_Series['Round'] = round
	try:
		# The dictionary is updated with AUC and footprint
		# quantification information.
		dict_to_Series = get_AUC(folder,dict_to_Series)
		dict_to_Series = get_footprint_num(folder,dict_to_Series)
		# This dictionary is transformed to a pandas Series
		# and returned to form a dataframe.
		appendable = pd.Series(dict_to_Series)
		return appendable
	except:
		# We raise an error.
		raise RuntimeError('The Series could not be populated.')

# This is the main function, which generates the dataframe.
def generate_DF()
	# The dictionary of series that will form the dataframe
	# is initalized, as is the index that will serve as
	# keys.
	seriesdict= {}; index = 0
	# This is done for both ATAC and DNase.
	for experiment in ['ATAC','DNase']:
		# Folders are read for each experimental condition.
		folders = read_folders()
		# For each project folder...
		for folder in folders:
			# The index of course has to change.
			index += 1
			# The reads and rounds are calculated from the
			# folder name.
			reads,round = get_reads_rounds(folder)
			# If these are even (to avoid clutter)...
			if reads % 20 == 0:
				# The series is generated.
				series = generate_series(experiment,reads,rounds)
				# Assuming one exists...
				if series:
					# The dictionary is populated!
					seriesdict[index] = series
	# Finally, the dataframe is transformed and returned.
	DF = pd.DataFrame(seriesdict).transpose().fillna(value=0)
	return DF

# These functions plot the AUC and footprints against the
# read depth. Both of these functions save to data using savefig,
# so that this task can be done on the cluster.
def plot_1(DF):
	# Here we clear the figure to prevent confusing images.
	plt.clf()
	# A "pointplot" is what seaborn calls a line. We're using this
	# to show some trend as reads increase.
	sns.pointplot(x="Reads",y="Mean AUC",hue="Footprinting",palette="muted",data=DF)
	# This is a high-dpi pdf image, suitable for publication.
	plt.savefig('/data/pranzatellitj/JUBAL_figures/cotil_fig_readsAUC.pdf',format='pdf',dpi=600)

def plot_2(DF)
	# The figure is again cleared.
	plt.clf()
	# We produce a line plot of footprints over reads.
	sns.pointplot(x="Reads",y="Number of Footprints",hue="Footprinting",palette="muted",data=DF)
	# We save the figure at high resolution.
	plt.savefig('/data/pranzatellitj/JUBAL_figures/cotil_fig_readsFP#.pdf',format='pdf',dpi=600)

def __main__():
	DF = generate_DF()
	plot1(DF)
	plot2(DF)


# We're very limited in the things we can unit test for this
# code because everything is I/O. This project predominantly
# uses small test files and test cases. For unit testing,
# we'll just test the ability of the program to read a folder
# name and get the reads and round from it.
class TestCase_DownsamplingReads(unittest.TestCase):
	def test_readround(self):
		foldername = '10/70-2'
		reads,round = get_reads_rounds(foldername)
		self.assertEqual(reads,70)
		self.assertEqual(round,2)
		

if __name__ == '__main__':
	__main__()
	unittest.main()

