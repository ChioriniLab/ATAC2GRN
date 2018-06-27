# Thomas Pranzatelli
# This code plots the different trimmed reads and plots their alignment
# to the genome using standard out-of-the-box Bowtie2 alongside the
# expected proportion of reads aligning stochastically.

# This is a weird method that allows us to plot from a cluster.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# This is the import block.
from glob import glob
from math import pow
import numpy as np

# We can perform systems testing.
import unittest


# First we grab all of the slurm outputs on the cluster in this
# project folder. We'll be reading the percentage of reads that align
# from these files.
def read_cluster_files():
	# This is the list of all files that begin with "slurm-".
	files = glob('slurm-*')
	# Here we check and make sure that we actually have files.
	if len(files) == 0:
		raise IOError('No cluster files beginning with "slurm-*" were found.')
	# Otherwise we return what we have.
	else:
		return files

# Opens and reads each of the files to a list of the lines. We'll
# only end up using the seventh and seventeenth lines.
def read_files_to_lines(files):
	# We initialize the all_files list object.
	all_lines = []
	# We iterate through the list of files we've passed...
	for file in files:
		try:
			# ...opening them to RAM...
			openfile = open(file)
			# ...reading them to a string in RAM...
			readfile = openfile.read()
			# ...closing the opened file object...
			openfile.close()
			# ...splitting the RAM string by lines ('\n') and then
			# appending it to the output list.
			if len(readfile.split('\n')) < 17:
				raise RuntimeError('File '+file+' has fewer lines than expected.')
			all_lines.append(readfile.split('\n'))
		except:
			raise IOError('Somehow the file '+file+' has been lost.')
	# And then we return the all_files object!
	return all_lines

# This function reads the files at the seventh and seventeenth
# lines to produce data.
def turn_lines_to_data(all_lines):
	# We initialize the lists that store the data.
	x = []; y1 = []; y2 = []
	# We iterate through the files, taking only those files
	# that correspond to instances in which the files had 10 bases
	# trimmed from the front. These were a representative set, but
	# in the paper we end up using a HEADCROP of 20 to remove what
	# looks in fastqc like tagmentation+overhang.
	for lines in all_lines:
		if int(lines[6].split(' ')[0]) == 10:
			# The final length of the reads. This corresponds to
			# the CROP argument.
			length = int(lines[6].split(' ')[1])
			# We append the final length of the read to X.
			x.append(length)
			# The percentage alignment to the genome is stored
			# where bowtie2 prints it out. We'd like to
			# maximize this number.
			y1.append(float(lines[16].split('%')[0])/100)
			# Human genome is 3.2 billion base pairs. That number
			# divided by 4 to the power of the length of the read
			# is a conservative estimate of the chance each read
			# will be randomly mapped to the genome. We'd like
			# to minimize this number.
			y2.append(3234000000/pow(4,length))
	# Finally, we return the data.
	return x,y1,y2

# Here we sort the data so that we have a clean plot. We'd like
# to sort by the x axis, and have the corresponding values in y1
# and y2 sort with x.
def sort_data(x,y1,y2):
	# Zip makes sorting very easy by linking the values at
	# certain indices in each list.
	zipped = zip(x,y1,y2)
	# Then we have the cleverly named sort().
	zipped.sort()
	# We unzip the lists so that we can use them individually.
	x,y1,y2 = zip(*zipped)
	# Finally, they are returned.
	return x,y1,y2

# This is the final function, and in it we plot the data and see
# what final read lengths have an acceptably low chance of randomly
# aligning to the genome and still align at some rate.
def plot_data(x,y1,y2):
	# We form axes for two plots which will end up as one graph.
	fig,ax1 = plt.subplots()
	ax2 = ax1.twinx()
	# We plot the proportion of alignment to the genome in green.
	ax2.plot(x,y1,color='#1b9e77')
	ax2.set_ylabel('Proportion of Alignment to the Genome',color='#1b9e77')
	ax2.set_ylim([0,1])
	# We plot the proportion of random alignment in orange.
	ax1.plot(x,y2,color='#d95f02')
	ax1.set_ylabel('Expected Proportion of Random Alignment',color='#d95f02')
	ax1.set_ylim([0,1])
	# We set the x-ticks to correspond to the X-axis. This is why
	# sorting with links is so important.
	ax1.set_xticks(x)
	plt.xlim([15,25])
	ax1.set_xlabel('Final Length of Read')
	# Finally, we save the graph on the cluster.
	plt.savefig('GM12878_Trimming_Line.pdf',format='pdf',dpi=600,bbox_inches='tight')


def __main__():
	files = read_cluster_files()
	all_lines = read_files_to_lines(files)
	x,y1,y2 = turn_lines_to_data(all_lines)
	x,y1,y2 = sort_data(x,y1,y2)
	plot_data(x,y1,y2)


# It's important to run your code through tests. This code is
# difficult to test, because the units involve integrations with
# an environment and outside objects that may or may not function
# as intended. So, most of the code will be error code, and a few
# unit tests will be below.
class TestCase_TrimPlot(unittest.TestCase):
	def test_lines_to_data(self):
		all_lines = [[0,1,2,3,4,5,'10 18',7,8,9,10,11,12,13,14,15,'39.8% *']]
		x,y1,y2 = turn_lines_to_data(all_lines)
		self.assertEqual(x[0],18)
		self.assertAlmostEqual(y1[0],0.398)
	def test_lines_only_ten(self):
		all_lines = [[0,1,2,3,4,5,'9 18',7,8,9,10,11,12,13,14,15,'39.8% *']]
		x,y1,y2 = turn_lines_to_data(all_lines)
		self.assertEqual(len(x),0)
	def test_lines_float(self):
		all_lines = [[0,1,2,3,4,5,'10 18',7,8,9,10,11,12,13,14,15,'39.8% *']]
		x,y1,y2 = turn_lines_to_data(all_lines)
		self.assertTrue(type(y1[0]) == type(1.0))
		self.assertTrue(type(y2[0]) == type(1.0))
	def test_lines_equal_size(self):
		all_lines = [[0,1,2,3,4,5,'9 18',7,8,9,10,11,12,13,14,15,'39.8% *']]
		x,y1,y2 = turn_lines_to_data(all_lines)
		self.assertTrue(len(x) == len(y1) == len(y2))

	def test_sort_x(self):
		x = [1,0]; y1 = ['B','A']; y2 = [1.0,0.0]
		out_x,out_y1,out_y2 = sort_data(x,y1,y2)
		self.assertEqual(out_x[0],0)
		self.assertEqual(out_x[1],1)
	def test_sort_y1(self):
		x = [1,0]; y1 = ['B','A']; y2 = [1.0,0.0]
		out_x,out_y1,out_y2 = sort_data(x,y1,y2)
		self.assertSequenceEqual(out_y1[0],'A')
		self.assertSequenceEqual(out_y1[1],'B')
	def test_sort_y2(self):
		x = [1,0]; y1 = ['B','A']; y2 = [1.0,0.0]
		out_x,out_y1,out_y2 = sort_data(x,y1,y2)
		self.assertAlmostEqual(out_y2[0],0.0)
		self.assertAlmostEqual(out_y2[1],1.0)



if __name__ == '__main__':
	__main__()
	unittest.main()
