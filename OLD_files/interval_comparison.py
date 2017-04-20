##python script - Genes within interval

from __future__ import division
import csv
import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--input',action='store', required=True, metavar='sample.bed',  help="Input list of intervals - Bed3 format")
parser.add_argument('--ref',action='store', required=True, metavar='ref.bed', help="Reference file of gene intervals - Bed4 format")
args = parser.parse_args()


def inputfunc():
	with open(args.input, 'rb') as i:
		infile = pd.read_table(i, sep='\t', header = None)
		infile.columns = ['Chr', 'Start', 'Stop', 'Sample']
		unichr = pd.unique(infile.Chr.ravel())
		#print(unichr)
		#print(infile)
		#Reads in inputfile in bed format and defines the col headers - setting the raw file as infile - generates a list of unique chromosome positions
	with open(args.ref, 'rb') as r:
		reffile = pd.read_table(r, sep='\t', header = None, low_memory = False)
		reffile.columns = ['Chr', 'Start', 'Stop', 'Gene']
		#print(reffile)
		#Reads in reference list of intervals for known genes - bed format and sets the local variable to reffile
	out = [[]]
	#output list which is returned to be appended to output dataframe
	for v in unichr:
		Chr = v
		inparse = (infile.loc[infile['Chr'] == v])
		#for loop for each item in unichar list - rows matching v (chr) are parsed as inparse to the next loop
		refparse = (reffile.loc[reffile['Chr'] == v])
		#same as previous comment but applied to the ref file
		#print(inparse)
		#print(refparse)
		for index, r in inparse.iterrows():	
			istart = r['Start']
			istop = r['Stop']
			samp = r['Sample']
			#print(r)
			#For loop for each row in each chr, define the start and stop points of the interval in the input file
			for index, row in refparse.iterrows():
				gstart = row['Start']
				gstop = row['Stop']
				info = row['Gene']
				add = calcoverlap(istart,istop,gstart,gstop)
				#print(add)
				#for loop for each row in the ref file and define start and stop - assess if overlap is true using calcoverlap function
				if add  == True:
					out.append([Chr, istart, istop, samp, Chr, gstart, gstop, info])
					#print('Added to Dataframe')
					#if calcoverlap returns true add value to out list
				else:
					#print('Not found')
					break
					#if calcoverlap returns false break loop and try next ref gene
	return out #out list returned to be appended to empty dataframe 



def calcoverlap(istart,istop,gstart,gstop):
	#series of logical comparisons of values derived from the interval positions of both input and ref
	#gene contained within interval
	if int(gstart) >= int(istart) and int(gstop) <= int(istop) and int(gstart) <= int(istop) and int(gstop) >= int(istart):
		return True
		print(true)
	#gene overlaps entire interval
	elif int(gstart) <= int(istart) and int(gstop) >= int(istop) and int(gstart) <= int(istop) and int(gstop) >= int(istart):
		return True
		print(true)
	#gene overlaps left side of interval
	elif int(gstart) <= int(istart) and int(gstop) <= int(istop) and int(gstart) <= int(istop) and int(gstop) >= int(istart):
		return True
		print(true)
	#gene overlaps right side of interval
	elif int(gstart) >= int(istart) and int(gstop) >= int(istop) and int(gstart) <= int(istop) and int(gstop) >= int(istart):
		return True
		print(true)
	else:
		return False
		print(false)

def format(write):
	wr2 = write.ix[1:]
	wr2.columns = cols
	wr2.to_csv('merged_interval.txt', sep='\t')



#core function calls and header assignments - varible that holds final dataframe after analysis
cols = ['chr', 'start', 'end', 'info', 'chr', 'start', 'end', 'info']
o = pd.DataFrame()
write = o.append(inputfunc())
format(write)


