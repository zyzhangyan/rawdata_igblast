#!/usr/bin/env python
# encoding: utf-8
"""
zhangyantools.py

Created by Zhang yan on 2016-06-12.
Copyright (c) 2016 Southern Medical University. All rights reserved.
"""
import sys, os, csv, re,collections,math
from collections import Counter
#from common_info import *
from Bio import Seq
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalwCommandline
from Bio.Blast import NCBIStandalone
from math import log
from itertools import imap
from numpy import mean, array, zeros, ones, nan, std, isnan
import getopt,sys 
#
#--	BEGUN -- Count length
#
def count_length(infile):
	fq = open(infile,"r")
	len_distribution = csv.writer(open("%s_length_distribution.txt"%(infile), "wb+"),delimiter="\t")
	seq_len = []
	for record in SeqIO.parse(fq,"fastq"):
		length = len(record)
		seq_len.append(length)
	group_dict = {}
	counts = Counter(seq_len)
	interval = 50
	for length in sorted(counts.keys()):
		len_distribution.writerow( [length, counts[length]] ) 
		group_num = int(math.floor(length/interval))   
		if group_num in group_dict.keys():        
			group_dict[group_num] += counts[length]
		else:
			group_dict[group_num] = counts[length]
	for indexs, x_label in enumerate(sorted(group_dict.keys())):
		label = "%s-%s"%((x_label)*interval,(x_label+1)*interval)
		len_distribution.writerow([label, group_dict[x_label]])
#
#--End -- Count length
#


#
# -- BEGIN -- Quality control
#
def quality_control(fqfile,phreds,lens):
	good_out = []
	bad_out = []
	for record in SeqIO.parse(fqfile,'fastq'):
		quality = record.letter_annotations["phred_quality"]
		new_record = record[retrive_lindex(quality,phreds) : retrive_rindex(quality,phreds)]	
		if len (new_record) >= lens:
			good_out.append(new_record)
		else:
			bad_out.append(record)
	return good_out,bad_out

def retrive_lindex(quality,phreds):
	for i, val in enumerate(quality):
  		if quality[i]>phreds:
			break
	return i
def retrive_rindex(quality,phreds):
	reverse_quality = quality[::-1]
	for j, val in enumerate(reverse_quality):
		if reverse_quality[j]>phreds:
			j = len(reverse_quality)-j
			break
	return j
#
# -- END -- Quality control
#

#
# -- BEGIN -- Getopt --
#	 
def process_parameters(parameters):
	#print parameters
	phred, lens = 20, 250
	try:
		opts,args = getopt.getopt(parameters,"hp:l:",["h","p","l"])
		for option,value in opts:
			#print option,value
			if option == "-p":
				phred = value
			if option == "-l":
				lens = value
			
			if option == "-h":	
				print "Usage: programname.py [option] [value] ..."
				print "-p: cut off phred in qulity control, default 20"
				print "-l: cut off length in qulity control,default 250 "
	except:
		print "ERROR!! Please use [python program_name.py -h] for help."	
	
	return (phred,lens)
#
# -- END -- Getopt --
#
	 
#
#-- BEGIN -- split FASTA file
#

def split_fasta(infile,lines):
	myrecord = SeqIO.parse(open(infile), "fasta")
	entry = True #Make sure we loop once
	while entry :
		batch = []
		while len(batch) < lines :
			try :
				entry = myrecord.next()
			except StopIteration :
				entry = None
			if entry is None :
				#End of file
				break
			batch.append(entry)
		if batch :
			yield batch
#
#-- END --
#

#
#--BEGIN--
#
#def all_igblast(infile)
