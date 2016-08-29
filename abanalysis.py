#!/usr/bin/python
# encoding: utf-8
"""
Usage:1 for cut the low quality bases(<=20) in both ends and make sure the final length >=280
	  
parameter: counts, basecount, finalcount, ncontent.
Created by Zhang Yan on 2016-06-06

"""
import sys, os, glob, subprocess
from Bio import SeqIO
from zhangyan_tools import *
#from mytools import *
#from common_info.py import *

#def myigblast(infiles,num,prj_folder):
#	#cmd = '%s <  %s'%('bsub', IgBLAST_job)
#	#prj_folder = os.getcwd()
#	thisfname = os.path.splitext(os.path.split(infiles[int(num)])[1])[0]
#	IgBLAST_run = subprocess.call("igblastn -num_threads 4 -germline_db_V /zzh_gpfs02/zhangyan/20160314-Donor45-zhangyan/IgBLAST_database/20150429-human-gl-v -germline_db_J /zzh_gpfs02/zhangyan/20160314-Donor45-zhangyan/IgBLAST_database/20150429-human-gl-j -germline_db_D /zzh_gpfs02/zhangyan/20160314-Donor45-zhangyan/IgBLAST_database/20150429-human-gl-d -organism human -domain_system imgt -query %s -auxiliary_data /zzh_gpfs02/zhangyan/20150429-IgBLAST_optional_files/human_gl.aux -outfmt '7 qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore qlen slen qseq sseq score frames qframe sframe positive ppos btop staxids stitle sstrand qcovs qcovhsp' -num_alignments_V 10 -num_alignments_D 10 -num_alignments_J 10 -out %s/1.6-IgBLAST_result/IgBLAST_result_%s &"%(infiles[int(num)],prj_folder,thisfname,),shell=True)


def main():
	prj_folder = os.getcwd()
	print prj_folder, type(prj_folder)
	fname = str(prj_folder.split("/")[-1])
	#print fname
	
	'''
	### -- BEGIN -- combine fna file and qual file 
	print "Combining fna file and qual file"
	fna = glob.glob("%s/rawdata/*.fna"%(prj_folder))
	qual = glob.glob("%s/rawdata/*.qual"%(prj_folder))
	if len(fna) == 1:
		fq = open("%s/rawdata/%s.fastq"%(prj_folder,fname), "w")
		records = PairedFastaQualIterator(fna[0],qual[0])
		counts = SeqIO.write(records, fq, "fastq")
		fq.close()
	else:
		print "No or more than FNA file"
		pass

	
	### -- BEGIN -- Count length
	print "Count length"
	infiles = glob.glob("%s/rawdata/*.fastq"%(prj_folder))
	for infile in infiles:
		print infile
		count_length(infile)
	

	
	### -- BEGIN -- QC
	filecount = 0	
	print "QC-ing,length:%s, phred:%s"%(lens,phreds)
	infiles = glob.glob("%s/rawdata/*.fastq"%(prj_folder))
	if len(infiles) == 2:
		mergeflag = "merge"
	if len(infiles) == 1:
		mergeflag = "not-merge"
	#infiles = gzip.GzipFile(infile)
	for infile in infiles:
		wholefname = os.path.splitext(os.path.split(infile)[1])[0]
		filenum = int(wholefname.split('-')[-1])
		good_out = open("%s/1.1-QC-fastq-file/%s_passed-%d.fastq" %(prj_folder,fname,filenum), 'w')
		bad_out = open("%s/1.1-QC-fastq-file/%s_unpassed-%d.fq" % (prj_folder,fname,filenum), 'w')
		goods, bads = quality_control(infile,phreds,lens)		## use zhangyan_tools
		for iterms in goods:
			SeqIO.write(iterms, good_out, "fastq")
		for iterms in bads:
			SeqIO.write(iterms, bad_out, "fastq")
		print "%s_%d_raw data has %d records"%(fname,filecount,len(goods)+len(bads))
		print "%s_%d_finally we get %d result"%(fname,filecount,len(goods))
		#print "%s_we trimmed %d records"%(fname,trimcount)
	#else:
	#	print "No or more than 1 fastq file in RAWDATA. Please check."
	#	sys.exit(0)
	
	
	##-- BEGIN -- intersection
	infiles = glob.glob("%s/1.1-QC-fastq-file/*.fastq"%(prj_folder))
	if len(infiles) == 2:
		leftfile = glob.glob("%s/1.1-QC-fastq-file/*-1.fastq"%(prj_folder))[0]
		rightfile = glob.glob("%s/1.1-QC-fastq-file/*-2.fastq"%(prj_folder))[0]
		lefts = open("%s/1.2-for-merge/%s-1.fastq"%(prj_folder,fname),'w')
		rights = open("%s/1.2-for-merge/%s-2.fastq"%(prj_folder,fname),'w')
		leftid, rightid = set(),set()
		for record in SeqIO.parse(leftfile,"fastq"):
			leftid.add(record.id)
		for record in SeqIO.parse(rightfile,"fastq"):
			rightid.add(record.id)
		intersectionid = leftid&rightid
		for record in SeqIO.parse(leftfile,"fastq"):
			if record.id in intersectionid:
				SeqIO.write(record,lefts, "fastq")
			else:
				continue
		for record in SeqIO.parse(rightfile,"fastq"):
			if record.id in intersectionid:
				SeqIO.write(record, rights, "fastq")
	else:
		print "not intersection"


	
	###-- BEGIN -- merge
	infiles = glob.glob("%s/1.1-QC-fastq-file/*.fastq"%(prj_folder))
	if len(infiles) == 2:
		leftfile = "%s/1.2-for-merge/%s-1.fastq"%(prj_folder,fname)
		rightfile = "%s/1.2-for-merge/%s-2.fastq"%(prj_folder,fname)
		print "Merge-ing %s and %s"%(leftfile,rightfile)
		merge = subprocess.call("pear -f %s -r %s -o %s/1.2-merged-file/merged -j 16"%(leftfile,rightfile,prj_folder),shell=True)
		mergeflag = "merge"
	elif len(infiles) == 1:
		print "SE sequencing, just 1 FASTQ file!"
		mergeflag = "not-merge"
		pass
	else:
		print "Not merge, because more than 2 files or no file."
		pass
	

	###-- BEGIN -- fastq to fasta
	print "Convert fastq into fasta-ing"
	if mergeflag == "merge":
		infile = glob.glob("%s/1.2-merged-file/merged.assembled.fastq"%prj_folder)[0]
		SeqIO.convert(infile, "fastq", "%s/1.3-convert/%s_merged.fasta"%(prj_folder,fname), "fasta")	
		infiles = glob.glob("%s/1.1-QC-fastq-file/*.fastq"%prj_folder)
		for infile in infiles:
			wholefname = os.path.splitext(os.path.split(infile)[1])[0]
			filenum = int(wholefname.split('-')[-1])
			SeqIO.convert(infile, "fastq", "%s/1.3-convert/%s-%d.fasta"%(prj_folder,fname,filenum), "fasta")
	if mergeflag == "not-merge":
		infile = glob.glob("%s/1.1-QC-fastq-file/*.fastq"%prj_folder)[0]
		SeqIO.convert(infile, "fastq", "%s/1.3-convert/%s.fasta"%(prj_folder,fname), "fasta")     
		
	
	
	### --BEGIN-- split file
	print "Split file for Igblast"
	infiles = glob.glob("%s/1.3-convert/*.fasta"%prj_folder)
	filenum = len(infiles)
	for flag in range(0,filenum):
	   	for indexs,rec in enumerate(split_fasta(infiles[flag],10000)):
			thisfname = os.path.splitext(os.path.split(infiles[flag])[1])[0]
			filename = "%s_%s.fa"%(thisfname,indexs+1)
			myfolder = open("%s/1.4-split/%s"%(prj_folder,filename),'w+')
			counts = SeqIO.write(rec,myfolder,"fasta")
			print "%d reads in %s"%(counts,filename)
	
	'''

    ###-- BEGIN -- Igblast
	print "Igblast-ing"
	#infiles = glob.glob("%s/1.4-split/*.fa"%prj_folder)
	#igblasts = "igblastn -germline_db_V /zzh_gpfs02/zhangyan/20160314-Donor45-zhangyan/IgBLAST_database/20150429-human-gl-v -germline_db_J /zzh_gpfs02/zhangyan/20160314-Donor45-zhangyan/IgBLAST_database/20150429-human-gl-j -germline_db_D /zzh_gpfs02/zhangyan/20160314-Donor45-zhangyan/IgBLAST_database/20150429-human-gl-d -organism human -domain_system imgt -query $i -auxiliary_data /zzh_gpfs02/zhangyan/20150429-IgBLAST_optional_files/human_gl.aux -outfmt '7 qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore qlen slen qseq sseq score frames qframe sframe positive ppos btop staxids stitle sstrand qcovs qcovhsp' -num_alignments_V 10 -num_alignments_D 10 -num_alignments_J 10 -out ../1.7-IgBLAST_result/IgBLAST_result_$i"
	#cmds = "for i in `ls *.fa`;do bsub -n 1 -q cpu -e errput_$i -o output_$i %s;done"%igblasts
	cd_cmd = "sh igblast.sh"
	#pre_IgBLAST = subprocess.call(cd_cmd, shell=True)
	run_IgBLAST = subprocess.call(cd_cmd, shell=True)
	#file_numbers =len(infiles)
	#for num in range(0,file_numbers):
	#	myigblast(infiles,num,prj_folder)
			

if __name__ == '__main__':
	phreds, lens = process_parameters(sys.argv[1:])		## use zhangyan_tools
	if "-h" in sys.argv[1:]:
		sys.exit(0)
	main()
	#print process_parameters(sys.argv[1:])
	#print "finished"

