#!/usr/bin/env python2

import sys, getopt
import re
import os
import time
import pandas as pd
import numpy as np


def get_option():
	opts, args = getopt.getopt(sys.argv[1:], "hi:o:f:s:d:l:")
	input_file = "" #genes info file
	output_file = "" #output name
	ref = "" #gff3
	sam = "" #sam file
	divide = 1
	l = 35
	h = ""
	for op, value in opts:
		if op == "-i":
			input_file = value
		elif op == "-o":
			output_file = value
		elif op == "-f":
			ref = value
		elif op == "-s":
			sam = value
		elif op == "-d":
			divide = value
		elif op == "-l":
			l = int(value)
		elif op == "-h":
			h = 'useages:\nChIP-seq data processing workflow.\n-i : genes info\n-o : output name\n-f : gene annotation file (eg : gff3, gtf)\n-s : sam file\n-d : divide genes into <int> (default : False)\n-l : reads length (default : 35)'
	return input_file,output_file,ref,sam,divide,l,h

def sort_sam(sam,l):
	time_start = time.time()
	print ("----------------- first step : sort sam file -------------------")
	wk_dir = sam + "-sort"
	os.makedirs(wk_dir)
	with open(sam) as f:
		for i in f:
			i = i.split("\t")
			try:
				if i[1] == "0":
					fh = open((wk_dir + r"/" + i[2] + ".txt"), 'a')
					fh.write(str(int(i[3]) + 73) + "\n")
					fh.close()
				elif i[1] == "16":
					fh = open((wk_dir + r"/" + i[2] + ".txt"), 'a')
					fh.write(str(int(i[3]) - (73 - l)) + "\n")
					fh.close()
			except Exception as e:
				print ("error:" + i[0] + "\t" + str(e))
	print ("----------------- sort sam file okay .. -------------------")
	print ("time: " + str (time.time()-time_start))
					
def genes_info_processing(input_file,ref,divide):
	time_start = time.time()
	print ("----------------- second step : gene info processing -------------------")
	wk_dir = input_file + "-sort"
	os.makedirs(wk_dir)
	if divide:
		cmd = "sort -nrk 2 " + input_file + " |awk '{if($2>=1){print $0}}' > " + wk_dir + "/gene-sort.txt"
		print (cmd)
		os.system(cmd)
		f = open(wk_dir + "/gene-sort.txt")
		total = len(f.readlines())
		f.close()
		pair = total // int(divide)
		
		print ("divide into " + str(pair) + " genes each pair..")

		#processing gff3/gtf file

		cmd = ("awk '{if($3==" + '"gene"' + "){print $0}}' /copy1/a-reference/1-maize/0-all-genome/genes.gtf|awk -F '[\\t" + '"' + "]' '{print $1,$4,$5,$7,$10}' > " + wk_dir + "/gtf-sort.txt")
		print (cmd)
		os.system(cmd)

		#divide genes list..
		
		print ("divide genes list..")
		count = 1
		di = 1
		f = open(wk_dir + "/gene-pair-" + str(di), 'a')
		with open(wk_dir + "/gene-sort.txt") as fh:
			for i in fh:
				if count == pair:
					if di == int(divide):
						out = os.popen("grep '" + i.split("\t| ")[0] + "' " + wk_dir + "/gtf-sort.txt").read()
						f.write(out)
					else:
						out = os.popen("grep '" + i.split("\t| ")[0] + "' " + wk_dir + "/gtf-sort.txt").read()
						f.write(out)
						f.close()
						count = 0
						di += 1
						f = open(wk_dir + "/gene-pair-" + str(di), 'a')
				else:
					out = os.popen("grep '" + i.split("\t| ")[0] + "' " + wk_dir + "/gtf-sort.txt").read()
					f.write(out)
				count += 1
		print ("divide genes list okay..")

		#divide into different chrs..

		print ("divide into different chrs..")
		for i in range(1, int(divide) + 1):
			os.makedirs(wk_dir + "/pair-" + str(i) + ".+")
			os.makedirs(wk_dir + "/pair-" + str(i) + ".-")
			with open(wk_dir + "/gene-pair-" + str(i)) as f:
				for e in f:
					each = e[:-1].split(" ")
					if each[3] == "+":
						ff = open(wk_dir + "/pair-" + str(i) + ".+/" + each[0] + ".txt", 'a')
						ff.write(each[1] + "\t" + each[2] + "\n")
						ff.close()
					else:
						ff = open(wk_dir + "/pair-" + str(i) + ".-/" + each[0] + ".txt", 'a')
						ff.write(each[1] + "\t" + each[2] + "\n")
						ff.close()
		
		print ("divide into different chrs okay..")
		print ("----------------- gene info processing okay .. -------------------")
		print ("time: " + str (time.time()-time_start))

def mapping_reads(input_file,sam,output_file,divide):
	time_start = time.time()
	print ("----------------- third step : mapping reads -------------------")
	po_dir = input_file + "-sort/"
	read_dir = sam + "-sort/"
	out_dir = output_file + "/"
	os.makedirs(output_file)
	
	for n in range(int(divide)):
		os.makedirs(out_dir + "pair-" + str(n + 1) + "/TSS")
		os.makedirs(out_dir + "pair-" + str(n + 1) + "/TTS")	

	reads_file_list = os.popen("ls " + read_dir).read()[:-1].split("\n")
	for i in reads_file_list[10:]:
		diff_list_map(i,divide,po_dir,read_dir,out_dir)

	cmd = "./new.sh " + str(divide) + " " + po_dir + " " + read_dir + " " + out_dir
	print (cmd)
	os.system(cmd)
	print ("okay..")
	print ("----------------- mapping reads okay .. -------------------")
	print ("time: " + str (time.time()-time_start))

def calculate_frequency(input_file,sam,output_file,divide):
	time_start = time.time()
	print ("----------------- fourth step : calculating frequency in each bin -------------------")
	po_dir = input_file + "-sort/"
	read_dir = sam + "-sort/"
	out_dir = output_file + "/"

	# combine mapping results

	print ("combine mapping results..")
	all_reads_number = os.popen("wc -l " + read_dir + "*.txt|tail -1").read()[:-1].split(" ")[-2]
	for n in range(int(divide)):
		cmd1 = "cat " + out_dir + "pair-" + str(n + 1) + "/TSS/*.txt > " + out_dir + "pair-" + str(n + 1) + "/TSS-original.txt"
		cmd2 = "cat " + out_dir + "pair-" + str(n + 1) + "/TTS/*.txt > " + out_dir + "pair-" + str(n + 1) + "/TTS-original.txt"
		os.system(cmd1)
		os.system(cmd2)
		
		da = pd.read_table(out_dir + "pair-" + str(n + 1) + "/TSS-original.txt",header = None)
		da = da[0]
		counts = da.value_counts()
		gene_number = os.popen("wc -l " + po_dir + "gene-pair-" + str(n + 1)).read()[:-1].split(" ")[-2]
		f = open(out_dir + str(n + 1) + "-pair-TSS.txt", 'w')		

		for each in range(len(counts.index)):
			f.write((str(counts.index[each]) + "\t" + str(int(counts.values[each]) / (int(gene_number) * int(all_reads_number)))) + "\t" + "expression level " + str(n + 1) + "\n")
		f.close()

		da = pd.read_table(out_dir + "pair-" + str(n + 1) + "/TTS-original.txt",header = None)
		da = da[0]
		counts = da.value_counts()
		gene_number = os.popen("wc -l " + po_dir + "gene-pair-" + str(n + 1)).read()[:-1].split(" ")[-2]
		f = open(out_dir + str(n + 1) + "-pair-TTS.txt", 'w')           
    
		for each in range(len(counts.index)):
			f.write((str(counts.index[each]) + "\t" + str(int(counts.values[each]) / (int(gene_number) * int(all_reads_number)))) + "\t" + "expression level " + str(n + 1) + "\n")
		f.close()		

		print ("\tpair " + str(n + 1) + " is okay..")

	# combine calculating results

	cmd3 = "cat " + out_dir + "*TSS.txt > " + out_dir + "TSS.result.txt"
	print (cmd3)
#	f3 = open(out_dir + "TSS.result.txt", 'w')
	os.system(cmd3)
		
	cmd4 = "cat " + out_dir + "*TTS.txt > " + out_dir + "TTS.result.txt"
	print (cmd4)
#	f4 = open(out_dir + "TTS.result.txt", 'w')
	os.system(cmd4)
#	time.sleep(2)

	print ("----------------- calculating frequency in each bin okay .. -------------------")
	


def diff_list_map(i,divide,po_dir,read_dir,out_dir):
		print ("mapping chr " + str(i) + " ...")
		for n in range(int(divide)):

			print ("\tmapping to genes pair " + str(n + 1) + " ...")
			positive = po_dir + "pair-" + str(n + 1) + ".+/"
			negative = po_dir + "pair-" + str(n + 1) + ".-/"
			if i in os.popen("ls " + positive).read()[:-1].split("\n"):
				with open(read_dir + i) as f:
					for each in f:
						read_po = int(each[:-1])
						with open(positive + i) as fh:
							for t in fh:
								t = t[:-1].split("\t")
								calculate_bin_positive(read_po, int(t[0]), out_dir + "pair-" + str(n + 1) + "/TSS/" + i)# TSS position
								calculate_bin_positive(read_po, int(t[1]), out_dir + "pair-" + str(n + 1) + "/TTS/" + i)# TTS position
				print ("\t......positive is okay..")
			if i in os.popen("ls " + negative).read()[:-1].split("\n"):
				with open(read_dir + i) as f:
					for each in f:
						read_po = int(each[:-1])
						with open(negative + i) as fh:
							for t in fh:
								t = t[:-1].split("\t")
								calculate_bin_negative(read_po, int(t[1]), out_dir + "pair-" + str(n + 1) + "/TSS/" + i)# TSS position
								calculate_bin_negative(read_po, int(t[0]), out_dir + "pair-" + str(n + 1) + "/TTS/" + i)# TTS position
				print ("\t......negative is okay..")

def calculate_bin_positive(a, b, c):

	#   calculate reads in which bin (positive)  #
	#   a : reads position                       #
	#   b : TSS/TTS position                     #
	#   c : output dir                           #
	
	if a < b + 2000 and a >= b - 2000:
		f = open(c, 'a')
		i = str((a - b) // 20)
		f.write(i + "\n")

def calculate_bin_negative(a, b, c):

	#   calculate reads in which bin (negative)  #
	#   a : reads position                       #
	#   b : TSS/TTS position                     #
	#   c : output dir                           #

	if b < a + 2000 and b >= a - 2000:
		f = open(c, 'a')
		i = str((b - a) // 20)
		f.write(i + "\n")


def main(input_file,output_file,ref,sam,divide,l):
	sort_sam(sam,l)
	genes_info_processing(input_file,ref,divide)
	mapping_reads(input_file,sam,output_file,divide)
	calculate_frequency(input_file,sam,output_file,divide)

if __name__ == "__main__":
	time_start = time.time()

	input_file,output_file,ref,sam,divide,l,h = get_option()
	if str(h) == "":
		main(input_file,output_file,ref,sam,divide,l)
		print ("time: " + str (time.time()-time_start))
	else:
		print (h)
