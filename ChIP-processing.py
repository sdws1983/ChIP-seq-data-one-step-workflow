#!/usr/bin/env python2

import sys, getopt
import re
import os
import time
import pandas as pd
import numpy as np


def get_option():
	opts, args = getopt.getopt(sys.argv[1:], "hi:o:f:s:d:l:")
	input_file = "" #genes info file ( ps : split by '\t' )
	output_file = "" #output name
	ref = "" #gtf
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
			h = 'useages:\nChIP-seq data processing workflow.\n-i : genes info (split by tab)\n-o : output name\n-f : gene annotation file (only support gtf file)\n-s : sam file\n-d : divide genes into <int> (default : 1)\n-l : reads length (default : 35)'
	return input_file,output_file,ref,sam,divide,l,h

def sort_sam(sam,l):
	time_start = time.time()
	print ("----------------- first step : sort sam file -------------------")
	wk_dir = sam + "-tmp"
	wk_dir_sort = sam + "-sort"
	if os.path.exists(wk_dir):
		print ('tmp file existed..')
		pass
	else:
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
	if os.path.exists(wk_dir_sort):
		pass
	else:
		os.makedirs(wk_dir_sort)
		sam_list = os.popen("ls " + wk_dir + "/").read()[:-1].split("\n")
		print (sam_list)
		for i in sam_list:
			os.system("sort -nk 1 " + wk_dir + "/" + i + " > " + wk_dir_sort + "/" + i)

	print ("----------------- sort sam file okay .. -------------------")
	print ("time: " + str (time.time()-time_start))
					
def genes_info_processing(input_file,ref,divide):
	time_start = time.time()
	print ("----------------- second step : gene info processing -------------------")
	wk_dir = input_file + "-sort"
	if os.path.exists(wk_dir):
		print ("----------------- gene info file exist .. -------------------")
		print ("time: " + str (time.time()-time_start))
		return 1

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
						out = os.popen("grep '" + i.split("\t")[0] + "' " + wk_dir + "/gtf-sort.txt").read()
						f.write(out)
					else:
						out = os.popen("grep '" + i.split("\t")[0] + "' " + wk_dir + "/gtf-sort.txt").read()
						f.write(out)
						f.close()
						count = 0
						di += 1
						f = open(wk_dir + "/gene-pair-" + str(di), 'a')
				else:
					out = os.popen("grep '" + i.split("\t")[0] + "' " + wk_dir + "/gtf-sort.txt").read()
					f.write(out)
				count += 1
		f.close()

		print ("divide genes list okay..")

		#divide into different chrs..

		print ("divide into different chrs..")
		for i in range(1, int(divide) + 1):
			os.makedirs(wk_dir + "/pair-" + str(i) + ".+")
			os.makedirs(wk_dir + "/pair-" + str(i) + ".-")

			cmd = "sort -nk 2 " + wk_dir + "/gene-pair-" + str(i) + "|awk '{print $2" + '"\\t"' + '$3 > "' + wk_dir + "/pair-" + str(i) + '."$4"/' + '"$1".txt"' + "}'"
			print (cmd)
			os.system(cmd)

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
		os.makedirs(out_dir + "pair-" + str(n + 1) + "/body")
		os.makedirs(out_dir + "pair-" + str(n + 1) + "/TTS")	

	reads_file_list = os.popen("ls " + read_dir).read()[:-1].split("\n")
	for i in reads_file_list[10:]:
		diff_list_map(i,divide,po_dir,read_dir,out_dir)

	cmd = "/home/yumin/0-script/2-NGS/ChIP-workflow/scripts/new.sh " + str(divide) + " " + po_dir + " " + read_dir + " " + out_dir
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
			f.write((str(counts.index[each]) + "\t" + str(int(counts.values[each]) / (int(gene_number) * int(all_reads_number)))) + "\t" + "expression level " + str(n + 1) + "\tTSS\n")
		f.close()

		da = pd.read_table(out_dir + "pair-" + str(n + 1) + "/TTS-original.txt",header = None)
		da = da[0]
		counts = da.value_counts()
		gene_number = os.popen("wc -l " + po_dir + "gene-pair-" + str(n + 1)).read()[:-1].split(" ")[-2]
		f = open(out_dir + str(n + 1) + "-pair-TTS.txt", 'w')           
    
		for each in range(len(counts.index)):
			f.write((str(counts.index[each]) + "\t" + str(int(counts.values[each]) / (int(gene_number) * int(all_reads_number)))) + "\t" + "expression level " + str(n + 1) + "\tTTS\n")
		f.close()		

		print ("\tpair " + str(n + 1) + " is okay..")

	# combine calculating results

	cmd3 = "cat " + out_dir + "*TSS.txt > " + out_dir + "TSS.result.txt"
	print (cmd3)
	os.system(cmd3)
	os.system("sed -i '1 iposition\tnormalized-read-counts\tgroup\tdiff' " + out_dir + "TSS.result.txt")
		
	cmd4 = "cat " + out_dir + "*TTS.txt > " + out_dir + "TTS.result.txt"
	print (cmd4)
	os.system(cmd4)
	os.system("sed -i '1 iposition\tnormalized-read-counts\tgroup\tdiff' " + out_dir + "TTS.result.txt")

	cmd5 = "cat " + out_dir + "*TTS.txt " + out_dir + "*TSS.txt > " + out_dir + "all.result.txt"
	print (cmd5)
	os.system(cmd5)
	os.system("sed -i '1 iposition\tnormalized-read-counts\tgroup\tdiff' " + out_dir + "all.result.txt")

	print ("----------------- calculating frequency in each bin okay .. -------------------")
	

def diff_list_map(i,divide,po_dir,read_dir,out_dir):
		print ("mapping chr " + str(i) + " ...")
		for n in range(int(divide)):

			print ("\tmapping to genes pair " + str(n + 1) + " ...")
			positive = po_dir + "pair-" + str(n + 1) + ".+/"
			negative = po_dir + "pair-" + str(n + 1) + ".-/"
			if i in os.popen("ls " + positive).read()[:-1].split("\n"):

				# TSS

				with open(read_dir + i) as f:
					gi_list = open(positive + i, 'r').read()[:-1].split('\n')
#					print (gi_list[:50])
					for each in f:
						read_po = int(each[:-1])
						if read_po < ( int(gi_list[0].split('\t')[0]) - 2000 ):
							continue
						elif read_po < ( int(gi_list[0].split('\t')[0]) + 2000 ):
							calculate_bin_positive(read_po, int(gi_list[0].split('\t')[0]), "/TSS/", out_dir + "pair-" + str(n + 1), i) #!!!!
							for u in gi_list[1:]:
								if read_po < ( int(u.split('\t')[0]) - 2000 ):
									break
								elif read_po < ( int(u.split('\t')[0]) + 2000 ):
									calculate_bin_positive(read_po, int(u.split('\t')[0]), "/TSS/", out_dir + "pair-" + str(n + 1), i) #!!!!

						else:
							try:
								while read_po >= ( int(gi_list[0].split('\t')[0]) + 2000 ):
									gi_list.pop(0)
							except:
								break

							for u in gi_list:
								if read_po < ( int(u.split('\t')[0]) - 2000 ):
									break
								elif read_po < ( int(u.split('\t')[0]) + 2000 ):
									calculate_bin_positive(read_po, int(u.split('\t')[0]), "/TSS/", out_dir + "pair-" + str(n + 1), i) #!!!!

				# TTS

				with open(read_dir + i) as f:
					gi_list = open(positive + i, 'r').read()[:-1].split('\n')
					#print (gi_list[:50])
					for each in f:
						read_po = int(each[:-1])
						if read_po < ( int(gi_list[0].split('\t')[1]) - 2000 ):
							continue
						elif read_po < ( int(gi_list[0].split('\t')[1]) + 2000 ):
							calculate_bin_positive(read_po, int(gi_list[0].split('\t')[1]), "/TTS/", out_dir + "pair-" + str(n + 1), i) #!!!!
							for u in gi_list[1:]:
								if read_po < ( int(u.split('\t')[1]) - 2000 ):
									break
								elif read_po < ( int(u.split('\t')[1]) + 2000 ):
									calculate_bin_positive(read_po, int(u.split('\t')[1]), "/TTS/", out_dir + "pair-" + str(n + 1), i) #!!!!

						else:
							try:
								while read_po >= ( int(gi_list[0].split('\t')[1]) + 2000 ):
									gi_list.pop(0)
							except:
								break

							for u in gi_list:
								if read_po < ( int(u.split('\t')[1]) - 2000 ):
									break
								elif read_po < ( int(u.split('\t')[1]) + 2000 ):
									calculate_bin_positive(read_po, int(u.split('\t')[1]), "/TTS/", out_dir + "pair-" + str(n + 1), i) #!!!!

				print ("\t......positive is okay..")

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------------------------------------------ #

			if i in os.popen("ls " + negative).read()[:-1].split("\n"):

				# TTS

				with open(read_dir + i) as f:#open reads file
					gi_list = open(negative + i, 'r').read()[:-1].split('\n')
				#	print (gi_list[:50])
					for each in f:
						read_po = int(each[:-1])
						if read_po <= ( int(gi_list[0].split('\t')[0]) - 2000 ):
							continue
						elif read_po <= ( int(gi_list[0].split('\t')[0]) + 2000 ):
							calculate_bin_negative(read_po, int(gi_list[0].split('\t')[0]), "/TTS/", out_dir + "pair-" + str(n + 1), i)
							for u in gi_list[1:]:
								if read_po <= ( int(u.split('\t')[0]) - 2000 ):
									break
								elif read_po <= ( int(u.split('\t')[0]) + 2000 ):
									calculate_bin_negative(read_po, int(u.split('\t')[0]), "/TTS/", out_dir + "pair-" + str(n + 1), i)

						else:
							try:
								while read_po > ( int(gi_list[0].split('\t')[0]) + 2000 ):
									gi_list.pop(0)
							except:
								break

							for u in gi_list:
								if read_po <= ( int(u.split('\t')[0]) - 2000 ):
									break
								elif read_po <= ( int(u.split('\t')[0]) + 2000 ):
									calculate_bin_negative(read_po, int(u.split('\t')[0]), "/TTS/", out_dir + "pair-" + str(n + 1), i)

				# TSS

				with open(read_dir + i) as f:#open reads file
					gi_list = open(negative + i, 'r').read()[:-1].split('\n')
				#   print (gi_list[:50])
					for each in f:
						read_po = int(each[:-1])
						if read_po <= ( int(gi_list[0].split('\t')[1]) - 2000 ):
							continue
						elif read_po <= ( int(gi_list[0].split('\t')[1]) + 2000 ):
							calculate_bin_negative(read_po, int(gi_list[0].split('\t')[1]), "/TSS/", out_dir + "pair-" + str(n + 1), i)
							for u in gi_list[1:]:
								if read_po <= ( int(u.split('\t')[1]) - 2000 ):
									break
								elif read_po <= ( int(u.split('\t')[1]) + 2000 ):
									calculate_bin_negative(read_po, int(u.split('\t')[1]), "/TSS/", out_dir + "pair-" + str(n + 1), i)

						else:
							try:
								while read_po > ( int(gi_list[0].split('\t')[1]) + 2000 ):
									gi_list.pop(0)
							except:
								break

							for u in gi_list:
								if read_po <= ( int(u.split('\t')[1]) - 2000 ):
									break
								elif read_po <= ( int(u.split('\t')[1]) + 2000 ):
									calculate_bin_negative(read_po, int(u.split('\t')[1]), "/TSS/", out_dir + "pair-" + str(n + 1), i)

				print ("\t......negative is okay..")


def calculate_bin_positive(a, b, c, d, e):

	#   calculate reads in which bin (positive)  #
	#   a : reads position                       #
	#   b : TSS/TTS position                     #
	#   c : TSS/TTS status                       #
	#   d : output dir                           #
	#   e : chr number                           #	

	f = open(d + c + str(e), 'a')
	f.write(str (( a - b ) // 20 ) + "\n")

def calculate_bin_negative(a, b, c, d, e):

	#   calculate reads in which bin (negative)  #
	#   a : reads position                       #
	#   b : TSS/TTS position                     #
	#   c : TSS/TTS status                       #
	#   d : output dir                           #
	#   e : chr number                           #
	
	f = open(d + c + str(e), 'a')
	f.write(str (( b - a ) // 20 ) + "\n")


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
