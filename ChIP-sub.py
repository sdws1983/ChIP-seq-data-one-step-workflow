#!/usr/bin/env python2

import sys, getopt
import re
import os
import time

def get_option():
	opts, args = getopt.getopt(sys.argv[1:], "hi:o:f:s:d:l:")
	input_file = "" #genes info file
	output_file = "" #output name
	ref = "" #gff3
	sam = "" #sam file
	divide = False
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

	return input_file,output_file,ref,sam,divide

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
	#   a : reads position		       #
	#   b : TSS/TTS position		     #
	#   c : output dir			   #
	
	if a < b + 2000 and a >= b - 2000:
		f = open(c, 'a')
		i = str((a - b) // 20)
		f.write(i + "\n")

def calculate_bin_negative(a, b, c):

	#   calculate reads in which bin (negative)  #
	#   a : reads position		       #
	#   b : TSS/TTS position		     #
	#   c : output dir			   #

	if b < a + 2000 and b >= a - 2000:
		f = open(c, 'a')
		i = str((b - a) // 20)
		f.write(i + "\n")

if __name__ == "__main__":

	i,divide,po_dir,read_dir,out_dir = get_option()
	diff_list_map(i,divide,po_dir,read_dir,out_dir)


