#!/usr/bin/env python2

import sys, getopt
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

def diff_list_map(i, divide, po_dir, read_dir, out_dir):
		print("mapping chr " + str(i) + " ...")
		for n in range(int(divide)):

			print("\tmapping to genes pair " + str(n + 1) + " ...")
			positive = po_dir + "pair-" + str(n + 1) + ".+/"
			negative = po_dir + "pair-" + str(n + 1) + ".-/"
			if i in os.popen("ls " + positive).read()[:-1].split("\n"):

				# TSS

				with open(read_dir + i) as f:
					gi_list = open(positive + i, 'r').read()[:-1].split('\n')
					#					print (gi_list[:50])
					for each in f:
						read_po = int(each[:-1])
						if read_po < (int(gi_list[0].split('\t')[0]) - 2000):
							continue
						elif read_po < (int(gi_list[0].split('\t')[0]) + 2000):
							calculate_bin_positive(read_po, int(gi_list[0].split('\t')[0]), "/TSS/",
							                       out_dir + "pair-" + str(n + 1), i)  # !!!!
							for u in gi_list[1:]:
								if read_po < (int(u.split('\t')[0]) - 2000):
									break
								elif read_po < (int(u.split('\t')[0]) + 2000):
									calculate_bin_positive(read_po, int(u.split('\t')[0]), "/TSS/",
									                       out_dir + "pair-" + str(n + 1), i)  # !!!!

						else:
							try:
								while read_po >= (int(gi_list[0].split('\t')[0]) + 2000):
									gi_list.pop(0)
							except:
								break

							for u in gi_list:
								if read_po < (int(u.split('\t')[0]) - 2000):
									break
								elif read_po < (int(u.split('\t')[0]) + 2000):
									calculate_bin_positive(read_po, int(u.split('\t')[0]), "/TSS/",
									                       out_dir + "pair-" + str(n + 1), i)  # !!!!

				# TTS

				with open(read_dir + i) as f:
					gi_list = open(positive + i, 'r').read()[:-1].split('\n')
					# print (gi_list[:50])
					for each in f:
						read_po = int(each[:-1])
						if read_po < (int(gi_list[0].split('\t')[1]) - 2000):
							continue
						elif read_po < (int(gi_list[0].split('\t')[1]) + 2000):
							calculate_bin_positive(read_po, int(gi_list[0].split('\t')[1]), "/TTS/",
							                       out_dir + "pair-" + str(n + 1), i)  # !!!!
							for u in gi_list[1:]:
								if read_po < (int(u.split('\t')[1]) - 2000):
									break
								elif read_po < (int(u.split('\t')[1]) + 2000):
									calculate_bin_positive(read_po, int(u.split('\t')[1]), "/TTS/",
									                       out_dir + "pair-" + str(n + 1), i)  # !!!!

						else:
							try:
								while read_po >= (int(gi_list[0].split('\t')[1]) + 2000):
									gi_list.pop(0)
							except:
								break

							for u in gi_list:
								if read_po < (int(u.split('\t')[1]) - 2000):
									break
								elif read_po < (int(u.split('\t')[1]) + 2000):
									calculate_bin_positive(read_po, int(u.split('\t')[1]), "/TTS/",
									                       out_dir + "pair-" + str(n + 1), i)  # !!!!

				print("\t......positive is okay..")

			# ------------------------------------------------------------------------------------------------------------------------------------------ #
			# ------------------------------------------------------------------------------------------------------------------------------------------ #

			if i in os.popen("ls " + negative).read()[:-1].split("\n"):

				# TTS

				with open(read_dir + i) as f:  # open reads file
					gi_list = open(negative + i, 'r').read()[:-1].split('\n')
					#	print (gi_list[:50])
					for each in f:
						read_po = int(each[:-1])
						if read_po <= (int(gi_list[0].split('\t')[0]) - 2000):
							continue
						elif read_po <= (int(gi_list[0].split('\t')[0]) + 2000):
							calculate_bin_negative(read_po, int(gi_list[0].split('\t')[0]), "/TTS/",
							                       out_dir + "pair-" + str(n + 1), i)
							for u in gi_list[1:]:
								if read_po <= (int(u.split('\t')[0]) - 2000):
									break
								elif read_po <= (int(u.split('\t')[0]) + 2000):
									calculate_bin_negative(read_po, int(u.split('\t')[0]), "/TTS/",
									                       out_dir + "pair-" + str(n + 1), i)

						else:
							try:
								while read_po > (int(gi_list[0].split('\t')[0]) + 2000):
									gi_list.pop(0)
							except:
								break

							for u in gi_list:
								if read_po <= (int(u.split('\t')[0]) - 2000):
									break
								elif read_po <= (int(u.split('\t')[0]) + 2000):
									calculate_bin_negative(read_po, int(u.split('\t')[0]), "/TTS/",
									                       out_dir + "pair-" + str(n + 1), i)

				# TSS

				with open(read_dir + i) as f:  # open reads file
					gi_list = open(negative + i, 'r').read()[:-1].split('\n')
					#   print (gi_list[:50])
					for each in f:
						read_po = int(each[:-1])
						if read_po <= (int(gi_list[0].split('\t')[1]) - 2000):
							continue
						elif read_po <= (int(gi_list[0].split('\t')[1]) + 2000):
							calculate_bin_negative(read_po, int(gi_list[0].split('\t')[1]), "/TSS/",
							                       out_dir + "pair-" + str(n + 1), i)
							for u in gi_list[1:]:
								if read_po <= (int(u.split('\t')[1]) - 2000):
									break
								elif read_po <= (int(u.split('\t')[1]) + 2000):
									calculate_bin_negative(read_po, int(u.split('\t')[1]), "/TSS/",
									                       out_dir + "pair-" + str(n + 1), i)

						else:
							try:
								while read_po > (int(gi_list[0].split('\t')[1]) + 2000):
									gi_list.pop(0)
							except:
								break

							for u in gi_list:
								if read_po <= (int(u.split('\t')[1]) - 2000):
									break
								elif read_po <= (int(u.split('\t')[1]) + 2000):
									calculate_bin_negative(read_po, int(u.split('\t')[1]), "/TSS/",
									                       out_dir + "pair-" + str(n + 1), i)

				print("\t......negative is okay..")


def calculate_bin_positive(a, b, c, d, e):
	#   calculate reads in which bin (positive)  #
	#   a : reads position                       #
	#   b : TSS/TTS position                     #
	#   c : TSS/TTS status                       #
	#   d : output dir                           #
	#   e : chr number                           #

	f = open(d + c + str(e), 'a')
	f.write(str((a - b) // 20) + "\n")


def calculate_bin_negative(a, b, c, d, e):
	#   calculate reads in which bin (negative)  #
	#   a : reads position                       #
	#   b : TSS/TTS position                     #
	#   c : TSS/TTS status                       #
	#   d : output dir                           #
	#   e : chr number                           #

	f = open(d + c + str(e), 'a')
	f.write(str((b - a) // 20) + "\n")


if __name__ == "__main__":

	i,divide,po_dir,read_dir,out_dir = get_option()
	diff_list_map(i,divide,po_dir,read_dir,out_dir)


