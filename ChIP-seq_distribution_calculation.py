#!/usr/bin/env python2

import re
import time
import pandas as pd
import numpy as np
import sys
import argparse
import os
import textwrap

def version():
        v = "----"
        return v

def warning(*objs):
        print("WARNING: ", *objs, end='\n', file=sys.stderr)
        sys.exit()

def get_parser():
        parser = argparse.ArgumentParser(
                formatter_class=argparse.RawDescriptionHelpFormatter,
                description=textwrap.dedent(version())
        )

        parser.add_argument('-i','--input', help='gene info', type=str)
        parser.add_argument('-f','--gtf', help='gene annotation file (gtf)', type=str)
        parser.add_argument('-b','--bam', help='bam file', type=str)
        parser.add_argument('-d','--divide', help='divide genes into', type=str)
        parser.add_argument('-o', '--output', help='output name', type=str)
        
        return parser


def sort_bam(bam):
	time_start = time.time()
	print ("----------------- first step : sort bam file -------------------")
		
	command = "bamToBed -i " + bam + " | awk '{a=($2+$3)/2; printf $1" + '"\\t";printf "%d", a;printf "\\t";printf "%d\\n",a' + "}'|sort -k1,1 -k2,2n > " + bam + ".bed" 
	print ("run : " + command)
	os.system(command)

	print ("----------------- sort bam file okay .. -------------------")
	print ("time: " + str (time.time()-time_start))
					
def genes_info_processing(input_file,ref,divide,bam):
	time_start = time.time()
	print ("----------------- second step : gene info processing -------------------")
	wk_dir = input_file + "-sort"
	bam = bam +".bed"
	if os.path.exists(wk_dir):
		print ("----------------- gene info file exist .. -------------------")
		print ("time: " + str (time.time()-time_start))
		return 1

	os.makedirs(wk_dir)
	if divide:
		cmd1 = "sort -nrk 2 " + input_file + " |awk '{if($2>0){print $0}}' > " + wk_dir + "/gene-sort.txt"
		cmd2 = "cat " + input_file + " |awk '{if($2==0){print $0}}' > " + wk_dir + "/non-expressed-gene.txt"
		print (cmd1)
		print (cmd2)
		os.system(cmd1)
		os.system(cmd2)

		f = open(wk_dir + "/gene-sort.txt")
		total = len(f.readlines())
		f.close()
		pair = total // int(divide)
		
		print ("divide into " + str(pair) + " genes each pair..")

		#processing gff3/gtf file

		cmd = ("awk '{if($3==" + '"gene"' + "){print $0}}' " + ref + "|awk -F '[\\t" + '"' + "]' '{print " + '$1"\\t"$4"\\t"$5"\\t"$10"\\t"$7' + "}' > " + wk_dir + "/gtf-sort.txt")
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

		print ("process non-expressed genes...")		

		f = open(wk_dir + "/gene-pair-0", 'a')
		with open(wk_dir + "/non-expressed-gene.txt") as fh:
			for i in fh:
				out = os.popen("grep '" + i.split("\t")[0] + "' " + wk_dir + "/gtf-sort.txt").read()
				f.write(out)
		f.close()

		print ("----------------- gene info processing okay .. -------------------")
		print ("time: " + str (time.time()-time_start))
	

def calculate_distribution(input_file,ref,divide,bam,output):
		time_start = time.time()
		print ("----------------- thrid step : calculate distribution -------------------")

		bam = bam + ".bed"
		wk_dir = input_file + "-sort"	
		os.makedirs(wk_dir + "/spliting")#sort -k1,1 -k2,2n
		os.makedirs(wk_dir + "/spliting-sort")
		os.makedirs(wk_dir + "/output")
		os.makedirs(output)
		
		size = 2190000000
		read_num = int(os.popen("wc -l " + bam).read()[:-1].split(" ")[0])

		for a in range(0, int(divide) + 1):
			f1 = open(wk_dir + "/spliting/pair-" + str(a) + ".TSS.bed", "a")
			f2 = open(wk_dir + "/spliting/pair-" + str(a) + ".TTS.bed", "a")
			f3 = open(wk_dir + "/spliting/pair-" + str(a) + ".body.bed", "a")

			gene_num = int(os.popen("wc -l " + wk_dir + "/gene-pair-" + str(a)).read()[:-1].split(" ")[0])
			print ("gene number : " + str(gene_num))

			with open(wk_dir + "/gene-pair-" + str(a)) as f:
				for i in f:
					i = i[:-1].split("\t")
					if int(i[1]) >= 1000:
						if i[4] == "+":
							n = 1
							for each in range(int(i[1])-1000,int(i[1]),10):
								f1.write(i[0] + "\t" + str(each) + "\t" + str(each + 10) + "\t" + i[3] + "_tss_" + str(n) + "\n")
								n += 1
							n = 1
							for each in range(int(i[2]),int(i[2])+1000,10):
								f2.write(i[0] + "\t" + str(each) + "\t" + str(each + 10) + "\t" + i[3] + "_tts_" + str(n) + "\n")
								n += 1
						if i[4] == "-":
							n = 1
							for each in range(int(i[2])+1000,int(i[2]),-10):
								f1.write(i[0] + "\t" + str(each - 10) + "\t" + str(each) + "\t" + i[3] + "_tss_" + str(n) + "\n")
								n += 1
							n = 1
							for each in range(int(i[1]),int(i[1])-1000,-10):
								f2.write(i[0] + "\t" + str(each - 10) + "\t" + str(each) + "\t" + i[3] + "_tts_" + str(n) + "\n")
								n += 1
					else:
						print (i)
						if i[4] == "+":
							n = 100
							for each in range(int(i[1]),0,-10):
								each1 = (each -10) if (each -10 > 0) else 0
								f1.write(i[0] + "\t" + str(each1) + "\t" + str(each) + "\t" + i[3] + "_tss_" + str(n) + "\n")
								n -= 1
							n = 1
							for each in range(int(i[2]),int(i[2])+1000,10):
								f2.write(i[0] + "\t" + str(each) + "\t" + str(each + 10) + "\t" + i[3] + "_tts_" + str(n) + "\n")
								n += 1
						else:
							n = 1
							for each in range(int(i[2])+1000,int(i[2]),-10):
								f1.write(i[0] + "\t" + str(each - 10) + "\t" + str(each) + "\t" + i[3] + "_tss_" + str(n) + "\n")
								n += 1
							n = 1
							for each in range(int(i[1]),int(i[1])//100,-10):
								each1 = (each -10) if (each -10 > 0) else 0
								f2.write(i[0] + "\t" + str(each1) + "\t" + str(each) + "\t" + i[3] + "_tts_" + str(n) + "\n")
								n += 1
					if i[4] == "+":
						n = 1
						for each in np.arange(int(i[1]), int(i[2]), (int(i[2]) - int(i[1]))/100):
							f3.write(i[0] + "\t" + str(int(round(each))) + "\t" + str(int(round(each + (int(i[2]) - int(i[1]))/100))) + "\t" + i[3] + "_body_" + str(n) + "\n")
							n += 1
							if n == 101:
								break
					else:
						n = 100
						for each in np.arange(int(i[1]), int(i[2]), (int(i[2]) - int(i[1]))/100):
							f3.write(i[0] + "\t" + str(int(round(each))) + "\t" + str(int(round(each + (int(i[2]) - int(i[1]))/100))) + "\t" + i[3] + "_body_" + str(n) + "\n")
							n -= 1
							if n == 0:
								break
			f1.close()
			f2.close()
			f3.close()
			os.system("sort -k1,1 -k2,2n " + wk_dir + "/spliting/pair-" + str(a) + ".TSS.bed > "  + wk_dir + "/spliting-sort/pair-" + str(a) + ".TSS.bed")
			os.system("sort -k1,1 -k2,2n " + wk_dir + "/spliting/pair-" + str(a) + ".TTS.bed > "  + wk_dir + "/spliting-sort/pair-" + str(a) + ".TTS.bed")
			os.system("sort -k1,1 -k2,2n " + wk_dir + "/spliting/pair-" + str(a) + ".body.bed > "  + wk_dir + "/spliting-sort/pair-" + str(a) + ".body.bed")

			cmd1 = "bedtools coverage -sorted -a " + wk_dir + "/spliting-sort/pair-" + str(a) + ".TSS.bed" + " -b " + bam + " > " + wk_dir + "/output/pair-" + str(a) + ".TSS.out"
			print (cmd1)
			os.system(cmd1)
			cmd2 = "bedtools coverage -sorted -a " + wk_dir + "/spliting-sort/pair-" + str(a) + ".TTS.bed" + " -b " + bam + " > " + wk_dir + "/output/pair-" + str(a) + ".TTS.out"
			print (cmd2)
			os.system(cmd2)
			cmd3 = "bedtools coverage -sorted -a " + wk_dir + "/spliting-sort/pair-" + str(a) + ".body.bed" + " -b " + bam + " > " + wk_dir + "/output/pair-" + str(a) + ".body.out"
			print (cmd3)
			os.system(cmd3)
			
			d1 = {}
			d2 = {}
			d3 = {}
			body_len = {}
			
			stat1 = open(output + "/pair-" + str(a) + ".TSS.stat", 'a')
			stat2 = open(output + "/pair-" + str(a) + ".TTS.stat", 'a')
			stat3 = open(output + "/pair-" + str(a) + ".body.stat", 'a')

			for i in range(1,101):
				d1[i] = 0
				d2[i] = 0
				d3[i] = 0
				body_len[i] = 0

			with open(wk_dir + "/output/pair-" + str(a) + ".TSS.out") as f:
				for e in f:
					e = e.split("\t")
					d1[int(e[3].split("_")[2])] += int(e[4])
			with open(wk_dir + "/output/pair-" + str(a) + ".TTS.out") as f:
				for e in f:
					e = e.split("\t")
					d2[int(e[3].split("_")[2])] += int(e[4])
			with open(wk_dir + "/output/pair-" + str(a) + ".body.out") as f:
				for e in f:
					e = e.split("\t")
					d3[int(e[3].split("_")[2])] += int(e[4])
					body_len[int(e[3].split("_")[2])] += int(e[6])
			
			for i in range(1,101):
				stat1.write(str(i-100) + "\t" + str(d1[i] * 1000000000/(gene_num * read_num * 10)) + "\texpression level " + str(a) + "\n")
				stat3.write(str(i) + "\t" + str(((d3[i]/(body_len[i]/gene_num)) * 1000000000)/(gene_num * read_num)) + "\texpression level " + str(a) + "\n")
				stat2.write(str(i+100) + "\t" + str(d2[i] * 1000000000/(gene_num * read_num * 10)) + "\texpression level " + str(a) + "\n")

			stat1.close()
			stat2.close()
			stat3.close()

		cmd = "cat " + output + "/*.stat > " + output + "/all.results.txt"
		os.system(cmd)

		print ("----------------- calculate distribution processing okay .. -------------------")
		print ("time: " + str (time.time()-time_start))

def combine(input_file,output_file,ref,bam,divide):
#	sort_bam(bam)
#	genes_info_processing(input_file,ref,divide,bam)
	calculate_distribution(input_file,ref,divide,bam,output_file)

def main():
        parser = get_parser()
        args = vars(parser.parse_args())

        if args['input'] is not None:
                print(version())

        combine(args['input'], args['output'], args['gtf'], args['bam'], args['divide'])


if __name__ == "__main__":
	main()
