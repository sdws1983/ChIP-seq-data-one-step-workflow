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
        parser.add_argument('-o', '--output', help='output name', type=str)
        
        return parser


def make_annotation(ref, output):
	time_start = time.time()
	print ("----------------- first step : make annotation -------------------")
	
	rscript = open(output + "/Rscript.r", "a")
	command1 = 'library(GenomicFeatures)\nmaize.db <- makeTxDbFromGFF("' + ref + '",format = "gtf")\n'
	command2 = 'a<-promoters(maize.db,upstream = 3000,downstream = 0)\nb<-as.data.frame(fiveUTRsByTranscript(maize.db))\nc<-as.data.frame(threeUTRsByTranscript(maize.db))\nd<-as.data.frame(exonicParts(maize.db))\ne<-as.data.frame(intronicParts(maize.db))\n'
	command3 = 'write.table(a, file="' + output + "/promoters.xls" + '", sep = "\\t",col.names =F, row.names=F,quote=F)\n' + 'write.table(b[3:10], file="' + output + "/5utr.xls" + '", sep = "\\t",col.names =F, row.names=F,quote=F)\n' + 'write.table(c[3:10], file="' + output + "/3utr.xls" + '", sep = "\\t", col.names =F, row.names=F,quote=F)\n' + 'write.table(d[1:5], file="' + output + "/exons.xls" + '", sep = "\\t", col.names =F, row.names=F,quote=F)\n' + 'write.table(e[1:5], file="' + output + "/introns.xls" + '", sep = "\\t", col.names =F, row.names=F,quote=F)\n'
	command4 = "Rscript " + output + "/Rscript.r" + " && cut -f1-7 " + output + "/3utr.xls" + " |awk '" + '{if($5=="+"){print $1"\\t"$3+1"\\t"$3+3000"\\t"$4"\\t"$5"\\t"$6"\\t"$7}else{if($5=="-"){print $1"\\t"$2-3000"\\t"$2-1"\\t"$4"\\t"$5"\\t"$6"\\t"$7}' + "else{print}}}' > " + output + "/downstream.xls" 
	rscript.write(command1)
	rscript.write(command2)
	rscript.write(command3)
	rscript.close()

	os.system(command4)


	print ("----------------- make annotation done .. -------------------")
	print ("time: " + str (time.time()-time_start))
					
def annotation(input_file,ref,output):
	time_start = time.time()
	print ("----------------- second step : annotation processing -------------------")
	
	command1 = "intersectBed -a " + output + "/promoters.xls -b " + input_file + " -wb |cut -f11|sort|uniq > " + output + "/promoters.peak && grep -w -v -f " + output + "/promoters.peak " + input_file + " > " + output + "/1.bed"
	command2 = "intersectBed -a " + output + "/5utr.xls -b " + output + "/1.bed" + " -wb |cut -f12|sort|uniq > " + output + "/5utr.peak && grep -w -v -f " + output + "/5utr.peak " + output + "/1.bed" + " > " + output + "/2.bed"
	command3 = "intersectBed -a " + output + "/3utr.xls -b " + output + "/2.bed" + " -wb |cut -f12|sort|uniq > " + output + "/3utr.peak && grep -w -v -f " + output + "/3utr.peak " + output + "/2.bed" + " > " + output + "/3.bed"
	command4 = "intersectBed -a " + output + "/exons.xls -b " + output + "/3.bed" + " -wb |cut -f9|sort|uniq > " + output + "/exons.peak && grep -w -v -f " + output + "/exons.peak " + output + "/3.bed" + " > " + output + "/4.bed"
	command5 = "intersectBed -a " + output + "/introns.xls -b " + output + "/4.bed" + " -wb |cut -f9|sort|uniq > " + output + "/introns.peak && grep -w -v -f " + output + "/introns.peak " + output + "/4.bed" + " > " + output + "/5.bed"
	command6 = "intersectBed -a " + output + "/downstream.xls -b " + output + "/5.bed" + " -wb |cut -f11|sort|uniq > " + output + "/downstream.peak && grep -w -v -f " + output + "/downstream.peak " + output + "/5.bed" + " > " + output + "/Intergenic.peak"
	
	print (command1)
	os.system(command1)
	print (command2)
	os.system(command2)
	print (command3)
	os.system(command3)
	print (command4)
	os.system(command4)
	print (command5)
	os.system(command5)
	print (command6)
	os.system(command6)

	print ("----------------- annotation done .. -------------------")
	print ("time: " + str (time.time()-time_start))


def combine(input_file,output_file,ref):
	make_annotation(ref,output_file)
	annotation(input_file,ref,output_file)
def main():
        parser = get_parser()
        args = vars(parser.parse_args())

        if args['input'] is not None:
                print(version())

        combine(args['input'], args['output'], args['gtf'])


if __name__ == "__main__":
	main()
