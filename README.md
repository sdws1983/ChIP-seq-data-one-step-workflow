# ChIP-seq-data-one-step-workflow

### A self-built workflow for processing ChIP-seq data.
```
usage: Reads_distribution_calculation.py [-h] [-i INPUT] [-f GTF] [-b BAM]
                                         [-d DIVIDE] [-o OUTPUT]

Environment : 

Python 3+, pandas, numpy, pysam
Bedtools, Samtools


optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        gene info
  -f GTF, --gtf GTF     gene annotation file (gtf)
  -b BAM, --bam BAM     bam file
  -d DIVIDE, --divide DIVIDE
                        divide genes into
  -o OUTPUT, --output OUTPUT
                        output name


usage: Peak_annotation.py [-h] [-i INPUT] [-f GTF] [-o OUTPUT]


Environment : 

Python 3+, R 
library("GenomicFeatures")
Bedtools


optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        gene info
  -f GTF, --gtf GTF     gene annotation file (gtf)
  -o OUTPUT, --output OUTPUT
                        output name

```

  
