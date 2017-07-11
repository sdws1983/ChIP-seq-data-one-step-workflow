#!/bin/bash
set -e

for ((i=1; i<=10; i ++))
do
{
A=$i".txt"
python2 ChIP-sub.py -i $A -o $1 -f $2 -s $3 -d $4
} &
done
wait  ##等待所有子后台进程结束
