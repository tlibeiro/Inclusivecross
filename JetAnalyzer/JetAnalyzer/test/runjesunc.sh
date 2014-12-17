#!/bin/bash

for n in `seq 6 1 10`; 
#for n in `seq 0 1 5`; 
do
	echo $n
  ./runjesunc $n >& jesunclog$n & 
done

