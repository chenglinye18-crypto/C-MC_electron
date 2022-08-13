#!/bin/bash

FILENAME="config.txt"
TEMPLATE_FILE="template.gnuplot"
WORKING_FILE="working.gnuplot"

read gridfile  
read fig_num 

echo $fig_num

for ((i=1 ; i <=$fig_num; i++))
do
  for ((j=1 ; j <= 4; j++))
  do
    read LINE
    echo $j, $LINE
  done
  read data_name
  echo $data_name
  read fig_name
  rm -f $WORKING_FILE
  echo "set output \"$fig_name\"" | cat >> $WORKING_FILE 
  cat $TEMPLATE_FILE >> $WORKING_FILE
  echo "splot \"$data_name\" u 1:2:3 ti 'Vg=Vd=1V, noqc'  with line palette" | cat >> $WORKING_FILE
  gnuplot $WORKING_FILE
done
rm -f $WORKING_FILE
