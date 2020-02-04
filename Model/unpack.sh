#!/bin/bash

for a in *.txz
do
    a_dir=`expr $a : '\(.*\).txz'`
    mkdir $a_dir 2>/dev/null
    tar -xvf $a -C $a_dir
done

