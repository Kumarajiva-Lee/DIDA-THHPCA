#!/bin/bash 
namelist=$1
atm_group=`grep -i atm_group $namelist |grep -v num| awk '{print $3}' | sed $'s/\'//g'`
echo $atm_group

#atm_group_line=`cat $namelist | grep -n atm_group |grep -v num| awk -F ":" '{print $1}'`
atm_group_line=`grep -n atm_group $namelist |grep -v num| awk -F ":" '{print $1}'`
echo $atm_group_line

sed -i "${atm_group_line}s/$atm_group/aaa/g" $namelist


#sed -i "${atm_group_line}s/$atm_group/8888/g" $namelist


#`grep -i atm_group $namelist |grep -v num| sed -i "s/$atm_group/333/g" $namelist`
#echo $atm_group


#atm_group=`grep -i atm_group $namelist |grep -v num| awk '{print $3}' | sed $'s/\'//g'`
#echo $atm_group
