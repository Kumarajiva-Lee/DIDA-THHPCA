#/bin/bash
string="Hello,Hi,How are you"
array=(${string//,/ })
echo ${array[@]}
for i in ${array[@]}
do
    echo $i
done
