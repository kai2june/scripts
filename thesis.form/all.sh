#!/bin/bash

path=''
transcript=''
number=''
pro=''
while getopts "p:t:r:n:" argv
do
    case "${argv}" in
        p)
            path=${OPTARG};
            ;;
        t)
            transcript=${OPTARG};
            ;;
        r)
            pro=${OPTARG}
            ;;
        n)
            number=${OPTARG}
            ;;
        *) 
            echo "exception: getopts problem."
            ;;
    esac
done

transcript_id=`grep -n ${transcript} ${path} | awk -F":" '{print $1-3}'`
echo "transcript_id" ${transcript_id}
# grep "\<${transcript_id}\>" ${path}
new_transcript=`awk "NR==$((${number}+3)){print;}" ${path}`
grep ${new_transcript} ${pro}