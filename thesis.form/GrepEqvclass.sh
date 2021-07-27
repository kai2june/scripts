#!/bin/bash

path=''
transcript=''
while getopts "p:t:" argv
do
    case "${argv}" in
        p)
            path=${OPTARG};
            ;;
        t)
            transcript=${OPTARG};
            ;;
        *) 
            echo "exception: getopts problem."
            ;;
    esac
done

transcript_id=`grep -n ${transcript} ${path} | awk -F":" '{print $1-3}'`
echo "transcript_id" ${transcript_id}
grep "\<${transcript_id}\>" ${path} | awk -F"\n" '{print $1,"\n"}'