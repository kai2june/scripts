#!/bin/bash

path=''
number=''
while getopts "p:n:" argv
do
    case "${argv}" in
        p)
            path=${OPTARG}
            ;;
        n)
            number=${OPTARG}
            ;;
        *)
            echo "exception of getopts"
            ;;
    esac
done

head -n $((${number}+3)) ${path} | tail -n 1