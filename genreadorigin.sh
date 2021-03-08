#!/bin/bash

awk '{print $4}' Drosophila_melanogaster.BDGP6.80.bed | awk '{split($0, a, ":"); print a[3]}' > read_origin.txt
