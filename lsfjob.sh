#!/bin/bash

docker run --rm --privileged --name jimmy_lsfjob -v /mnt/es/ness/jimmy:/home -v /mnt/es/mammoth/jimmy:/mammoth combinelab/salmon:latest /home/quant.sh
