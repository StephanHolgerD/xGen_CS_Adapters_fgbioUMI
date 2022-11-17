#!/bin/bash
samtools view -F 4 $1 | wc -l > $2
samtools view -f 2048 $1 | wc -l > $3
