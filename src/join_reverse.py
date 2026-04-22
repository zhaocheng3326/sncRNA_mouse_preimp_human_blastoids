#!/usr/bin/python3

import os
import sys

nd=[]
with open (sys.argv[1],"r") as not_include_file:
    for line in not_include_file:
        nd.append(line.strip("\n"))
for line in sys.stdin:
    ID=line.strip("\n").split("\t")[0]
    if ID not in nd:
        print(line.strip("\n"))

