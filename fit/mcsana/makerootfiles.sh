#!/bin/bash
python read_lst.py $1 $2
root -b -q "mcs_make_rootfilev0.C((char*)\"$2.csv\",(char*)\"$2.v0\")"
root -b -q "mcs_make_rootfilev1.C((char*)\"$2.v0\",(char*)\"$2\")"
