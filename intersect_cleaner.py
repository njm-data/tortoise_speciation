#!/usr/bin/env python3

#script to clean up bedtools intersect results, remove duplicates

import argparse
import sys, os
import subprocess
import shutil
import pandas as pd

parser = argparse.ArgumentParser(description="Script to clean bedtools intersect")

parser.add_argument("-f","--file", type=str, help="vcf to clean")

args = parser.parse_args()
file = args.file

if file ==None:
	print("Missing input file")
	quit()
	
og = pd.read_csv(file, sep="\t", header=None)

#remove sliding windows so that they're all the same
mapping = {og.columns[0]:"scaffold", og.columns[1]: "locus"}
df = og.rename(columns=mapping)

#drop duplicates, write out new master
df = df.drop_duplicates(subset = ["scaffold","locus"],keep = 'last').reset_index()
df = df.iloc[:,1:]
newname = file+"_NODUPS.txt"
df.to_csv(newname, sep="\t", header = False, index = False)
