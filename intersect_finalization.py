#!/usr/bin/env python3

#script to get the gene (and genic region) for each mutation, identify mutations within a specified window of a TSS for potential promoter/enhancer mutations.

import argparse
import sys, os
import subprocess
import shutil
import pandas as pd
import numpy as np
import gffpandas.gffpandas as gffpd

parser = argparse.ArgumentParser(description="Script to clean bedtools intersect output")

parser.add_argument("-f","--file", type=str, help="Cleaned intersect output")
parser.add_argument("-g", "--gff", type=str, help="gff file")
parser.add_argument("-o", "--output", type=str, help="Prefix for output")
parser.add_argument("-w", "--window", type= , help="Window to determine if SNP is proximal to the gene or in a potential promoter region. Default is 1000")
parser.add_argument("-c", "--distance", type="store_true", help="Flag if you want distance from TSS")
parser.add_argument("-fc", "--fcount", type="store_true", help="Flag if you want number of flanking mutations per gene")

args = parser.parse_args()
file = args.file
gff = args.gff
out = args.output
window = args.output

gffd = gffpd.read_gff3(gff)
gene_gff = gffd.filter_feature_of_type(['gene']) #assumes all the genes are annotated as genes in the gff

df = pd.read_csv(file,sep="\t", header = None)

#drop unnecessary columns and rename them. This will have to be edited based on what you are working with. Just make sure that scaffold and locus remain.
df2 = df.drop(df.iloc[:,2:30],axis = 1) #this will need to be changed
df2.drop(columns=df2.columns[-2:], axis=1, inplace=True) #this assumes that the last two columns are the bed coordinates
coln = ["scaffold", "locus", "annotator","feature_type", "gff_start-ish", "gff_end-ish", "score", "strand", "frame", "attribute"] #columns that must be kept are scaffold, locus, feature_type, strand, and attribute. if these change, will need to change again at 263
df2.columns=coln

#Pull out the gene and scaffold-level SNPS
gene_raw = df2.loc[df2['feature_type'] == "gene"].reset_index()
gene_raw = gene_raw.iloc[:,1:]
contig_raw = df2.loc[df2['feature_type'] == "contig"].reset_index()
contig_raw = contig_raw.iloc[:,1:]

#Take out contig and gene from the looping feature set, mRNA, exon, and protein/esp matches can be removed for expediency as well.
# leaves behind CDS, 3' UTR and 5' UTR. Introns are calculated later.
looping = df2[df2.feature_type != 'mRNA']
looping =looping[looping.feature_type != 'contig']
looping =looping[looping.feature_type != 'exon']
looping = looping[looping.feature_type != 'gene']
looping = looping[looping.feature_type != 'expressed_sequence_match']
looping = looping[looping.feature_type != 'protein_match']
looping = looping[looping.feature_type != 'match_part']
looping = looping.reset_index()
looping = looping.iloc[:,1:]
looping

#runs in about 10-15 minutes. This creates a df for the geneic feature types. IDs introns by seeing what is in the 'gene_raw' df but hasn't already been assigned to CDS or UTR 
print("Identifying Geneic feature types")
three_prime_filtered = pd.DataFrame(columns = coln)
five_prime_filtered = pd.DataFrame(columns = coln)
cds_filtered = pd.DataFrame(columns = coln)
introns_filtered = pd.DataFrame(columns = coln)
for index,row in gene_raw.iterrows():
    #print(row["attribute"])
    n = 0
    for x,y in looping.iterrows():
        #print(row["attribute"])
        if (row["scaffold"] == y["scaffold"] and row["locus"] == y["locus"]):
            if (y["feature_type"] == "three_prime_UTR"):
                temp = looping[x:(x+1)]
                three_prime_filtered= pd.concat([three_prime_filtered, temp], ignore_index = True)
                n = 1
            elif (y["feature_type"] == "five_prime_UTR"):
                temp = looping[x:(x+1)]
                five_prime_filtered= pd.concat([five_prime_filtered, temp], ignore_index = True)
                n = 1
            elif (y["feature_type"] == "CDS"):
                temp = looping[x:(x+1)]
                cds_filtered= pd.concat([cds_filtered, temp], ignore_index = True)
                n = 1
    
    if n == 0:
        temp = gene_raw[index:(index+1)]
        temp.loc[temp.feature_type == "gene", 'feature_type'] = "intron"
        introns_filtered= pd.concat([introns_filtered, temp], ignore_index = True)
        
        
        
keep = pd.concat([five_prime_filtered, three_prime_filtered, cds_filtered, introns_filtered], ignore_index=True) #this is now the gene related mutations

#Identify the flanking region mutations by seeing what is in the 'contig_raw' df but hasn't already been kept.
print("Identifying flanking region mutations")
flanking_filtered  = contig_raw
for index,row in flanking_filtered.iterrows():
    for x,y in keep.iterrows():
        if (row["scaffold"] == y["scaffold"] and row["locus"] == y["locus"]):
            flanking_filtered=flanking_filtered.drop([index])  
            
            
# This takes the longest. This goes through the gff, and for each gene sees if the mutation is within X bps
if window == None:
	window = 1000000
print("Identifying proximal flanking genes, with window set to" + window + ". Go stretch.")

flanking_filtered["flanking_genes"] = 'Nothing'
flanking_filtered["percent_long"] = 'Nothing'
flanking_filtered["flank_count"] = 'Nothing'
for index,row in flanking_filtered.iterrows():
    n = 0
    flc = 0
    flanking_filtered.at[index,"Gene_Name"] = "Flanking"
    flanking_filtered.at[index, 'percentage'] = -1
    for x,y in gene_gff.df.iterrows():
        if row["scaffold"] == y["seq_id"]:
            if y["strand"] == "+":
                if (y["start"] - window) <= row["locus"]  <= (y["start"] + window):
                    if n == 1:
                        t1 = (flanking_filtered.at[index,'flanking_genes'] + "#" + y["attributes"])
                        flanking_filtered.at[index,'flanking_genes'] = t1
                        locusdiff = abs(row["locus"] - y["start"])
                        absdiff = abs(y["start"] - window)
                        percent = (locusdiff / absdiff)* 100
                        percent = str(percent)
                        t2 = (flanking_filtered.at[index,'percent_long'] + "#" + percent)
                        flanking_filtered.at[index,'percent_long'] = t2
                        locusdiff = str(locusdiff)
                        t3 = (flanking_filtered.at[index,'bp_long'] + "#" + locusdiff)
                        flanking_filtered.at[index,'bp_long'] = t3
                        flc = flc + 1
                    if n == 0:
                        n = 1
                        flanking_filtered.loc[index,'feature_type'] = "Proximal Flanking"
                        flanking_filtered.at[index,'flanking_genes'] = y["attributes"]
                        locusdiff = abs(row["locus"] - y["start"])
                        absdiff = abs(y["start"] - window)
                        percent = (locusdiff / absdiff)* 100
                        percent = str(percent)
                        flanking_filtered.at[index,'percent_long'] = percent
                        locusdiff = str(locusdiff)
                        flanking_filtered.at[index,'bp_long'] = locusdiff
                        flc = flc + 1
            if y["strand"] == "-":
                if (y["end"] - window) <= row["locus"]  <= (y["end"] + window):
                    if n == 1:
                        t1 = (flanking_filtered.at[index,'flanking_genes'] + "#" + y["attributes"])
                        flanking_filtered.at[index,'flanking_genes'] = t1
                        locusdiff = abs(row["locus"] - y["end"])
                        absdiff = abs(y["end"] - window)
                        percent = (locusdiff / absdiff)* 100
                        percent = str(percent)
                        t2 = (flanking_filtered.at[index,'percent_long'] + "#" + percent)
                        flanking_filtered.at[index,'percent_long'] = t2
                        locusdiff = str(locusdiff)
                        t3 = (flanking_filtered.at[index,'bp_long'] + "#" + locusdiff)
                        flanking_filtered.at[index,'bp_long'] = t3
                        flc = flc + 1
                    if n == 0:
                        n = 1
                        flanking_filtered.loc[index,'feature_type'] = "Proximal Flanking"
                        flanking_filtered.at[index,'flanking_genes'] = y["attributes"]
                        locusdiff = abs(row["locus"] - y["end"])
                        absdiff = abs(y["end"] - window)
                        percent = (locusdiff / absdiff)* 100
                        percent = str(percent)
                        flanking_filtered.at[index,'percent_long'] = percent
                        locusdiff = str(locusdiff)
                        flanking_filtered.at[index,'bp_long'] = locusdiff
                        flc = flc + 1
                
    flanking_filtered.loc[index,"flank_count"] = flc
    if n == 0:
        flanking_filtered.loc[index, 'feature_type'] = "Flanking_region"
    

    
flanking_filtered.loc[flanking_filtered["percentage"] == -1, "percentage"] = np.nan
flanking_filtered.loc[flanking_filtered["flanking_genes"] == 'Nothing', "flanking_genes"] = np.nan
flanking_filtered.loc[flanking_filtered["percent_long"] == 'Nothing', "percent_long"] = np.nan


# What % of the gene is the mutation found in? and Flanking?
print("Checking geneic mutations for flanking possibility. Go for another stretch")
keep["flanking_genes"] = 'Nothing'
keep["percent_long"] = 'Nothing'
keep["bp_long"] = 'Nothing'
keep["flank_count"] = 'Nothing'
for index,row in keep.iterrows():
    n = 0
    flc = 0
    for x,y in gene_gff.df.iterrows():
        if row["scaffold"] == y["seq_id"]:
            #first get % for the actual gene
            if y["strand"] == "+":
                if y["start"] <= row["locus"]  <= y["end"]:
                    keep.at[index, 'Gene_Name'] = y["attributes"]
                    absdiff = abs(y["start"] - y["end"])
                    locusdiff = abs(row["locus"] - y["start"])
                    percent = (locusdiff / absdiff)* 100
                    keep.at[index, 'percentage'] = percent
                    keep.at[index, 'TSS_dist'] = locusdiff
                if (y["start"] - window) <= row["locus"]  <= (y["start"] + window):
                    if n == 1:
                        t1 = (keep.at[index,'flanking_genes'] + "#" + y["attributes"])
                        keep.at[index,'flanking_genes'] = t1
                        locusdiff = abs(row["locus"] - y["start"])
                        absdiff = abs(y["start"] - window)
                        percent = (locusdiff / absdiff)* 100
                        percent = str(percent)
                        t2 = (keep.at[index,'percent_long'] + "#" + percent)
                        keep.at[index,'percent_long'] = t2
                        locusdiff = str(locusdiff)
                        t3 = (keep.at[index,'bp_long'] + "#" + locusdiff)
                        keep.at[index,'bp_long'] = t3
                        flc = flc + 1
                    if n == 0:
                        n = 1
                        keep.at[index,'flanking_genes'] = y["attributes"]
                        locusdiff = abs(row["locus"] - y["start"])
                        absdiff = abs(y["start"] - window)
                        percent = (locusdiff / absdiff)* 100
                        percent = str(percent)
                        keep.at[index,'percent_long'] = percent
                        locusdiff = str(locusdiff)
                        keep.at[index,'bp_long'] = locusdiff
                        flc = flc + 1
            if y["strand"] == "-":
                if y["end"] >= row["locus"]  >= y["start"]:
                    keep.at[index, 'Gene_Name'] = y["attributes"]
                    absdiff = abs(y["start"] - y["end"])
                    locusdiff = abs(row["locus"] - y["end"]) ##because neg strand
                    percent = (locusdiff / absdiff)* 100
                    keep.at[index, 'percentage'] = percent
                    keep.at[index, 'TSS_dist'] = locusdiff
            #Next check for flanking genes
                if (y["end"] - window) <= row["locus"]  <= (y["end"] + window):
                    if n == 1:
                        t1 = (keep.at[index,'flanking_genes'] + "#" + y["attributes"])
                        keep.at[index,'flanking_genes'] = t1
                        locusdiff = abs(row["locus"] - y["end"])
                        absdiff = abs(y["end"] - window)
                        percent = (locusdiff / absdiff)* 100
                        percent = str(percent)
                        t2 = (keep.at[index,'percent_long'] + "#" + percent)
                        keep.at[index,'percent_long'] = t2
                        locusdiff = str(locusdiff)
                        t3 = (keep.at[index,'bp_long'] + "#" + locusdiff)
                        keep.at[index,'bp_long'] = t3
                        flc = flc + 1
                    if n == 0:
                        n = 1
                        keep.at[index,'flanking_genes'] = y["attributes"]
                        locusdiff = abs(row["locus"] - y["end"])
                        absdiff = abs(y["end"] - window)
                        percent = (locusdiff / absdiff)* 100
                        percent = str(percent)
                        keep.at[index,'percent_long'] = percent
                        locusdiff = str(locusdiff)
                        keep.at[index,'bp_long'] = locusdiff
                        flc = flc + 1

    keep.loc[index,"flank_count"] = flc
    
keep = keep.loc[:,['scaffold','locus','annotator', 'feature_type', 'gff_start-ish','gff_end-ish','score','strand','frame','attribute','Gene_Name','TSS_dist', 'flank_count','percentage','flanking_genes','percent_long', 'bp_long']]
flanking_filtered["TSS_dist"] = np.nan
flanking_filtered = flanking_filtered.loc[:,['scaffold','locus','annotator', 'feature_type', 'gff_start-ish','gff_end-ish','score','strand','frame','attribute','Gene_Name','TSS_dist','flank_count','percentage','flanking_genes','percent_long','bp_long']]

#Write it out
name = out + "_gene.csv"
name2 = out + "_flanking.csv"
keep.to_csv(name, sep="\t", index = False)
flanking_filtered.to_csv(name2, sep="\t", index = False) 

#Whole Flank - Counts how many mutations are in the flanking region for each gene.
if args.fcount == True:
	print("Whole Flank Tally")
	alldf = pd.concat([keep,flanking_filtered], ignore_index=True)
	fll = []
	coord
	for x,y in gene_gff.df.iterrows():
		for index, row in alldf.iterrows():
			if y["strand"] == "+":
				if (y["start"] - window) <= row["locus"] <= y["start"] and row["scaffold"] == y["seq_id"]:
					fll.append(y["attributes"])
			if y["strand"] == "-":
				if y["end"] <= row["locus"] <= (y["end"] + window) and row["scaffold"] == y["seq_id"]:
					fll.append(y["attributes"])

	gene_list=[]
	counts = []
	for x in set(fll):
		gene_list.append(x)
		counts.append(fll.counts(x))
	
	flanking_whole = pd.DataFrame()
	flanking_whole["gene_list"] = gene_list
	flanking_whole["count"] = counts
	name = out + "_flanking_count.csv"
	flanking_whote.to_csv(name, index=False)

#Coordinates - for each mutation within the window, it will write how far upstream from the TSS it is. (or technically downstream for - genes? you get the point)
if args.distance == True:
	print("Calculating flanking location in relation to TSS")
	contig1 = contig_raw

	for index,row in contig1.iterrows():
		n = 0
		for x,y in keep.iterrows():
			if (row["scaffold"] == y["scaffold"] and row["locus"] == y["locus"]):
				contig1.at[index,"Gene_Name"] = y["Gene_Name"] #label within gene snps with gene name
				n = 1
		if n == 0:
			contig1.at[index,"Gene_Name"] =  "Flanking"

	gene_list=[]
	coord = []
	for index,row in gene_gff.df.iterrows():
		for x,y in contig1.iterrows():
			if row["attributes"] != y["Gene_Name"] and row["seq_id"] == y["scaffold"]: #this way geneic snps are not counted.
				if row["strand"] == "+":
					if (row["start"] - window) <= y["locus"] <= (row["start"] + window):
						gene_list.append(row["attributes"])
						spot = y["locus"] - row["start"] #a snp at 100 with a TSS at 150 would be -50 since it is upstream
						coord.append(spot)
				if row["strand"] == "-":
					if (row["end"] - window) <= y["locus"] <= (row["end"] + window):
						gene_list.append(row["attributes"])
						spot = row["end"] - y["locus"] #a snp at 150 with a TSS at 100 would be -50 since that is upstream
						coord.append(spot)


	flanking_coords = pd.DataFrame()
	flanking_coords["gene_list"] = gene_list
	flanking_coords["count"] = coord
	flanking_coords
	name = out + "flanking _coords.csv"
	flanking_coords.to_csv(name, sep="\t", index = False)
