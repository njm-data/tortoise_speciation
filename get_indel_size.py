#!/usr/bin/env python3

#script to get indel sizes

import argparse
import sys, os
import subprocess
import shutil
import pandas as pd
import gffpandas.gffpandas as gffpd

parser = argparse.ArgumentParser(description="Script to make indel sizes")

parser.add_argument("-f","--file", type=str, help="genes of interest output (all_limitedv2.csv)")
parser.add_argument("-i","--indels", type=str, help="VCF with all INDELs (All_INDELs.txt)")
parser.add_argument("-g","--gff", type=str, help="GFF file location")
parser.add_argument("-o","--output", type=str, help="Output file name")

args = parser.parse_args()
file = args.file
alli = args.indels
gff = args.gff

cleandf = pd.read_csv(file,sep=",")
raw_indels = pd.read_csv(alli, sep="\t", header = None)
df2 = raw_indels.drop(raw_indels.iloc[:,5:30],axis = 1)
df2.drop(columns=df2.columns[-2:], axis=1, inplace=True)
coln = ["scaffold", "locus","something", "reference", "alt", "annotator","feature_type", "gff_start-ish", "gff_end-ish", "score", "strand", "frame", "attribute"]
df2.columns=coln
#df2
test = cleandf.loc[cleandf["Mutation_type"] == "INDELs"].loc[cleandf["Count"] != 0].reset_index().iloc[:,2:]
#test
df3 = df2[df2.feature_type != 'mRNA']
df3 =df3[df3.feature_type != 'contig']
df3 =df3[df3.feature_type != 'exon']
df3 = df3[df3.feature_type != 'gene']
df3 = df3[df3.feature_type != 'expressed_sequence_match']
df3 = df3[df3.feature_type != 'protein_match']
df3 = df3[df3.feature_type != 'match_part']
df3 = df3.reset_index().iloc[:,1:]

tidy = gffpd.read_gff3(gff)
gene_gff = tidy.filter_feature_of_type(['gene'])
window = 1000

df4 = df2[df2.feature_type == 'contig']
#df4

tss_dist = []
for x in range(len(position)):
    tss_dist.append(None)
    
promoter_only = test[test.Feature_Type == "Promoter"]
for i2,r2 in gene_gff.df.iterrows():
    for index,row in promoter_only.iterrows():
        if row["Gene_Name"] == r2["attributes"]: # so we have the gene
            for i,r in df4.iterrows(): #now we go through the indel detail file
                if r2["attributes"] != r["attribute"] and r2["seq_id"] == r["scaffold"]:#avoid geneic indels but also make sure it's on the same scaffold
                    n = 0
                    if r2["strand"] == "+":
                        if (r2["start"] - window) <= r["locus"] <= r2["start"]:
                            for l in r["alt"]:
                                if l == ",":
                                    n = n + 1
                            if n == 0:
                                d = len(r["alt"]) - len(r["reference"])
                                difference.append(d)
                                long_gene.append(r["attribute"]) #this is done as a secondary check. make sure the long gene and the short gene match!
                                short_gene.append(row["gene"])
                                feature_type.append(row["Feature_Type"])
                                position.append(r["locus"])
                                scaffold.append(r["scaffold"])
                                tss_dist.append(r["locus"] - r2["start"])
                            elif n == 1:
                                pt1 = r["alt"].split(",")[0]
                                pt2 = r["alt"].split(",")[1]
                                d1 = len(pt1) - len(r["reference"])
                                d2 = len(pt2) - len(r["reference"])
                                d3 = (d1,d2)
                                difference.append(d3)
                                long_gene.append(r["attribute"]) #this is done as a secondary check. make sure the long gene and the short gene match!
                                short_gene.append(row["gene"])
                                feature_type.append(row["Feature_Type"])
                                position.append(r["locus"])
                                scaffold.append(r["scaffold"])
                                tss_dist.append(r["locus"] - r2["start"])
                            else:
                                print("more than two")
                                difference.append("more than two")
                                long_gene.append(r["attribute"]) #this is done as a secondary check. make sure the long gene and the short gene match!
                                short_gene.append(row["gene"])
                                feature_type.append(row["Feature_Type"])
                                position.append(r["locus"])
                                scaffold.append(r["scaffold"])
                                tss_dist.append(r["locus"] - r2["start"])
                    if r2["strand"] == "-":
                        if r2["end"] <= r["locus"] <= (r2["end"] + window) :
                            for l in r["alt"]:
                                if l == ",":
                                    n = n + 1
                            if n == 0:
                                d = len(r["alt"]) - len(r["reference"])
                                difference.append(d)
                                long_gene.append(r["attribute"]) #this is done as a secondary check. make sure the long gene and the short gene match!
                                short_gene.append(row["gene"])
                                feature_type.append(row["Feature_Type"])
                                position.append(r["locus"])
                                scaffold.append(r["scaffold"])
                                tss_dist.append(r2["end"] - r["locus"])
                            elif n == 1:
                                pt1 = r["alt"].split(",")[0]
                                pt2 = r["alt"].split(",")[1]
                                d1 = len(pt1) - len(r["reference"])
                                d2 = len(pt2) - len(r["reference"])
                                d3 = (d1,d2)
                                difference.append(d3)
                                long_gene.append(r["attribute"]) #this is done as a secondary check. make sure the long gene and the short gene match!
                                short_gene.append(row["gene"])
                                feature_type.append(row["Feature_Type"])
                                position.append(r["locus"])
                                scaffold.append(r["scaffold"])
                                tss_dist.append(r2["end"] - r["locus"])
                            else:
                                print("more than two")
                                difference.append("more than two")
                                long_gene.append(r["attribute"]) #this is done as a secondary check. make sure the long gene and the short gene match!
                                short_gene.append(row["gene"])
                                feature_type.append(row["Feature_Type"])
                                position.append(r["locus"])
                                scaffold.append(r["scaffold"])
                                tss_dist.append(r2["end"] - r["locus"])


indels_easy2 = pd.DataFrame()
indels_easy2["gene_name"] = long_gene
indels_easy2["gene"] = short_gene
indels_easy2["feature_type"] = feature_type
indels_easy2["size_change"] = difference
indels_easy2["scaffold"] = scaffold
indels_easy2["position"] = position
indels_easy2["tss_dist"] = tss_dist
indels_easy2

outputname = (args.output + ".csv")
indels_easy2.to_csv(outputname, sep="\t", index = False)
