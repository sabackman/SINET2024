#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 17:07:05 2024

@author: samuel
"""
import os

chr_arms = {}

ploidy=2

with open("chromosome_arms.txt", "r") as infile:
	infile.readline()
	for line in infile:
		words = line.strip().split("\t")
		chr_arms[words[0]] = [int(words[1]), int(words[2])]
		

samples = {}
for file in os.listdir("SINET_cf2"):
	if not "bam_ratio.txt" in file: continue
	if "png" in file: continue
	sample = file.split(".")[0]
	with open(os.path.join("SINET_cf2", file), "r") as infile:
		segs = {}
		for k in chr_arms.keys():
			segs[k] = {"loss": 0, "neutral": 0, "gain": 0}
		infile.readline()
		for line in infile:
			words = line.strip().split("\t")
			if words[0] == "X" or words[0] == "Y": continue
			if words[2] == "-1": continue #skip missing data			
			if int(words[1]) >= chr_arms[words[0]+"p"][0] and int(words[1]) <= chr_arms[words[0]+"p"][1]:
				seg = words[0]+"p"
			elif int(words[1]) >= chr_arms[words[0]+"q"][0] and int(words[1]) <= chr_arms[words[0]+"q"][1]:
				seg = words[0]+"q"
			else:
				print("ERROR")
			
			if int(words[4]) < ploidy:
				type_="loss"
			elif int(words[4]) == ploidy:
				type_="neutral"
			else:
				type_="gain"
				
			segs[seg][type_] += 1
			
			samples[sample] = segs
	
outfile = open("chr_arms_data.txt", "w+")
outfile.write("Sample\t"+"\t".join(chr_arms.keys()))
for sample in samples.keys():
	outline = [sample]
	for k in chr_arms.keys():
		type_ = "0"
		#if 50000*samples[sample][k]["loss"] > 0.8*(chr_arms[k][1]-chr_arms[k][0]):
		if samples[sample][k]["loss"] > 0.8*(samples[sample][k]["loss"]+samples[sample][k]["neutral"]+samples[sample][k]["gain"]):
			type_ = "-1"
		elif samples[sample][k]["gain"] > 0.8*(samples[sample][k]["loss"]+samples[sample][k]["neutral"]+samples[sample][k]["gain"]):
		#elif 50000*samples[sample][k]["gain"] > 0.8*(chr_arms[k][1]-chr_arms[k][0]):
			type_ = "1"
		outline.append(type_)
	outfile.write("\n"+"\t".join(outline))
outfile.close()		
			
		
		
