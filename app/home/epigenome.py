#!/usr/local/bin/python3
# -*- coding: utf-8 -*-

import requests, sys
import pybedtools
import gzip
import shutil
import itertools
import os
import os.path
from os import path

class Epigenome (object):

    def verificar_snps (self, tecido, snps):

            #E071
            url = 'https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/imputed12marks/jointModel/final/'+ tecido +'_25_imputed12marks_stateno.bed.gz'
            r = requests.get(url, allow_redirects=True)
            if r.status_code == 200:
                open('testando2.bed.gz', 'wb').write(r.content)
                
                with gzip.open('testando2.bed.gz', 'rb') as f_in:
                    with open('testando2.bed', 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                
                a = pybedtools.BedTool('testando2.bed')
            
            url2 = 'https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/'+ tecido +'_15_coreMarks_stateno.bed.gz'
            
            r1 = requests.get(url2, allow_redirects=True)
            if r1.status_code == 200:
                open('testando3.bed.gz', 'wb').write(r1.content)
                
                with gzip.open('testando3.bed.gz', 'rb') as f_in:
                    with open('testando3.bed', 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                b = pybedtools.BedTool('testando3.bed')
            
            
            url3 = 'https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/'+ tecido +'_18_core_K27ac_stateno.bed.gz'
            
            
            r2 = requests.get(url3, allow_redirects=True)
            if r2.status_code == 200:
                open('testando4.bed.gz', 'wb').write(r2.content)
                
                with gzip.open('testando4.bed.gz', 'rb') as f_in:
                    with open('testando4.bed', 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                c = pybedtools.BedTool('testando4.bed')
 
            state_15 = {
                         1 : "Active TSS",
                         2 : "Flanking Active TSS",
                         3 : "Transcr. at gene 5' and 3'",
                         4 : "Strong transcription",
                         5 : "Weak transcription",
                         6 : "Genic enhancers",
                         7 : "Enhancers",
                         8 : "ZNF genes & repeats",
                         9 : "Heterochromatin",
                         10 : "Bivalent/Poised TSS",
                         11 : "Flanking Bivalent TSS/Enh",
                         12 : "Bivalent Enhancer",
                         13 : "Repressed PolyComb",
                         14 : "Weak Repressed PolyComb",
                         15 : "Quiescent/Low",
                     }
            state_18 = {
                         1 : "Active TSS",
                         2 : "Flanking TSS",
                         3 : "Flanking TSS Upstream",
                         4 : "Flanking TSS Downstream",
                         5 : "Strong transcription",
                         6 : "Weak transcription",
                         7 : "Genic enhancer1",
                         8 : "Genic enhancer2",
                         9 : "Active Enhancer 1",
                         10 : "Active Enhancer 2",
                         11 : "Weak Enhancer",
                         12 : "ZNF genes & repeats",
                         13 : "Heterochromatin",
                         14 : "Bivalent/Poised TSS",
                         15 : "Bivalent Enhancer",
                         16 : "Repressed PolyComb",
                         17 : "Weak Repressed PolyComb",
                         18 : "Quiescent/Low",
                     }
            state_25 = {
                         1 : "Active TSS",
                         2 : "Promoter Upstream TSS",
                         3 : "Promoter Downstream TSS with DNase",
                         4 : "Promoter Downstream TSS",
                         5 : "Transcription 5'",
                         6 : "Transcription",
                         7 : "Transcription 3'",
                         8 : "Weak transcription",
                         9 : "Transcription Regulatory",
                         10 : "Transcription 5' Enhancer",
                         11 : "Transcription 3' Enhancer",
                         12 : "Transcription Weak Enhancer",
                         13 : "Active Enhancer 1",
                         14 : "Active Enhancer 2",
                         15 : "Active Enhancer Flank",
                         16 : "Weak Enhancer 1",
                         17 : "Weak Enhancer 2",
                         18 : "Enhancer Acetylation Only",
                         19 : "DNase only",
                         20 : "ZNF genes & repeats",
                         21 : "Heterochromatin",
                         22 : "Poised Promoter",
                         23 : "Bivalent Promoter",
                         24 : "Repressed PolyComb",
                         25 : "Quiescent/Low",
                     }
                     
                 
            #for i in a[2:10]:
            #    print (i)
            
            #25,15,18
            for o in snps:
                print ("nome da Snp: " + o.name)
                if 'a' in locals():
                    for m in a:
                        if m[0] == 'chr'+str(o.chrom)  and (  (int(m[1]) <= o.location and int(m[2]) >= o.location))  :
                            print (" | Elemento Regulatório do STATE MODEL 25 do tecido "+ tecido + " :"  + state_25[int(m[3])])
                if 'b' in locals():
                    for n in b:    
                        if n[0] == 'chr'+str(o.chrom)  and (  (int(n[1]) <= o.location and int(n[2]) >= o.location)):
                            print (" | Elemento Regulatório do STATE MODEL 15 do tecido "+ tecido + " :" + state_15[int(n[3])])
                if 'c' in locals():
                    for l in c:    
                        if l[0] == 'chr'+str(o.chrom)  and (  (int(l[1]) <= o.location and int(l[2]) >= o.location)):
                            print (" | Elemento Regulatório do STATE MODEL 18 do tecido "+ tecido + " :" + state_18[int(l[3])])

            if path.exists("testando2.bed"):
                os.remove('testando2.bed')
                os.remove('testando2.bed.gz')
            if path.exists("testando3.bed"):
                os.remove('testando3.bed')
                os.remove('testando3.bed.gz')
            if path.exists("testando4.bed"):
                os.remove('testando4.bed')
                os.remove('testando4.bed.gz')
            