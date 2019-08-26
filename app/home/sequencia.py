#!usr/local/bin/python3
# -*- coding: utf-8 -*-
#---------
#
#@params
#chrom
#local
#

class Sequence (object):

    #object constructor
    def __init__(self,seq,chrom,alelles,start,end,snps_nomes,snps_pos):
        self.seq = seq
        self.chrom = chrom
        self.alleles = alelles
        self.start = start
        self.end = end
        self.snps_nomes = snps_nomes
        self.snps_pos = snps_pos

        