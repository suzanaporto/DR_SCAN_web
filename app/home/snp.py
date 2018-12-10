#!usr/local/bin/python3
# -*- coding: utf-8 -*-
#-----------Information about snp class-------------------///
#object class to maintain information about snps in the program
#
#@params
#name = <String>
#location = <int> //snp`s location in the genome
#chrom = <int> //chromossome
#charact = <string> // this attribute is needed to inform if snp is indel or snv
#ancestral_al = <allele >//
#minor_al = <allele []>//
#
#
class Snp(object):

    #object constructor
    def __init__(self, name, location, chrom, charact, ancestral_al, minor_al):
        self.name = name
        self.location = location
        self.chrom = chrom
        self.charact = charact
        self.ancestral_al = ancestral_al
        self.minor_al = minor_al

    def __repr__(self):
        return self.name
