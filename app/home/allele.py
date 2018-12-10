#!usr/local/bin/python3
# -*- coding: utf-8 -*-
#
#@params
#nome = <String>
#local = <int>
#cromossomo = <int>
#is_comum = <Boolean>
#snp_pos = <int>
#
#
class Allele(object):

    #object constructor
    def __init__(self, nome,local,cromossomo,is_comum,snp_pos):
        self.nome = nome
        self.local = local
        self.cromossomo = cromossomo
        self.is_comum = is_comum
        self.snp_pos = snp_pos

    def __repr__(self):
        return self.nome
