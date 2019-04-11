#generate sequences
import app.home.main3 as m3
from app.home.snp import Snp
from app.home.allele import Allele
import pandas as pd
import argparse

def gen_sequence(snp_names, info_list,genome_version):
	#parser = argparse.ArgumentParser(description='')
	#parser.add_argument('--snp_name', default="rs268", type=list)
	#parser.add_argument('--verify', default='n', type=str)
	#parser.add_argument('--delete_snp', default='n', type=str)
	#parser.add_argument('--return_list', default=False, type=bool)
	#parser.add_argument('--genome_version', default=GRCh37, type=str)
	#args = vars(parser.parse_args())
	snp_list_object = []
	minor_allele_list = []
	res = []

	# index = 0

	for i,snp_name in zip(info_list,snp_names):

		# snp_location = int(i[index]['location'])
		# print("Identifying problem: ", i[index]['chrom'])
		# snp_chrom = int(i[index]['chrom'])
		# snp_al = i[index]['allele_wt']

		#from table. Snp location is in first position
		snp_location = int(i[1])
		snp_chrom = int(i[2])
		snp_allele = i[3]

		list_minor = []

		#split allele and insert in a list
		minor_allele = i[4].split('|')
		for al in minor_allele:
			list_minor.append(al)

		#check if there is more than one minor allele
		# if (len(i) > 1):
		# 	for dic in i:
		# 		list_minor.append(dic['allele_v'])
		# else:
		# 	list_minor.append(i[index]['allele_v'])

		allele_comum_insert = Allele(nome=snp_allele,local=snp_location,cromossomo=snp_chrom,is_comum=True,snp_pos=0)

		for m_allele in list_minor:
			minor_allele_insert = Allele(nome=m_allele,local=snp_location,cromossomo=snp_chrom,is_comum=False,snp_pos=0)
			minor_allele_list.append(minor_allele_insert)

		#cria um objeto do tipo snp e insere as informações da snp nele
		snp_insert = Snp(name=snp_name, location=snp_location, chrom=snp_chrom,
				 charact="SNV",ancestral_al= allele_comum_insert,
				 minor_al=minor_allele_list)

		#coloca a snp na lista de snps
		snp_list_object.append(snp_insert)

		minor_allele_list = []
		list_minor = []

	args = {"snp_list":snp_list_object,
		"return_list":True,
		"verify":'n',
		"delete_snp":'n',
		"genome_version":genome_version}
	res.append(m3.main(args))
	
	#print(res)
	return res
