#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
#
#main class
########Make snp sequences
#
#Description: program receives as input snp name and size, process the information 
#and return as output sequences
# with snps in middle position
#com as combinações das snps se necessário. Está modificado para python3
#snps assembly: GRCh37
#
########Algorithm:
#1)Define variables
#2)Pegar informações da snp atraves do servidor e colocar a snp na lista
#3)Organizar a lista com base no local
#4)Fazer um loop na lista para ver se existe uma distancia maenor de 50 nts entre uma snp e outra
#5)Se não tiver, fazer as combinacões com uma snp e se tiver juntas snps e fazer as sequencias com as snps e as suas
#combinações
#
#
#
from app.home.snp import Snp
from app.home.allele import Allele
from app.home.graphDAO import GraphDAO
from app.home.snpDAO import SnpDAO
from app.home.sequencia import Sequence
from app.home.epigenome import Epigenome 
import pandas as pd
import numpy as np
import requests, sys
import argparse

#----------------------------------------------------Métodos-----------------------------------------------------///

#----------------------------------------------Programa começa aqui----------------------------------------------///

def main(args):

	#cria uma lista vazia de snps
	#lista_snp = []

	#declara o servidor
	#server = "http://rest.ensembl.org"

	#usuario define o tamanho da lista
	#lista_tamanho = int(input("Defina o tamanho da sua lista(25 ou 50): "))

	#Usuário insere quantas snps ele que fazer as sequencia
	#numero_de_snps = int(input("coloque com quantas snps você deseja fazer as sequencias: "))

	#lista de minor allele_list
	#minor_allele_list = []

	#----------------------Isso aqui dentro de um while para poder pegar as várias snps do usuario------------------///

	#for i in range(numero_de_snps):
	genome_version = args['genome_version']
	
	lista_snp = args['snp_list']

	
	#allele1 = Allele(nome='C',local=6704583,cromossomo=1,is_comum=True,snp_pos=0)
	#allele2 = Allele(nome='A',local=6704583,cromossomo=1,is_comum=False,snp_pos=0)
	#allele7 = Allele(nome='T',local=6704583,cromossomo=1,is_comum=False,snp_pos=0)
	#allele3 = Allele(nome='T',local=6704603,cromossomo=1,is_comum=True,snp_pos=0)
	#allele4 = Allele(nome='C',local=6704603,cromossomo=1,is_comum=False,snp_pos=0)
	#allele5 = Allele(nome='T',local=135536917,cromossomo=5,is_comum=True,snp_pos=0)
	#allele6 = Allele(nome='C',local=135536917,cromossomo=5,is_comum=False,snp_pos=0)
	#allele8 = Allele(nome='G',local=85215584,cromossomo=1,is_comum=True,snp_pos=0)
	#allele9 = Allele(nome='A',local=85215584,cromossomo=1,is_comum=False,snp_pos=0)
	#allele10 = Allele(nome='C',local=85215584,cromossomo=1,is_comum=False,snp_pos=0)

	#snp1 = Snp(name='rs278020',location=6704583,chrom=1,charact='snv',ancestral_al=allele1,
				#minor_al=[allele2,allele7])
	#snp2 = Snp(name='rs278021',location=6704603,chrom=1,charact='snv',ancestral_al=allele3,
				#minor_al=[allele4])
	#snp3 = Snp(name='rs123121',location=135536917,chrom=5,charact='snv',ancestral_al=allele5,
				#minor_al=[allele6])
	#snp4 = Snp(name='rs12117219',location=85215584,chrom=1,charact='snv',ancestral_al=allele8,
				#minor_al=[allele9,allele10])
	
	#lista_snp = [snp1,snp2,snp3,snp4]

	#---------------------sort list by local and chromossome-------------------------///

	# timsort on chromossome and location
	lista_snp.sort(key=lambda x: (x.chrom, x.location), reverse=False)

	#-----------------------------------------------------------------verificar o tamanho----------------------------------------------------///

	lista_comb = []

	snp_stuff = SnpDAO()

	first_write = True
	#create snps dataframes
	columns_df = ['Name',
				'Chromossome',
				'Location',
				'Characteristic',
				'Allele Wild Type',
				'Allele Variation']
	
	#insert snp information into dataframe
	snp_df = pd.DataFrame(([i.name,
							i.chrom,
							i.location,
							i.charact,
							i.ancestral_al,
							i.minor_al] for i in lista_snp), columns=columns_df)
	#make sequences
	#group snps by their chromossomes
	snp_df_2 = snp_df.groupby('Chromossome')
	#divide them in chromossome groups
	groups = dict(list(snp_df_2))
	#iterate chromossomes groups
	for chrom_id in groups:
		#get dataframe for whole group
		chrom_df = snp_df_2.get_group(chrom_id)
		#calculate difference between locations
		chrom_df['Difference'] = chrom_df['Location'].diff()
		#null values is zero
		chrom_df.fillna(0,inplace = True)
		#get size until the end of iterations
		size = len(chrom_df)-1
		#make sequence if there is only one element
		if (size == 0):
			snp = Snp(	name=chrom_df.iloc[0]['Name'],
						location=chrom_df.iloc[0]['Location'],
						chrom=chrom_df.iloc[0]['Chromossome'],
						charact=chrom_df.iloc[0]['Characteristic'],
						ancestral_al=chrom_df.iloc[0]['Allele Wild Type'],
						minor_al=chrom_df.iloc[0]['Allele Variation'])
			print("135 name: " + chrom_df.iloc[0]['Name'])
			first_write = snp_stuff.request_sequence(snp,genome_version,first_write)
		#check if there is more than one element
		elif (size > 0):
			#iterate chrmossome list
			for i in range(len(chrom_df)-1):
				next_smaller = chrom_df.iloc[i + 1]['Difference'] <50
				bigger = chrom_df.iloc[i]['Difference'] >=50
				current_0 = chrom_df.iloc[i]['Difference'] == 0
				smaller = chrom_df.iloc[i]['Difference'] <50
				#check if it needs to be in list of combinations
				if ( (bigger and next_smaller) or (smaller and not(current_0)) or (current_0 and next_smaller) ):
					
					snp = Snp(	name=chrom_df.iloc[i]['Name'],
								location=chrom_df.iloc[i]['Location'],
								chrom=chrom_df.iloc[i]['Chromossome'],
								charact=chrom_df.iloc[i]['Characteristic'],
								ancestral_al=chrom_df.iloc[i]['Allele Wild Type'],
								minor_al=chrom_df.iloc[i]['Allele Variation'])
					lista_comb.append(snp)
					#if its last iteration
					if i == size-1:
						if (next_smaller):
							snp = Snp(	name=chrom_df.iloc[i+1]['Name'],
										location=chrom_df.iloc[i+1]['Location'],
										chrom=chrom_df.iloc[i+1]['Chromossome'],
										charact=chrom_df.iloc[i+1]['Characteristic'],
										ancestral_al=chrom_df.iloc[i+1]['Allele Wild Type'],
										minor_al=chrom_df.iloc[i+1]['Allele Variation'])
							lista_comb.append(snp)
							first_write = snp_stuff.request_sequence_combinations(lista_comb,genome_version,first_write)
							#delete combinations list
							del lista_comb[:]
						#last iteration is not in combinations list
						else:
							snp = Snp(	name=chrom_df.iloc[i+1]['Name'],
										location=chrom_df.iloc[i+1]['Location'],
										chrom=chrom_df.iloc[i+1]['Chromossome'],
										charact=chrom_df.iloc[i+1]['Characteristic'],
										ancestral_al=chrom_df.iloc[i+1]['Allele Wild Type'],
										minor_al=chrom_df.iloc[i+1]['Allele Variation'])
							print("176 name: " + chrom_df.iloc[i+1]['Name'])
							first_write = snp_stuff.request_sequence(snp,genome_version,first_write)
							if (len(lista_comb) > 1 ):
								first_write = snp_stuff.request_sequence_combinations(lista_comb,genome_version,first_write)
                    			#delete combinations list
								del lista_comb[:]
				#if its not in combinations list make sequence
				else:
					#if there is more than one in combinations list make sequence
					if (len(lista_comb) > 1 ):
						print("------DEBUG------")
						first_write = snp_stuff.request_sequence_combinations(lista_comb,genome_version,first_write)
                    	#delete combinations list
						del lista_comb[:]
					snp = Snp(	name=chrom_df.iloc[i]['Name'],
								location= chrom_df.iloc[i]['Location'],
								chrom=chrom_df.iloc[i]['Chromossome'],
								charact=chrom_df.iloc[i]['Characteristic'],
								ancestral_al=chrom_df.iloc[i]['Allele Wild Type'],
								minor_al=chrom_df.iloc[i]['Allele Variation'])
					print("196 name: " + chrom_df.iloc[i]['Name'])
					first_write = snp_stuff.request_sequence(snp,genome_version,first_write)
					if i == size-1:
						snp = Snp(	name=chrom_df.iloc[i+1]['Name'],
									location= chrom_df.iloc[i+1]['Location'],
									chrom=chrom_df.iloc[i+1]['Chromossome'],
									charact=chrom_df.iloc[i+1]['Characteristic'],
									ancestral_al=chrom_df.iloc[i+1]['Allele Wild Type'],
									minor_al=chrom_df.iloc[i+1]['Allele Variation'])
						print("205 name: " + chrom_df.iloc[i+1]['Name'])
						first_write = snp_stuff.request_sequence(snp,genome_version,first_write)
	lista_comb = []

	return txt_2_list(args['return_list'], filename = "sequenciasdef.fna")
		
				
def txt_2_list(return_list, filename = "sequenciasdef.fna"):
  if(return_list):
    with open(filename,"r") as f:
      txt_to_list =[]
      for line in f:
        txt_to_list.append(line)
    f.close
    return txt_to_list
  else:
    return False		
	

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='')
	parser.add_argument('--snp_list', default=[], type=list)
	parser.add_argument('--verify', default='n', type=str)
	parser.add_argument('--delete_snp', default='n', type=str)
	parser.add_argument('--return_list', default=False, type=bool)
	parser.add_argument('--genome_version', default='GRCh37', type=str)
	args = vars(parser.parse_args())
	main(args)
