from app.home import blueprint
from flask import render_template
from flask_login import login_required
from flask import request

import requests, sys
import numpy as np
import pandas as pd
import re
from flask import jsonify
from flask import request
import pybedtools
import gzip
import shutil
import itertools
import os
import os.path
from os import path
import json
<<<<<<< HEAD
import app.home.sequence_gen as sg
=======
>>>>>>> 11f9506f8fd474fb447f3b18b6eb81739dcc6a49

#step_1---------function
def get_snp_info(snp_id):
    ### Download Snp Info based on its id
    def request_info_by_id(snp_id):
        server = "https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/"

        r = requests.get(server+snp_id, headers={ "Content-Type" : "application/json"})

        if not r.ok:
            r.raise_for_status()
            sys.exit()

        decoded = r.json()

        return decoded

    ### Parse information about the SNP requested
    def parse_json(res_json):
        #get num of genomic placements (versions)
        num_genomic_placements = len(res_json['primary_snapshot_data']['placements_with_allele'])  

        # sometimes there is empty fields in the genomic placements... Detect and remove them! 
        temp_assertion = np.array([res_json['primary_snapshot_data']['placements_with_allele'][idx]['placement_annot']['seq_id_traits_by_assembly'] for idx in range(num_genomic_placements)])
        ban_idexs = np.array([len(x)!=0 for x in temp_assertion])
        temp_assertion = temp_assertion[ban_idexs]

        # Get the genomic placements (versions) names
        gnenome_versions = [ x[0]["assembly_name"] for x in temp_assertion]

        # get the number of variations
        allele_variations = len([res_json['primary_snapshot_data']['placements_with_allele'][idx]['alleles'] for idx in range(1)][0])

        # get the snps idexes
        snp_idxs = np.array([[res_json['primary_snapshot_data']['placements_with_allele'][x]['alleles'][y]['hgvs'] for y in range(allele_variations)][1:] for x in range(num_genomic_placements)])

        # remove those without valid genomic placements (versions) names 
        snp_idxs = snp_idxs[ban_idexs].tolist()

        return {"gnenome_versions": gnenome_versions, "snp_idxs_list": snp_idxs}

    ### Parse snp information
    def snp_info_input(snp_id):
        base = re.compile("[^(\d)\w+]").split(snp_id)[3]
        base_chrom = re.compile("[^(\d)\w+]").split(snp_id)[0]

        #  print (base)
        dict = {
            #"chrom" : re.compile("(\d+).0").split(base_chrom)[2],
            "chrom" : re.compile(".0{2,}").split(base_chrom)[1],
            "location" : re.compile("[^(\d)]").split(base)[0],
            "allele_wt" : re.compile("[(\d)]").split(base)[-1],
            "allele_v": re.compile("[^(\d)\w+]").split(snp_id)[4]
        }
        return dict

    ### Download Snp Info based on its id
    res_json = request_info_by_id(snp_id)
        
    ### Parse information about the SNP requested
    snp_info = parse_json(res_json)

    df_snp_info = pd.DataFrame(snp_info)
        
    gn_version = []

    for snp_idx_list in df_snp_info['snp_idxs_list']:
        res = []
        for snp_idx in snp_idx_list:
            res.append(snp_info_input(snp_idx))  
        gn_version.append(res)

    df_snp_info['snp_info_dict'] = gn_version  
            
    return gn_version	

#step_2---------function
def verificar_snps (tecido, snps,snp_name,chrom):

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
        cont = 1
        #for o in snps:
        print ("nome da Snp: " + snp_name)
        if 'a' in locals():
          for m in a:
            if ((m[0] == 'chr'+str(chrom))  and (  (int(m[1]) <= int(snps['Location']) and int(m[2]) >= int(snps['Location'])))):
              print (" | Elemento Regulatório do STATE MODEL 25 do tecido "+ tecido + " :"  + state_25[int(m[3])])
              row25 = {'snp_name':snp_name,'state_model':str(25),'tissue':tecido,'reg_elemnt':state_25[int(m[3])]}
        if 'b' in locals():
          for n in b:    
            if (n[0] == 'chr'+str(chrom)  and (  (int(n[1]) <= int(snps['Location']) and int(n[2]) >= int(snps['Location'])))):
              print (" | Elemento Regulatório do STATE MODEL 15 do tecido "+ tecido + " :" + state_15[int(n[3])])
              row15 = {'snp_name':snp_name,'state_model':str(15),'tissue':tecido,'reg_elemnt':state_15[int(n[3])]}       
        if 'c' in locals():
          for l in c:    
            if (l[0] == 'chr'+str(chrom)  and (  (int(l[1]) <= int(snps['Location']) and int(l[2]) >= int(snps['Location'])))):
              print (" | Elemento Regulatório do STATE MODEL 18 do tecido "+ tecido + " :" + state_18[int(l[3])])
              row18 = {'snp_name':snp_name,'state_model':str(18),'tissue':tecido,'reg_elemnt':state_25[int(l[3])]}

        if path.exists("testando2.bed"):
          os.remove('testando2.bed')
          os.remove('testando2.bed.gz')
        if path.exists("testando3.bed"):
          os.remove('testando3.bed')
          os.remove('testando3.bed.gz')
        if path.exists("testando4.bed"):
          os.remove('testando4.bed')
          os.remove('testando4.bed.gz')

        return row25,row15,row18

@blueprint.route('/index',methods=['GET','POST'])
@login_required
def index():
    return render_template('index.html')

@blueprint.route('/<template>')
@login_required
def route_template(template):
    return render_template(template + '.html')

@blueprint.route('/get_snp_info',methods=['GET','POST'])
@login_required
def teste():
    #if request method is post execute
    if request.method == 'POST':
        #get value from snp_name field from html form
        snp_input_form = str(request.form['snp_name'])
<<<<<<< HEAD
        genome_input = str(request.form['genome_v'])
        snp = snp_input_form[2:]
        print(snp)
        a = get_snp_info(snp)
        print(a)
        print(a[0][0])
        if "GRCh38" in genome_input:
            snp_info = a[0][0]
        if "GRCh37" in genome_input:
            snp_info = a[1][0]
        return jsonify(snp_info,snp_input_form)
=======
        snp = snp_input_form[2:]
        print(snp)
        a = get_snp_info(snp)
        print(a[0][0])
        return jsonify(a[0][0],snp_input_form)
>>>>>>> 11f9506f8fd474fb447f3b18b6eb81739dcc6a49
        
@blueprint.route('/verify_snps',methods=['GET','POST'])
@login_required
def verify_snps():
    dict_snps = []
    #get snp list information from the datatable
    snp_list_rows = request.form['snp_list']
    #transform it in a json file(dictionary?)
    snp_list_rows = json.loads(snp_list_rows)
    print(snp_list_rows) 
    if request.method == 'POST':
        tissue_list = request.form['tissue_field'].split("|")
        #iterate through tissues in tissues list
        for tissue in tissue_list:
            #iterate through every snp
            for snp_info in snp_list_rows:
                print(snp_info)
                #change chromossome type from number to 'X' or 'Y' string
                if str(snp_info['Chromossome']) == '23':
                    chrom = 'X'
                elif str(snp_info['Chromossome']) == '24':
                    chrom = 'Y'
                else:
                    chrom = snp_info['Chromossome']
                #TODO: change genome version for analysis (already done in the command line version)
                dict_snps += verificar_snps(tissue, snp_info, snp_info['Snp Name'], chrom)
<<<<<<< HEAD
    return jsonify(dict_snps)

@blueprint.route('/gen_sequence',methods=['GET','POST'])
@login_required
def gen_sequence():

    def get_snp_info(snp_id):
        ### Download Snp Info based on its id
        def request_info_by_id(snp_id):
            server = "https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/"

            r = requests.get(server+snp_id, headers={ "Content-Type" : "application/json"})

            if not r.ok:
                r.raise_for_status()
                sys.exit()

            decoded = r.json()

            return decoded

        ### Parse information about the SNP requested
        def parse_json(res_json):
            #get num of genomic placements (versions)
            num_genomic_placements = len(res_json['primary_snapshot_data']['placements_with_allele'])  

            # sometimes there is empty fields in the genomic placements... Detect and remove them! 
            temp_assertion = np.array([res_json['primary_snapshot_data']['placements_with_allele'][idx]['placement_annot']['seq_id_traits_by_assembly'] for idx in range(num_genomic_placements)])
            ban_idexs = np.array([len(x)!=0 for x in temp_assertion])
            temp_assertion = temp_assertion[ban_idexs]

            # Get the genomic placements (versions) names
            gnenome_versions = [ x[0]["assembly_name"] for x in temp_assertion]

            # get the number of variations
            allele_variations = len([res_json['primary_snapshot_data']['placements_with_allele'][idx]['alleles'] for idx in range(1)][0])

            # get the snps idexes
            snp_idxs = np.array([[res_json['primary_snapshot_data']['placements_with_allele'][x]['alleles'][y]['hgvs'] for y in range(allele_variations)][1:] for x in range(num_genomic_placements)])

            # remove those without valid genomic placements (versions) names 
            snp_idxs = snp_idxs[ban_idexs].tolist()

            return {"gnenome_versions": gnenome_versions, "snp_idxs_list": snp_idxs}

        ### Parse snp information
        def snp_info_input(snp_id):
            base = re.compile("[^(\d)\w+]").split(snp_id)[3]
            base_chrom = re.compile("[^(\d)\w+]").split(snp_id)[0]

            #  print (base)
            dict = {
                #"chrom" : re.compile("(\d+).0").split(base_chrom)[2],
                "chrom" : re.compile(".0{2,}").split(base_chrom)[1],
                "location" : re.compile("[^(\d)]").split(base)[0],
                "allele_wt" : re.compile("[(\d)]").split(base)[-1],
                "allele_v": re.compile("[^(\d)\w+]").split(snp_id)[4]
            }
            return dict

        ### Download Snp Info based on its id
        res_json = request_info_by_id(snp_id)
            
        ### Parse information about the SNP requested
        snp_info = parse_json(res_json)

        df_snp_info = pd.DataFrame(snp_info)
            
        gn_version = []

        for snp_idx_list in df_snp_info['snp_idxs_list']:
            res = []
            for snp_idx in snp_idx_list:
                res.append(snp_info_input(snp_idx))  
            gn_version.append(res)

        df_snp_info['snp_info_dict'] = gn_version  
                
        return df_snp_info	

    def create_alleles_dictionary(snp_names,info_list):
        dict_snp_allele = {}
        tuple_al = []
        
        index = 0
        #iterate through names and all dictionaies in snp
        for i,snp_name in zip(info_list,snp_names):
            snp_al = i[index]['allele_wt']
            tuple_al.append(snp_al)
            list_minor = []
            
            #check if there is more than one minor allele
            if (len(i) > 1):
                for dic in i:
                    #list_minor.append(dic['allele_v'])
                    tuple_al.append(dic['allele_v'])
                    print("nome da snp",snp_name)
            else:
                #list_minor.append(i[index]['allele_v'])
                tuple_al.append(i[index]['allele_v'])
            
        list_minor=[]
        
        print (tuple_al)
        
        dict_snp_allele[snp_name] = tuple_al 
        tuple_al = []
        
        return dict_snp_allele

    #if request method is post execute
    if request.method == 'POST':
        #get snp list information from the datatable
        filtered_snp = json.loads(request.form['snp_list'])
        #TODO: receive as a parameter
        gnenome_version = 'GRCh37.p13'

        #get snps ids from filtered list (stage 2)
        snp_ids_list = []
        for x in filtered_snp:
            snp_ids_list.append(x['Snp Name'])
        
        #only unique values
        snp_ids_list_unique = list(set(snp_ids_list))

        #download snp info
        info_snp_list = []
        for snp_id in snp_ids_list_unique:
            snp_id = snp_id[2:]
            info_snp_list.append(get_snp_info(snp_id))

        snp_names = []
        snp_dicts = []
        for snp_id,df_snp_info in zip(snp_ids_list_unique,info_snp_list):
        
            snp_name = snp_id
            sample_dict = df_snp_info[df_snp_info["gnenome_versions"] == gnenome_version]['snp_info_dict'].values[0]
            
            snp_names.append(snp_name)
            snp_dicts.append(sample_dict)

        genome_ver = gnenome_version.split(".")[0]
        res = sg.gen_sequence(snp_names, snp_dicts,genome_ver)
        # create dictionary with alleles when creating sequences
        dictionary_snp_allele = create_alleles_dictionary(snp_names,snp_dicts)

        #write list in txt
        fna_filename = './seq.fna'
        with open(fna_filename, 'w') as f:
            #for line in list_lines_file:
            for line in res:
                for element in line:
                    f.write("%s\n" % element)
        
        meme= "./meme.meme"
        os.system("export PATH=$HOME/meme/bin:$PATH ")
        os.system("~/meme/bin/fimo "+meme+" "+fna_filename )                    
    
    return jsonify(res)

@blueprint.route('/dif_tf',methods=['GET','POST'])
@login_required
def dif_tf():

    def row_out_of_range_single(seq_start,seq_stop,seq_name,idx):
        relative_position =  int(seq_name.split('|')[4])
        condition_3 = relative_position < seq_start 
        condition_2 = relative_position > seq_stop
        #index list
        index_list = []

        if  (condition_3 or condition_2) :
            return idx
        else:
            return None

    def row_out_of_range_multiple(seq_start,seq_stop,seq_name,idx):
        #check outside portion of positions
        seq_name_list =  seq_name.split('|')
        #index list
        index_list = []
        verify = False
        count = 0
        #count how many snps there are in the sequence
        for element in seq_name_list:
            if element[0:2] == 'rs':
                #print ( "Elemento: ",element )
                count += 1
        
        #starting position of relative positions in seq_name
        start_relative_list = (count*2) + 2
  
        relative_list = seq_name.split('|')[start_relative_list:]
  
        #print ("Relative list",relative_list)
  
        for position in relative_list:
    
            condition_3 = int(position) > seq_start 
            condition_2 = int(position) < seq_stop
            #condition tests if there is at least one snp in start-stop range 
            if (condition_3 and condition_2) :
                verify = True
        if verify:
            return None
        else:
            return idx

    def filter_range(table_output):
        print("Tamanho:",len(table_output))
        #table into data frame
        df_fimo_output = pd.read_csv(table_output, sep='\t')
        df_fimo_output = df_fimo_output.dropna()
        #start values into array
        array_start = df_fimo_output['start'].values
        #stop values into array
        array_stop = df_fimo_output['stop'].values
        #sequence name values
        array_sequence_name = df_fimo_output['sequence_name'].values
        #index list
        index_list = []
  
        #with pd.option_context('display.max_rows', 100, 'display.max_columns', 100):
        #display(df_fimo_output)
  
        for idx, (seq_start,seq_stop,seq_name) in enumerate( zip(array_start, array_stop, array_sequence_name) ):
    
            seq_type = seq_name.split('|')[0]
            condition_1 = seq_type == 'sequence_variation'
            condition_2 = seq_type == 'sequence_wild_type'
            condition_3 = seq_type == 'sequence_combinations'
    
            #first filter is out_of_range filter
            if (condition_1 or condition_2):
                index = row_out_of_range_single(seq_start,seq_stop,seq_name,idx)
                if (not(index == None)):
                    index_list.append(index)
            if (condition_3):
                index = row_out_of_range_multiple(seq_start,seq_stop,seq_name,idx)
                if (not(index == None)):
                    index_list.append(index)
      
        #print (str(index_list))
        print ("Tamanho da lista: " + str(len(index_list)))
        #print ("-------------------------Testando----------------------------")
        #filter by id
        #new_df = df_fimo_output.drop(df_fimo_output.index[index_list],axis=0)
        idx_list_selected = [ i for i in range(len(df_fimo_output)) if not i in index_list ]
  
        new_df = df_fimo_output.iloc[idx_list_selected]
  
        return new_df

    def filter_dataframe(table_output, log=False):
  
        new_df = filter_range(table_output)
  
        print ("Tamanho linha: ", len(new_df))  
        print ("DATAFRAME WITH RANGE FILTER")
        #display(new_df.sort_values(['motif_alt_id','start']))
  
        #display sorted motifs
        #display (new_df.sort_values(['motif_alt_id','start']))
        #display (new_df)
  
        new_df.sort_values(['motif_alt_id','start'],inplace = True)
  
        #start values(new_df) into array
        array_start = new_df['start'].values
        #stop values(new_df) into array
        array_stop = new_df['stop'].values
        #sequence name(new_df) values
        array_sequence_name = new_df['sequence_name'].values
        #array motifs
        array_motifs = new_df['motif_alt_id'].values
  
        #iterate through entire list/need treatment for combinations
        for idx, (seq_start,seq_stop,seq_name,motif) in enumerate( zip(array_start, array_stop, array_sequence_name, array_motifs) ):
            #current variables
            current_index = idx
            verify_index = False
            current_motif = motif
            current_seq_name = seq_name
            current_start = seq_start
            current_stop = seq_stop
      
            seq_type = seq_name.split('|')[0]
            condition_1 = seq_type == 'sequence_variation'
            condition_2 = seq_type == 'sequence_wild_type'
            condition_3 = seq_type == 'sequence_combinations'
      
            if condition_1 or condition_2:
                for idx, (seq_start,seq_stop,seq_name,motif) in enumerate( zip(array_start, array_stop, array_sequence_name, array_motifs) ):
                    snp_name = seq_name.split('|')[1]
                    #compare list with current values
                    condition_m = current_motif == motif
                    condition_sn = current_seq_name.split('|')[0] != seq_name.split('|')[0]
                    condition_sn2 = current_seq_name.split('|')[1] == seq_name.split('|')[1]
                    condition_csta = current_start == seq_start
                    condition_csta1 = seq_start == range(int(current_start-50),int(current_start))
                    condition_csto = current_stop == seq_stop
                    condition_csto1 = seq_stop == range(int(current_stop),int(current_stop + 50))
        
                if condition_m and condition_sn and condition_sn2 and (condition_csta or condition_csta1).any and (condition_csto or condition_csto1).any:
                    #new data frame droping all occurences of specified TF-snp pair
                    new_df = new_df.drop(new_df[(new_df.motif_alt_id == motif) & (new_df.sequence_name.str.contains(snp_name))].index)
            
            #missing tests/half completed/reset variables after processing    
            if condition_3:
                #snp list
                snp_name_list = [] 
                #variable to count the amount of snps
                count = 0
                #get number of all snps in combinations sequence
                for i in current_seq_name.split('|'):
                    if i[:2] == 'rs':
                        snp_name_list.append(i)
                        #print ("SNP NAME IN THE CURRENT SEQUENCE: ",i)
                        count = count + 1
                #print("COUNT",count)
                #transfer count to total ?
                total = count
                #print (total)
                multi_snps = 1
                count_al = 0
                #iterate through dictionary of snp and alleles
                for snp_name,snp_dict in zip(snp_name_list,dictionary_snp_allele):
                    #if snp in list of combinations is the same in the dictionary 
                    if snp_dict == snp_name:
                        #loop through alleles and count it 
                        for allele in dictionary_snp_allele[snp_dict]:
                            #print ("SNP DICT: ",snp_dict," ",allele)
                            count_al +=1
                        #multiply the amount of alleles to find the amount of combinations
                        multi_snps *= count_al
                        count_al = 0
        
                print ("THE AMOUNT OF COMBINATIONS: ",multi_snps)
                #position_in_list = pos_in_TF(current_start,current_stop)
                #snp_in_list = snp_in_TF(current_start,current_stop)
                count_comb = 0
        
                for idx, (seq_start,seq_stop,seq_name,motif) in enumerate( zip(array_start, array_stop, array_sequence_name, array_motifs) ):
          
                    #same motif
                    condition_m = current_motif == motif
                    #same snps
                    condition_sn2 = snp_name_list == seq_name.split('|')[1:total+1]
                    #range search needs maintain in snp start stop in snp range.
                    condition_csta = current_start == seq_start
                    condition_csto = current_stop == seq_stop
          
                    #print("MOTIF CURRENT",current_motif)
                    #print("MOTIF",motif)
                    #print("seq_name_list: ",snp_name_list)
                    #print("seq_name: ",seq_name.split('|')[1:total+1])
                    #print("condition_m: ",condition_m)
                    #print("condition_sn2: ",condition_sn2)
          
          
                    if condition_m and condition_sn2 and condition_csta and condition_csto:
                        print(motif)
                        print(count_comb)
                        count_comb+=1
                        #new data frame droping all occurences of specified TF-snps pair
                        #new_df = new_df.drop(new_df[(new_df.motif_alt_id == motif) & (new_df.sequence_name.str.contains(snp_list))].index)
                #print("COMBINATIONS",count_comb)
                if multi_snps <= count_comb:
                    print("COMBINATIONS",count_comb)
                    print("MOTIF",current_motif)
                    #only the first snp
                    new_df = new_df.drop(new_df[ (new_df.motif_alt_id == current_motif) & (new_df.sequence_name.str.contains(snp_name_list[0]))].index)
          
        #display (new_df.sort_values(['motif_alt_id','start']))
        print ("Tamanho DataFrame: ",len(new_df))
        return new_df
    
    if request.method == 'POST':
        #insert in a variable fimo path for fimo results
        fimo_res = request.form['fimo_tsv']
        df_fimo_output = pd.read_csv(fimo_res, sep='\t') 
        print(df_fimo_output)

        f_dataframe = filter_dataframe (fimo_res,True)

        f_dataframe = f_dataframe.drop_duplicates(keep="last")

        dataframe_out = f_dataframe.to_dict(orient='records')
        print("DATA FRAME DICT")
        print(dataframe_out)

    return jsonify(dataframe_out)
=======
    return jsonify(dict_snps)
>>>>>>> 11f9506f8fd474fb447f3b18b6eb81739dcc6a49
