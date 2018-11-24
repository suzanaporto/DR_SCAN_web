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
        snp = snp_input_form[2:]
        print(snp)
        a = get_snp_info(snp)
        print(a[0][0])
        return jsonify(a[0][0],snp_input_form)
        
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
    return jsonify(dict_snps)