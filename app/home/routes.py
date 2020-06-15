from app.home import blueprint
from flask import render_template, url_for
from flask_login import login_required
from flask import request
# imports
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
import time
import json
import app.home.sequence_gen as sg
from app.base.models import User
from app.base.models import Workflow
from app import db
from datetime import datetime, timedelta
import smtplib
from flask_cors import CORS, cross_origin
import threading
from flask import copy_current_request_context
import allel
from celery import Celery
from app import celery

#step_1---------function
def get_snp_info(snp_id):

  ### Download Snp Info based on its id/OK
  def request_info_by_id(snp_id):

    server = "https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/"
    try:
        r = requests.get(server+snp_id, headers={ "Content-Type" : "application/json"})
    except OSError:
        print("Connection Error")
        time.sleep(1)
    try:
        r
    except NameError:
        print('var not defined')
        r = requests.get(server+snp_id, headers={ "Content-Type" : "application/json"})
    
    if r in locals() and not r.ok:
      print("r: " + str(r.status_code))
      if r.status_code == 404:
          r.raise_for_status()
          sys.exit()
      for i in range(10):
          r = requests.get(server+snp_id, headers={ "Content-Type" : "application/json"})
          if r.status_code == 200:
              break
          time.sleep(5)
      if not r.ok:
        print("STATUS CODE: " + str(r.status_code))
        r.raise_for_status()
        sys.exit()

    decoded = r.json()

    r.connection.close()

    return decoded

  ### Parse information about the SNP requested/ not OK
  def parse_json(res_json):
    
    ##placements with allele idx
    num_seq_id = len(res_json['primary_snapshot_data']['placements_with_allele'])
    ##ban indexes of some placements with allele idx/verifying seq id
    true_hgvs = np.array([res_json['primary_snapshot_data']['placements_with_allele'][idx]['seq_id'] for idx in range(num_seq_id)])
    ban = np.array([x[:2] == 'NC' for x in true_hgvs])
    ###OLD VERSION
    #get num of genomic placements (versions)
    num_genomic_placements = len(res_json['primary_snapshot_data']['placements_with_allele']) 

    # sometimes there is empty fields in the genomic placements... Detect and remove them! 
    temp_assertion = np.array([res_json['primary_snapshot_data']['placements_with_allele'][idx]['placement_annot']['seq_id_traits_by_assembly'] for idx in range(num_genomic_placements)])
    ban_idexs = np.array([len(x)!=0 for x in temp_assertion])
    temp_assertion = temp_assertion[ban]

    # Get the genomic placements (versions) names
    gnenome_versions = [ x[0]["assembly_name"] for x in temp_assertion]

    # get the number of variations
    allele_variations = len([res_json['primary_snapshot_data']['placements_with_allele'][idx]['alleles'] for idx in range(1)][0])

    # get the snps idexes
    snp_idxs = np.array([[res_json['primary_snapshot_data']['placements_with_allele'][x]['alleles'][y]['hgvs'] for y in range(allele_variations)][1:] for x in range(num_genomic_placements)])

    # remove those without valid genomic placements (versions) names 
    snp_idxs = snp_idxs[ban].tolist()


    return {"gnenome_versions": gnenome_versions, "snp_idxs_list": snp_idxs}

  ### Parse snp information
  def snp_info_input(snp_id):
    base = re.compile("[^(\d)\w+]").split(snp_id)[3]
    base_chrom = re.compile("[^(\d)\w+]").split(snp_id)[0]
    dict = {
        "chrom" : re.compile(".0{2,}").split(base_chrom)[1],
        "location" : re.compile("[^(\d)]").split(base)[0],
        "allele_wt" : re.compile("[(\d)]").split(base)[-1],
        "allele_v": re.compile("[^(\d)\w+]").split(snp_id)[4]
    }

    return dict

  ### Download Snp Info based on its id
  
  res_json = request_info_by_id(snp_id)
  
  # sleep for 1 second
  time.sleep(3)
  # TODO When request is 404

  ##verify if snp has another name
  if "merged_snapshot_data" in res_json:
    new_snp = res_json['merged_snapshot_data']['merged_into'][0]
    res_json = request_info_by_id(new_snp)
    
  empty_df = pd.DataFrame()  
  ##verify if snp is actually an indel
  if res_json['primary_snapshot_data']['variant_type'] != 'snv':
    return empty_df
  
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

## step_2---------function
# state model verification
def apply_state_model(tissue, snp_list, snp_id, chrom):

        #Initializing dictionaries for state model files
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

        #Declaring paths
        STATE_15_PATH = "./app/home/state_models/ROADMAP_15_statemodel/"
        STATE_18_PATH = "./app/home/state_models/ROADMAP_18_statemodel/"
        STATE_25_PATH = "./app/home/state_models/ROADMAP_25_statemodel/"

        #Getting bed files from path
        state15_bed = STATE_15_PATH + tissue + '_15_coreMarks_stateno.bed.gz'
        state18_bed = STATE_18_PATH + tissue + '_18_core_K27ac_stateno.bed.gz'
        state25_bed = STATE_25_PATH + tissue + '_25_imputed12marks_stateno.bed.gz'

        #Using pybedtools to open .bed.gz files
        if ( os.path.exists(state15_bed) ):
            state15_bed_f = pybedtools.BedTool(state15_bed)
        if ( os.path.exists(state18_bed) ):
            state18_bed_f = pybedtools.BedTool(state18_bed)
        if ( os.path.exists(state25_bed) ):
            state25_bed_f = pybedtools.BedTool(state25_bed)

        # loop through state models and get regulatory element
        if 'state15_bed_f' in locals():
            for element in state15_bed_f:
                eq_chrom = element[0] == 'chr'+str(chrom)
                start_loc = int(element[1]) <= int(snp_list[1])
                stop_loc = int(element[2]) >= int(snp_list[1])
                if( (eq_chrom)  and (  (start_loc) and (stop_loc) ) ):
                    print(" | STATE MODEL 15 regulatory element of tissue "+ tissue + " :"  + state_15[int(element[3])])
                    row15 = {'snp_name':snp_id,'state_model':str(15),'tissue':tissue,'reg_elemnt':state_15[int(element[3])]}
                    break
        
        # loop through state models and get regulatory element
        if 'state18_bed_f' in locals():
            for element in state18_bed_f:
                eq_chrom = element[0] == 'chr'+str(chrom)
                start_loc = int(element[1]) <= int(snp_list[1])
                stop_loc = int(element[2]) >= int(snp_list[1])
                if( (eq_chrom)  and (  (start_loc) and (stop_loc) ) ):
                    print(" | STATE MODEL 18 regulatory element of tissue "+ tissue + ":"  + state_18[int(element[3])])
                    row18 = {'snp_name':snp_id,'state_model':str(18),'tissue':tissue,'reg_elemnt':state_18[int(element[3])]}
                    break
        
        # loop through state models and get regulatory element
        if 'state25_bed_f' in locals():
            for element in state25_bed_f:
                eq_chrom = element[0] == 'chr'+str(chrom)
                start_loc = int(element[1]) <= int(snp_list[1])
                stop_loc = int(element[2]) >= int(snp_list[1])
                if( (eq_chrom)  and (  (start_loc) and (stop_loc) ) ):
                    print(" | STATE MODEL 25 regulatory element of tissue "+ tissue + " :"  + state_25[int(element[3])])
                    row25 = {'snp_name':snp_id,'state_model':str(25),'tissue':tissue,'reg_elemnt':state_25[int(element[3])]}
                    break
        #Making sure variables do not exist so it has something to return
        if 'row25' not in locals():
            row25 = {'snp_name':"DO NOT EXIST",'state_model':str(25),'tissue':tissue,'reg_elemnt':"DO NOT EXIST"}
        if 'row15' not in locals():
            row15 = {'snp_name':"DO NOT EXIST",'state_model':str(15),'tissue':tissue,'reg_elemnt':"DO NOT EXIST"}
        if 'row18' not in locals():
            row18 = {'snp_name':"DO NOT EXIST",'state_model':str(18),'tissue':tissue,'reg_elemnt':"DO NOT EXIST"}
        
        #create row with all 3 rows
        row = {'snp_name':snp_id,
                'state_model':str(25)+" | "+str(18)+" | "+str(15),
                'tissue':tissue,
                'reg_elemnt':row25['reg_elemnt']+" | "+row18['reg_elemnt']+" | "+row15['reg_elemnt']
        }
        #return info to javascript
        return row

# step 5 function
#enhancer promoter interaction
def epi_search(snp_location,snp_chrom,snp_name,tissue,tissue_name,tissue_id):
    result = []
    file_type = []
    for i in range(2):
        if i == 0 and tissue_id[:1] == 'E':
            csv_name = './app/home/JEME/EncodeRoadMap/encoderoadmap_elasticnet.'+tissue+'.csv'
        elif i == 1 and tissue_id[:1] == 'E':
            csv_name = './app/home/JEME/EncodeRoadMap/encoderoadmap_lasso.'+tissue+'.csv'
        elif i == 0 and tissue_id[:1] == 'C':
            csv_name = './app/home/JEME/Fantom5/fantom5_elasticnet/fantom5_elasticnet.'+tissue+'.csv'
        elif i == 1 and tissue_id[:1] == 'C':
            csv_name = './app/home/JEME/Fantom5/fantom5_lasso/fantom5_lasso.'+tissue+'.csv'
        df_csv = pd.read_csv(csv_name,names=['Location','Gene','Score'],usecols=['Location','Gene','Score'])
        locations = df_csv['Location']
        idx_list = []
        for idx,row in enumerate(locations):
            ##creating variables
            new_row_list = re.compile("[\W]").split(row)
            list_chrom = new_row_list[0]
            list_enh_start = new_row_list[1]
            list_enh_end = new_row_list[2]
            ##creating conditions
            same_chrom = list_chrom == 'chr'+str(snp_chrom)
            bigger_start = int(list_enh_start) <= int(snp_location) 
            smaller_end = int(list_enh_end) >= int(snp_location)
            if same_chrom and (bigger_start and smaller_end):
                idx_list.append(idx)
                if i == 0:
                    file_type.append('Elasticnet')
                elif i == 1:
                    file_type.append('Lasso')
    
        rows = df_csv.loc[idx_list]
        df2 = pd.DataFrame(rows,columns=['Location','Gene','Score','Tissue']
                    ).reset_index().drop(columns=['index'])
        df2['Tissue'] = tissue_name
        df2['SNP'] = snp_name
        df2.insert(loc=5, column='File_Type', value=file_type)
        df2 = df2[['SNP', 'Location', 'Gene', 'Score', 'Tissue','File_Type']]
        result.append(df2)
#       display(df2)
        file_type = []
    df_final = pd.DataFrame(columns= ['SNP','Location','Gene','Score','Tissue','File_Type'])
    df_final = pd.concat(result)
    df_final = df_final.reset_index().drop(columns=['index'])
    return df_final

def epi_function(snps,tissues):
    ##dictionary roadmap epigenomics
    roadmap = {
    
        1:['E001','ES-I3 Cells'],
        2:['E002','ES-WA7 Cells'],
        3:['E003','H1 Cells'],
        4:['E004','H1 BMP4 Derived Mesendoderm Cultured Cells'],
        5:['E005','H1 BMP4 Derived Trophoblast Cultured Cells'],
        6:['E006','H1 Derived Mesenchymal Stem Cells'],
        7:['E007','H1 Derived Neuronal Progenitor Cultured Cells'],
        8:['E008','H9 Cells'],
        9:['E009','H9 Derived Neuronal Progenitor Cultured Cells'],
        10:['E010','H9 Derived Neuron Cultured Cells'],
        11:['E011','hESC Derived CD184+ Endoderm Cultured Cells'],
        12:['E012','hESC Derived CD56+ Ectoderm Cultured Cells'],
        13:['E013','hESC Derived CD56+ Mesoderm Cultured Cells'],
        14:['E014','HUES48 Cells'],
        15:['E015','HUES6 Cells'],
        16:['E016','HUES64 Cells'],
        17:['E017','IMR90 fetal lung fibroblasts Cell Line'],
        18:['E018','iPS-15b Cells'],
        19:['E019','iPS-18 Cells'],
        20:['E020','iPS-20b Cells'],
        21:['E021','iPS DF 6.9 Cells'],
        22:['E022','iPS DF 19.11 Cells'],
        23:['E023','Mesenchymal Stem Cell Derived Adipocyte Cultured Cells'],
        24:['E024','ES-UCSF4 Cells'],
        25:['E025','Adipose Derived Mesenchymal Stem Cell Cultured Cells'],
        26:['E026','Bone Marrow Derived Cultured Mesenchymal Stem Cells'],
        27:['E027','Breast Myoepithelial Primary Cells'],
        28:['E028','Breast variant Human Mammary Epithelial Cells (vHMEC)'],
        29:['E029','Primary monocytes from peripheral blood'],
        30:['E030','Primary neutrophils from peripheral blood'],
        31:['E031','Primary B cells from cord blood'],
        32:['E032','Primary B cells from peripheral blood'],
        33:['E033','Primary T cells from cord blood'],
        34:['E034','Primary T cells from peripheral blood'],
        35:['E035','Primary hematopoietic stem cells'],
        36:['E036','Primary hematopoietic stem cells short term culture'],
        37:['E037','Primary T helper memory cells from peripheral blood 2'],
        38:['E038','Primary T helper naive cells from peripheral blood'],
        39:['E039','Primary T helper naive cells from peripheral blood'],
        40:['E040','Primary T helper memory cells from peripheral blood 1'],
        41:['E041','Primary T helper cells PMA-I stimulated'],
        42:['E042','Primary T helper 17 cells PMA-I stimulated'],
        43:['E043','Primary T helper cells from peripheral blood'],
        44:['E044','Primary T regulatory cells from peripheral blood'],
        45:['E045','Primary T cells effector/memory enriched from peripheral blood'],
        46:['E046','Primary Natural Killer cells from peripheral blood'],
        47:['E047','Primary T CD8+ naive cells from peripheral blood'],
        48:['E048','Primary T CD8+ memory cells from peripheral blood'],
        49:['E049','Mesenchymal Stem Cell Derived Chondrocyte Cultured Cells'],
        50:['E050','Primary hematopoietic stem cells G-CSF-mobilized Female'],
        51:['E051','Primary hematopoietic stem cells G-CSF-mobilized Male'],
        52:['E052','Muscle Satellite Cultured Cells'],
        53:['E053','Cortex derived primary cultured neurospheres'],
        54:['E054','Ganglion Eminence derived primary cultured neurospheres'],
        55:['E055','Foreskin Fibroblast Primary Cells skin01'],
        56:['E056','Foreskin Fibroblast Primary Cells skin02'],
        57:['E057','Foreskin Keratinocyte Primary Cells skin02'],
        58:['E058','Foreskin Keratinocyte Primary Cells skin03'],
        59:['E059','Foreskin Melanocyte Primary Cells skin01'],
        60:['E061','Foreskin Melanocyte Primary Cells skin03'],
        61:['E062','Primary mononuclear cells from peripheral blood'],
        62:['E063','Adipose Nuclei'],
        63:['E065','Aorta'],
        64:['E066','Liver'],
        65:['E067','Brain Angular Gyrus'],
        66:['E068','Brain Anterior Caudate'],
        67:['E069','Brain Cingulate Gyrus'],
        68:['E070','Brain Germinal Matrix'],
        69:['E071','Brain Hippocampus Middle'],
        70:['E072','Brain Inferior Temporal Lobe'],
        71:['E073','Brain_Dorsolateral_Prefrontal_Cortex'],
        72:['E074','Brain Substantia Nigra'],
        73:['E075','Colonic Mucosa'],
        74:['E076','Colon Smooth Muscle'],
        75:['E077','Duodenum Mucosa'],
        76:['E078','Duodenum Smooth Muscle'],
        77:['E079','Esophagus'],
        78:['E080','Fetal Adrenal Gland'],
        79:['E081','Fetal Brain Male'],
        80:['E082','Fetal Brain Female'],
        81:['E083','Fetal Heart'],
        82:['E084','Fetal Intestine Large'],
        83:['E085','Fetal Intestine Small'],
        84:['E086','Fetal Kidney'],
        85:['E087','Pancreatic Islets'],
        86:['E088','Fetal Lung'],
        87:['E089','Fetal Muscle Trunk'],
        88:['E090','Fetal Muscle Leg'],
        89:['E091','Placenta'],
        90:['E092','Fetal Stomach'],
        91:['E093','Fetal Thymus'],
        92:['E094','Gastric'],
        93:['E095','Left Ventricle'],
        94:['E096','Lung'],
        95:['E097','Ovary'],
        96:['E098','Pancreas'],
        97:['E099','Placenta Amnion'],
        98:['E100','Psoas Muscle'],
        99:['E101','Rectal Mucosa Donor 29'],
        100:['E102','Rectal Mucosa Donor 31'],
        101:['E103','Rectal Smooth Muscle'],
        102:['E104','Right Atrium'],
        103:['E105','Right Ventricle'],
        104:['E106','Sigmoid Colon'],
        105:['E107','Skeletal Muscle Male'],
        106:['E108','Skeletal Muscle Female'],
        107:['E109','Small Intestine'],
        108:['E110','Stomach Mucosa'],
        109:['E111','Stomach Smooth Muscle'],
        110:['E112','Thymus'],
        111:['E113','Spleen'],
        112:['E114','A549 EtOH 0.02pct Lung Carcinoma Cell Line'],
        113:['E115','Dnd41 TCell Leukemia Cell Line'],
        114:['E116','GM12878 Lymphoblastoid Cells'],
        115:['E117','HeLa-S3 Cervical Carcinoma Cell Line'],
        116:['E118','HepG2 Hepatocellular Carcinoma Cell Line'],
        117:['E119','HMEC Mammary Epithelial Primary Cells'],
        118:['E120','HSMM Skeletal Muscle Myoblasts Cells'],
        119:['E121','HSMM cell derived Skeletal Muscle Myotubes Cells'],
        120:['E122','HUVEC Umbilical Vein Endothelial Primary Cells'],
        121:['E123','K562 Leukemia Cells'],
        122:['E124','Monocytes-CD14+ RO01746 Primary Cells'],
        123:['E125','NH-A Astrocytes Primary Cells'],
        124:['E126','NHDF-Ad Adult Dermal Fibroblast Primary Cells'],
        125:['E127','NHEK-Epidermal Keratinocyte Primary Cells'],
        126:['E128','NHLF Lung Fibroblast Primary Cells'],
        127:['E129','Osteoblast Primary Cells'],
    }
    ##dictionary fantom 5
    fantom_5 = {
    
        1:['CNhs10608','Clontech Human Universal Reference Total RNA pool1'],
        2:['CNhs10610','SABiosciences XpressRef Human Universal Total RNA pool1'],
        3:['CNhs10612','Universal RNA - Human Normal Tissues Biochain pool1'],
        4:['CNhs10615','adipose tissue adult pool1'],
        5:['CNhs10616','bladder adult pool1'],
        6:['CNhs10617','brain adult pool1'],
        7:['CNhs10618','cervix adult pool1'],
        8:['CNhs10619','colon adult pool1'],
        9:['CNhs10620','esophagus adult pool1'],
        10:['CNhs10621','heart adult pool1'],
        11:['CNhs10622','kidney adult pool1'],
        12:['CNhs10624','liver adult pool1'],
        13:['CNhs10625','lung adult pool1'],
        14:['CNhs10626','ovary adult pool1'],
        15:['CNhs10627','placenta adult pool1'],
        16:['CNhs10628','prostate adult pool1'],
        17:['CNhs10629','skeletal muscle adult pool1'],
        18:['CNhs10630','small intestine adult pool1'],
        19:['CNhs10631','spleen adult pool1'],
        20:['CNhs10632','testis adult pool1'],
        21:['CNhs10633','thymus adult pool1'],
        22:['CNhs10634','thyroid adult pool1'],
        23:['CNhs10635','trachea adult pool1'],
        24:['CNhs10636','retina adult pool1'],
        25:['CNhs10637','temporal lobe adult pool1'],
        26:['CNhs10638','postcentral gyrus adult pool1'],
        27:['CNhs10640','pons adult pool1'],
        28:['CNhs10641','parietal lobe adult pool1'],
        29:['CNhs10642','paracentral gyrus adult pool1'],
        30:['CNhs10643','occipital pole adult pool1'],
        31:['CNhs10644','nucleus accumbens adult pool1'],
        32:['CNhs10645','medulla oblongata adult pool1'],
        33:['CNhs10646','insula adult pool1'],
        34:['CNhs10647','frontal lobe adult pool1'],
        35:['CNhs10648','dura mater adult donor1'],
        36:['CNhs10649','corpus callosum adult pool1'],
        37:['CNhs10650','thymus fetal pool1'],
        38:['CNhs10651','spleen fetal pool1'],
        39:['CNhs10652','kidney fetal pool1'],
        40:['CNhs10653','heart fetal pool1'],
        41:['CNhs10654','tonsil adult pool1'],
        42:['CNhs10722','acute myeloid leukemia (FAB M5) cell line:THP-1 (fresh)'],
        43:['CNhs10723','acute myeloid leukemia (FAB M5) cell line:THP-1 (revived)'],
        44:['CNhs10724','acute myeloid leukemia (FAB M5) cell line:THP-1 (thawed)'],
        45:['CNhs10726','lung adenocarcinoma cell line:PC-14'],
        46:['CNhs10727','chronic myelogenous leukemia cell line:KU812'],
        47:['CNhs10728','extraskeletal myxoid chondrosarcoma cell line:H-EMC-SS'],
        48:['CNhs10729','renal cell carcinoma cell line:OS-RC-2'],
        49:['CNhs10730','malignant trichilemmal cyst cell line:DJM-1'],
        50:['CNhs10731','glioma cell line:GI-1'],
        51:['CNhs10732','maxillary sinus tumor cell line:HSQ-89'],
        52:['CNhs10733','gall bladder carcinoma cell line:TGBC2TKB'],
        53:['CNhs10734','papillotubular adenocarcinoma cell line:TGBC18TKB'],
        54:['CNhs10735','transitional-cell carcinoma cell line:5637'],
        55:['CNhs10736','breast carcinoma cell line:MDA-MB-453'],
        56:['CNhs10737','colon carcinoma cell line:COLO-320'],
        57:['CNhs10738','adult T-cell leukemia cell line:ATN-1'],
        58:['CNhs10739','Burkitt\'s lymphoma cell line:DAUDI'],
        59:['CNhs10740','choriocarcinoma cell line:BeWo'],
        60:['CNhs10741','splenic lymphoma with villous lymphocytes cell line:SLVL'],
        61:['CNhs10742','astrocytoma cell line:TM-31'],
        62:['CNhs10743','epidermoid carcinoma cell line:A431'],
        63:['CNhs10744','b cell line:RPMI1788'],
        64:['CNhs10745','anaplastic carcinoma cell line:8305C'],
        65:['CNhs10746','acute lymphoblastic leukemia (T-ALL) cell line:HPB-ALL'],
        66:['CNhs10747','non T non B acute lymphoblastic leukemia (ALL) cell line:P30/OHK'],
        67:['CNhs10748','epidermoid carcinoma cell line:Ca Ski'],
        68:['CNhs10750','bile duct carcinoma cell line:HuCCT1'],
        69:['CNhs10751','giant cell carcinoma cell line:Lu99B'],
        70:['CNhs10752','oral squamous cell carcinoma cell line:Ca9-22'],
        71:['CNhs10753','signet ring carcinoma cell line:Kato III'],
        72:['CNhs10837','Endothelial Cells - Aortic donor0'],
        73:['CNhs10838','Smooth Muscle Cells - Aortic donor0'],
        74:['CNhs10839','Smooth Muscle Cells - Umbilical artery donor0'],
        75:['CNhs10842','Retinal Pigment Epithelial Cells donor0'],
        76:['CNhs10843','Urothelial cells donor0'],
        77:['CNhs10844','Mesenchymal stem cells - adipose donor0'],
        78:['CNhs10845','Mesenchymal stem cells - hepatic donor0'],
        79:['CNhs10847','Sebocyte donor1'],
        80:['CNhs10848','Fibroblast - Gingival donor4 (GFH2)'],
        81:['CNhs10850','Mesothelial Cells donor1'],
        82:['CNhs10851','Sertoli Cells donor1'],
        83:['CNhs10852','CD14+ Monocytes donor1'],
        84:['CNhs10853','CD4+ T Cells donor1'],
        85:['CNhs10854','CD8+ T Cells donor1'],
        86:['CNhs10855','Dendritic Cells - monocyte immature derived donor1 tech_rep1'],
        87:['CNhs10857','Dendritic Cells - plasmacytoid donor1'],
        88:['CNhs10858','CD14+ monocyte derived endothelial progenitor cells donor1'],
        89:['CNhs10859','Natural Killer Cells donor1'],
        90:['CNhs10860','Peripheral Blood Mononuclear Cells donor1'],
        91:['CNhs10861','Macrophage - monocyte derived donor1'],
        92:['CNhs10862','Neutrophils donor1'],
        93:['CNhs10863','Smooth Muscle Cells - Brain Vascular donor1'],
        94:['CNhs10864','Astrocyte - cerebral cortex donor1'],
        95:['CNhs10865','Endothelial Cells - Lymphatic donor1'],
        96:['CNhs10866','Fibroblast - Gingival donor1'],
        97:['CNhs10867','Fibroblast - Periodontal Ligament donor1'],
        98:['CNhs10868','Smooth Muscle Cells - Colonic donor1'],
        99:['CNhs10869','Skeletal Muscle Satellite Cells donor1'],
        100:['CNhs10870','Myoblast donor1'],
        101:['CNhs10871','Ciliary Epithelial Cells donor1'],
        102:['CNhs10872','Endothelial Cells - Umbilical vein donor1'],
        103:['CNhs10874','Fibroblast - Aortic Adventitial donor1'],
        104:['CNhs10875','Intestinal epithelial cells (polarized) donor1'],
        105:['CNhs10876','Anulus Pulposus Cell donor1'],
        106:['CNhs10877','Pancreatic stromal cells donor1'],
        107:['CNhs10878','Fibroblast - Pulmonary Artery donor1'],
        108:['CNhs10879','Keratinocyte - oral donor1'],
        109:['CNhs10881	','Nucleus Pulposus Cell donor1'],
        110:['CNhs10882','Prostate Epithelial Cells (polarized) donor1'],
        111:['CNhs10883','Prostate Stromal Cells donor1'],
        112:['CNhs10884','Small Airway Epithelial Cells donor1'],
        113:['CNhs11045','cord blood derived cell line:COBL-a untreated'],
        114:['CNhs11046','embryonic kidney cell line: HEK293/SLAM untreated'],
        115:['CNhs11047','embryonic kidney cell line: HEK293/SLAM infection 24hr'],
        116:['CNhs11049','cord blood derived cell line:COBL-a 24h infection(-C)'],
        117:['CNhs11050','cord blood derived cell line:COBL-a 24h infection'],
        118:['CNhs11051','Adipocyte - breast donor1'],
        119:['CNhs11052','Preadipocyte - breast donor1'],
        120:['CNhs11054','Adipocyte - omental donor1'],
        121:['CNhs11057','Mesenchymal Stem Cells - Wharton\'s Jelly donor1'],
        122:['CNhs11061','Gingival epithelial cells donor1 (GEA11)'],
        123:['CNhs11062','Dendritic Cells - monocyte immature derived donor1 tech_rep2'],
        124:['CNhs11063','Neural stem cells donor1'],
        125:['CNhs11064','Keratinocyte - epidermal donor1'],
        126:['CNhs11065','Preadipocyte - omental donor1'],
        127:['CNhs11073','Mast cell - stimulated donor1'],
        128:['CNhs11074','Fibroblast - skin spinal muscular atrophy donor1'],
        129:['CNhs11075','Whole blood (ribopure) donor090325 donation1'],
        130:['CNhs11076','Whole blood (ribopure) donor090325 donation2'],
        131:['CNhs11077','Mammary Epithelial Cell donor1'],
        132:['CNhs11079','Placental Epithelial Cells donor1'],
        133:['CNhs11082','Preadipocyte - visceral donor1'],
        134:['CNhs11083','Skeletal Muscle Cells donor1'],
        135:['CNhs11084','Skeletal muscle cells differentiated into Myotubes - multinucleated donor1'],
        136:['CNhs11085','Smooth Muscle Cells - Aortic donor1'],
        137:['CNhs11086','Smooth Muscle Cells - Brachiocephalic donor1'],
        138:['CNhs11087','Smooth Muscle Cells - Carotid donor1'],
        139:['CNhs11088','Smooth Muscle Cells - Coronary Artery donor1'],
        140:['CNhs11090','Smooth Muscle Cells - Subclavian Artery donor1'],
        141:['CNhs11091','Smooth Muscle Cells - Umbilical Artery donor1'],
        142:['CNhs11092','Tracheal Epithelial Cells donor1'],
        143:['CNhs11100','ductal cell carcinoma cell line:KLM-1'],
        144:['CNhs11183','schwannoma cell line:HS-PSS'],
        145:['CNhs11185','glioblastoma cell line:A172'],
        146:['CNhs11243','prostate cancer cell line:PC-3'],
        147:['CNhs11244','synovial sarcoma cell line:HS-SY-II'],
        148:['CNhs11245','schwannoma cell line:HS-PSS tech_rep2'],
        149:['CNhs11247','epithelioid sarcoma cell line:HS-ES-1'],
        150:['CNhs11248','glioblastoma cell line:A172 tech_rep2'],
        151:['CNhs11249','endometrial stromal sarcoma cell line:OMC-9'],
        152:['CNhs11250','chronic myelogenous leukemia cell line:K562'],
        153:['CNhs11251','acute lymphoblastic leukemia (B-ALL) cell line:BALL-1'],
        154:['CNhs11252','squamous cell carcinoma cell line:EC-GI-10'],
        155:['CNhs11253','acute lymphoblastic leukemia (T-ALL) cell line:Jurkat'],
        156:['CNhs11254','melanoma cell line:G-361'],
        157:['CNhs11255','rectal cancer cell line:TT1TKB'],
        158:['CNhs11256','gall bladder carcinoma cell line:TGBC14TKB'],
        159:['CNhs11257','renal cell carcinoma cell line:TUHR10TKB'],
        160:['CNhs11258','myeloma cell line:PCM6'],
        161:['CNhs11259','ductal cell carcinoma cell line:MIA Paca2'],
        162:['CNhs11260','prostate cancer cell line:DU145'],
        163:['CNhs11261','transitional-cell carcinoma cell line:JMSU1'],
        164:['CNhs11263','mesothelioma cell line:ACC-MESO-1'],
        165:['CNhs11264','mesothelioma cell line:ACC-MESO-4'],
        166:['CNhs11265','bile duct carcinoma cell line:TFK-1'],
        167:['CNhs11266','endometrial carcinoma cell line:OMC-2'],
        168:['CNhs11267','retinoblastoma cell line:Y79'],
        169:['CNhs11268','Burkitt\'s lymphoma cell line:RAJI'],
        170:['CNhs11269','rhabdomyosarcoma cell line:RMS-YM'],
        171:['CNhs11270','signet ring carcinoma cell line:NUGC-4'],
        172:['CNhs11271','hepatoma cell line:Li-7'],
        173:['CNhs11272','glioblastoma cell line:T98G'],
        174:['CNhs11273','squamous cell lung carcinoma cell line:EBC-1'],
        175:['CNhs11274','giant cell carcinoma cell line:LU65'],
        176:['CNhs11275','lung adenocarcinoma cell line:A549'],
        177:['CNhs11276','neuroblastoma cell line:CHP-134'],
        178:['CNhs11277','large cell lung carcinoma cell line:IA-LM'],
        179:['CNhs11279','osteosarcoma cell line:143B/TK^(-)neo^(R)'],
        180:['CNhs11280','colon carcinoma cell line:CACO-2'],
        181:['CNhs11281','melanoma cell line:COLO 679'],
        182:['CNhs11282','acute lymphoblastic leukemia (B-ALL) cell line:NALM-6'],
        183:['CNhs11283','cholangiocellular carcinoma cell line:HuH-28'],
        184:['CNhs11284','neuroblastoma cell line:NB-1'],
        185:['CNhs11285','small cell lung carcinoma cell line:LK-2'],
        186:['CNhs11286','gastric cancer cell line:AZ521'],
        187:['CNhs11287','oral squamous cell carcinoma cell line:HO-1-u-1'],
        188:['CNhs11288','cervical cancer cell line:D98-AH2'],
        189:['CNhs11289','cervical cancer cell line:ME-180'],
        190:['CNhs11290','osteosarcoma cell line:HS-Os-1'],
        191:['CNhs11303','Melanocyte - light donor1'],
        192:['CNhs11305','Smooth Muscle Cells - Aortic donor2'],
        193:['CNhs11309','Smooth Muscle Cells - Aortic donor3'],
        194:['CNhs11311','Osteoblast - differentiated donor1'],
        195:['CNhs11317','Pericytes donor1'],
        196:['CNhs11319','Fibroblast - Choroid Plexus donor1'],
        197:['CNhs11320','Meningeal Cells donor1'],
        198:['CNhs11321','Astrocyte - cerebellum donor1'],
        199:['CNhs11322','Fibroblast - Lymphatic donor1'],
        200:['CNhs11323','Esophageal Epithelial Cells donor1'],
        201:['CNhs11324','Smooth Muscle Cells - Esophageal donor1'],
        202:['CNhs11325','Alveolar Epithelial Cells donor1'],
        203:['CNhs11328','Smooth Muscle Cells - Bronchial donor1'],
        204:['CNhs11329','Smooth Muscle Cells - Tracheal donor1'],
        205:['CNhs11330','Renal Proximal Tubular Epithelial Cell donor1'],
        206:['CNhs11331','Renal Cortical Epithelial Cells donor1'],
        207:['CNhs11332','Renal Epithelial Cells donor1'],
        208:['CNhs11333','Renal Mesangial Cells donor1'],
        209:['CNhs11334','Urothelial Cells donor1'],
        210:['CNhs11335','Hepatic Stellate Cells (lipocyte) donor1'],
        211:['CNhs11336','Corneal Epithelial Cells donor1'],
        212:['CNhs11337','Keratocytes donor1'],
        213:['CNhs11338','Retinal Pigment Epithelial Cells donor1'],
        214:['CNhs11339','Fibroblast - Conjunctival donor1'],
        215:['CNhs11340','Trabecular Meshwork Cells donor1'],
        216:['CNhs11341','Amniotic Epithelial Cells donor1'],
        217:['CNhs11344','Mesenchymal Stem Cells - bone marrow donor1'],
        218:['CNhs11345','Mesenchymal Stem Cells - adipose donor1'],
        219:['CNhs11347','Mesenchymal Stem Cells - umbilical donor1'],
        220:['CNhs11349','Mesenchymal Stem Cells - amniotic membrane donor1'],
        221:['CNhs11350','Multipotent Cord Blood Unrestricted Somatic Stem Cells donor1'],
        222:['CNhs11351','Fibroblast - skin normal donor1'],
        223:['CNhs11352','Fibroblast - skin walker warburg donor1'],
        224:['CNhs11353','Fibroblast - skin dystrophia myotonica donor1'],
        225:['CNhs11354','Fibroblast - skin dystrophia myotonica donor2'],
        226:['CNhs11371','Adipocyte - subcutaneous donor2'],
        227:['CNhs11372','Chondrocyte - de diff donor2'],
        228:['CNhs11373','Chondrocyte - re diff donor2'],
        229:['CNhs11375','Endothelial Cells - Aortic donor2'],
        230:['CNhs11376','Endothelial Cells - Microvascular donor2'],
        231:['CNhs11377','Endothelial Cells - Vein donor2'],
        232:['CNhs11378','Fibroblast - Cardiac donor2'],
        233:['CNhs11379','Fibroblast - Dermal donor2'],
        234:['CNhs11381','Keratinocyte - epidermal donor2'],
        235:['CNhs11382','Mammary Epithelial Cell donor2'],
        236:['CNhs11383','Melanocyte - light donor2'],
        237:['CNhs11384','Neural stem cells donor2'],
        238:['CNhs11385','Osteoblast donor2'],
        239:['CNhs11386','Placental Epithelial Cells donor2'],
        240:['CNhs11671','Whole blood (ribopure) donor090309 donation2'],
        241:['CNhs11672','Whole blood (ribopure) donor090612 donation1'],
        242:['CNhs11673','Whole blood (ribopure) donor090612 donation2'],
        243:['CNhs11675','Whole blood (ribopure) donor090309 donation1'],
        244:['CNhs11676','uterus adult pool1'],
        245:['CNhs11677','salivary gland adult pool1'],
        246:['CNhs11680','lung fetal donor1'],
        247:['CNhs11714','chronic lymphocytic leukemia (T-CLL) cell line:SKW-3'],
        248:['CNhs11715','Hodgkin\'s lymphoma cell line:HD-Mar2'],
        249:['CNhs11716','papillary adenocarcinoma cell line:8505C'],
        250:['CNhs11717','oral squamous cell carcinoma cell line:HSC-3'],
        251:['CNhs11718','mesenchymal stem cell line:Hu5/E18'],
        252:['CNhs11722','leiomyoma cell line:10964C'],
        253:['CNhs11723','leiomyoma cell line:15242A'],
        254:['CNhs11724','leiomyoma cell line:15425'],
        255:['CNhs11725','argyrophil small cell carcinoma cell line:TC-YIK'],
        256:['CNhs11726','testicular germ cell embryonal carcinoma cell line:NEC8'],
        257:['CNhs11728','Wilms\' tumor cell line:HFWT'],
        258:['CNhs11729','myxofibrosarcoma cell line:MFH-ino'],
        259:['CNhs11731','embryonic pancreas cell line:1B2C6'],
        260:['CNhs11732','embryonic pancreas cell line:1C3D3'],
        261:['CNhs11733','embryonic pancreas cell line:1C3IKEI'],
        262:['CNhs11734','small-cell gastrointestinal carcinoma cell line:ECC4'],
        263:['CNhs11736','small cell gastrointestinal carcinoma cell line:ECC10'],
        264:['CNhs11737','gastric adenocarcinoma cell line:MKN1'],
        265:['CNhs11738','gastrointestinal carcinoma cell line:ECC12'],
        266:['CNhs11739','squamous cell carcinoma cell line:T3M-5'],
        267:['CNhs11740','granulosa cell tumor cell line:KGN'],
        268:['CNhs11741','diffuse large B-cell lymphoma cell line:CTB-1'],
        269:['CNhs11742','hepatoblastoma cell line:HuH-6'],
        270:['CNhs11744','neuroectodermal tumor cell line:FU-RPNT-1'],
        271:['CNhs11745','clear cell carcinoma cell line:JHOC-5'],
        272:['CNhs11746','serous adenocarcinoma cell line:JHOS-2'],
        273:['CNhs11747','carcinosarcoma cell line:JHUCS-1'],
        274:['CNhs11748','endometrioid adenocarcinoma cell line:JHUEM-1'],
        275:['CNhs11750','lens epithelial cell line:SRA 01/04'],
        276:['CNhs11752','mucinous adenocarcinoma cell line:JHOM-1'],
        277:['CNhs11753','neuroectodermal tumor cell line:FU-RPNT-2'],
        278:['CNhs11755','smooth muscle adult pool1'],
        279:['CNhs11756','pancreas adult donor1'],
        280:['CNhs11757','heart adult diseased post-infarction donor1'],
        281:['CNhs11758','heart adult diseased donor1'],
        282:['CNhs11760','aorta adult pool1'],
        283:['CNhs11761','blood adult pool1'],
        284:['CNhs11762','eye fetal donor1'],
        285:['CNhs11763','uterus fetal donor1'],
        286:['CNhs11764','spinal cord fetal donor1'],
        287:['CNhs11765','umbilical cord fetal donor1'],
        288:['CNhs11766','trachea fetal donor1'],
        289:['CNhs11768','tongue fetal donor1'],
        290:['CNhs11769','thyroid fetal donor1'],
        291:['CNhs11770','throat fetal donor1'],
        292:['CNhs11771','stomach fetal donor1'],
        293:['CNhs11772','temporal lobe fetal donor1 tech_rep1'],
        294:['CNhs11773','small intestine fetal donor1'],
        295:['CNhs11774','skin fetal donor1'],
        296:['CNhs11776','skeletal muscle fetal donor1'],
        297:['CNhs11777','rectum fetal donor1'],
        298:['CNhs11779','diaphragm fetal donor1'],
        299:['CNhs11780','colon fetal donor1'],
        300:['CNhs11781','duodenum fetal donor1 tech_rep1'],
        301:['CNhs11782','parietal lobe fetal donor1'],
        302:['CNhs11784','occipital lobe fetal donor1'],
        303:['CNhs11786','lung right lower lobe adult donor1'],
        304:['CNhs11787','occipital lobe adult donor1'],
        305:['CNhs11788','lymph node adult donor1'],
        306:['CNhs11789','left ventricle adult donor1'],
        307:['CNhs11790','left atrium adult donor1'],
        308:['CNhs11792','breast adult donor1'],
        309:['CNhs11794','colon adult donor1'],
        310:['CNhs11795','cerebellum adult pool1'],
        311:['CNhs11796','brain adult donor1'],
        312:['CNhs11797','brain fetal pool1'],
        313:['CNhs11798','liver fetal pool1'],
        314:['CNhs11810','oral squamous cell carcinoma cell line:SAS'],
        315:['CNhs11811','neuroblastoma cell line:NH-12'],
        316:['CNhs11812','small cell lung carcinoma cell line:WA-hT'],
        317:['CNhs11813','xeroderma pigentosum b cell line:XPL 17'],
        318:['CNhs11814','embryonic pancreas cell line:2C6'],
        319:['CNhs11818','neuroblastoma cell line:NBsusSR'],
        320:['CNhs11819','gastric adenocarcinoma cell line:MKN45'],
        321:['CNhs11820','choriocarcinoma cell line:T3M-3'],
        322:['CNhs11821','myxofibrosarcoma cell line:NMFH-1'],
        323:['CNhs11824','glassy cell carcinoma cell line:HOKUG'],
        324:['CNhs11825','large cell non-keratinizing squamous carcinoma cell line:SKG-II-SF'],
        325:['CNhs11827','serous cystadenocarcinoma cell line:HTOA'],
        326:['CNhs11828','tridermal teratoma cell line:HGRT'],
        327:['CNhs11829','sacrococcigeal teratoma cell line:HTST'],
        328:['CNhs11830','peripheral neuroectodermal tumor cell line:KU-SN'],
        329:['CNhs11832','pancreatic carcinoma cell line:NOR-P1'],
        330:['CNhs11834','carcinoid cell line:NCI-H1770'],
        331:['CNhs11835','osteoclastoma cell line:Hs 706.T'],
        332:['CNhs11836','Ewing\'s sarcoma cell line:Hs 863.T'],
        333:['CNhs11838','alveolar cell carcinoma cell line:SW 1573'],
        334:['CNhs11840','bronchioalveolar carcinoma cell line:NCI-H358'],
        335:['CNhs11841','bronchogenic carcinoma cell line:ChaGo-K-1'],
        336:['CNhs11842','fibrous histiocytoma cell line:GCT TIB-223'],
        337:['CNhs11843','hairy cell leukemia cell line:Mo'],
        338:['CNhs11844','acantholytic squamous carcinoma cell line:HCC1806'],
        339:['CNhs11845','biphenotypic B myelomonocytic leukemia cell line:MV-4-11'],
        340:['CNhs11846','carcinoid cell line:SK-PN-DW'],
        341:['CNhs11848','leiomyoblastoma cell line:G-402'],
        342:['CNhs11849','pharyngeal carcinoma cell line:Detroit 562'],
        343:['CNhs11851','liposarcoma cell line:SW 872'],
        344:['CNhs11852','lymphangiectasia cell line:DS-1'],
        345:['CNhs11853','neuroepithelioma cell line:SK-N-MC'],
        346:['CNhs11854','neurofibroma cell line:Hs 53.T'],
        347:['CNhs11856','pagetoid sarcoma cell line:Hs 925.T'],
        348:['CNhs11857','spindle cell sarcoma cell line:Hs 132.T'],
        349:['CNhs11858','mycosis fungoides T cell lymphoma cell line:HuT 102 TIB-162'],
        350:['CNhs11859','leukemia chronic megakaryoblastic cell line:MEG-01'],
        351:['CNhs11860','fibrosarcoma cell line:HT-1080'],
        352:['CNhs11861','medulloblastoma cell line:ONS-76'],
        353:['CNhs11862','bronchial squamous cell carcinoma cell line:KNS-62'],
        354:['CNhs11864','acute myeloid leukemia (FAB M2) cell line:NKM-1'],
        355:['CNhs11865','chronic myelogenous leukemia (CML) cell line:MEG-A2'],
        356:['CNhs11866','neuroectodermal tumor cell line:TASK1'],
        357:['CNhs11867','NK T cell leukemia cell line:KHYG-1'],
        358:['CNhs11868','hepatic mesenchymal tumor cell line:LI90'],
        359:['CNhs11869','somatostatinoma cell line:QGP-1'],
        360:['CNhs11870','liposarcoma cell line:KMLS-1'],
        361:['CNhs11872','thyroid carcinoma cell line:TCO-1'],
        362:['CNhs11873','mucinous cystadenocarcinoma cell line:MCAS'],
        363:['CNhs11875','choriocarcinoma cell line:SCH'],
        364:['CNhs11876','testicular germ cell embryonal carcinoma cell line:ITO-II'],
        365:['CNhs11877','rhabdomyosarcoma cell line:KYM-1'],
        366:['CNhs11878','teratocarcinoma cell line:NCC-IT-A3'],
        367:['CNhs11880','keratoacanthoma cell line:HKA-1'],
        368:['CNhs11881','anaplastic large cell lymphoma cell line:Ki-JK'],
        369:['CNhs11882','adenocarcinoma cell line:IM95m'],
        370:['CNhs11883','tubular adenocarcinoma cell line:SUIT-2'],
        371:['CNhs11884','teratocarcinoma cell line:NCR-G1'],
        372:['CNhs11885','small cell cervical cancer cell line:HCSC-1'],
        373:['CNhs11886','chronic myeloblastic leukemia (CML) cell line:KCL-22'],
        374:['CNhs11888','acute myeloid leukemia (FAB M7) cell line:MKPL-1'],
        375:['CNhs11889','anaplastic squamous cell carcinoma cell line:RPMI 2650'],
        376:['CNhs11890','teratocarcinoma cell line:PA-1'],
        377:['CNhs11891','hereditary spherocytic anemia cell line:WIL2-NS'],
        378:['CNhs11892','Wilms\' tumor cell line:G-401'],
        379:['CNhs11893','adrenal cortex adenocarcinoma cell line:SW-13'],
        380:['CNhs11894','normal embryonic palatal mesenchymal cell line:HEPM'],
        381:['CNhs11896','Gingival epithelial cells donor2 (GEA14)'],
        382:['CNhs11897','CD14+ monocyte derived endothelial progenitor cells donor2'],
        383:['CNhs11899','Macrophage - monocyte derived donor2'],
        384:['CNhs11900','Smooth Muscle Cells - Brain Vascular donor2'],
        385:['CNhs11901','Endothelial Cells - Lymphatic donor2'],
        386:['CNhs11902','Preadipocyte - omental donor2'],
        387:['CNhs11903','Gingival epithelial cells donor3 (GEA15)'],
        388:['CNhs11904','CD14+ monocyte derived endothelial progenitor cells donor3'],
        389:['CNhs11905','Neutrophils donor3'],
        390:['CNhs11906','Endothelial Cells - Lymphatic donor3'],
        391:['CNhs11907','Fibroblast - Periodontal Ligament donor3'],
        392:['CNhs11908','Myoblast donor3'],
        393:['CNhs11909','Fibroblast - Cardiac donor4'],
        394:['CNhs11911','Fibroblast - skin spinal muscular atrophy donor2'],
        395:['CNhs11912','Fibroblast - skin spinal muscular atrophy donor3'],
        396:['CNhs11913','Fibroblast - skin dystrophia myotonica donor3'],
        397:['CNhs11914','Fibroblast - skin normal donor2'],
        398:['CNhs11920','Smooth Muscle Cells - Prostate donor1'],
        399:['CNhs11923','Chondrocyte - de diff donor1'],
        400:['CNhs11925','Endothelial Cells - Microvascular donor1'],
        401:['CNhs11926','Endothelial Cells - Thoracic donor1'],
        402:['CNhs11927','Smooth Muscle Cells - Uterine donor3'],
        403:['CNhs11930','clear cell carcinoma cell line:TEN'],
        404:['CNhs11931','bone marrow stromal cell line:StromaNKtert'],
        405:['CNhs11932','basal cell carcinoma cell line:TE 354.T'],
        406:['CNhs11933','pleomorphic hepatocellular carcinoma cell line:SNU-387'],
        407:['CNhs11934','myelodysplastic syndrome cell line:SKM-1'],
        408:['CNhs11935','lymphoma malignant hairy B-cell cell line:MLMA'],
        409:['CNhs11943','breast carcinoma cell line:MCF7'],
        410:['CNhs11944','mixed mullerian tumor cell line:HTMMT'],
        411:['CNhs11945','meningioma cell line:HKBMM'],
        412:['CNhs11948','Whole blood (ribopure) donor090309 donation3'],
        413:['CNhs11949','Whole blood (ribopure) donor090612 donation3'],
        414:['CNhs11950','normal intestinal epithelial cell line:FHs 74 Int'],
        415:['CNhs11951','Sebocyte donor2'],
        416:['CNhs11952','Fibroblast - Gingival donor5 (GFH3)'],
        417:['CNhs11953','Fibroblast - Periodontal Ligament donor5 (PL30)'],
        418:['CNhs11954','CD14+ Monocytes donor2'],
        419:['CNhs11955','CD4+ T Cells donor2'],
        420:['CNhs11956','CD8+ T Cells donor2'],
        421:['CNhs11957','Natural Killer Cells donor2'],
        422:['CNhs11958','Peripheral Blood Mononuclear Cells donor2'],
        423:['CNhs11959','Neutrophils donor2'],
        424:['CNhs11960','Astrocyte - cerebral cortex donor2'],
        425:['CNhs11961','Fibroblast - Gingival donor2'],
        426:['CNhs11962','Fibroblast - Periodontal Ligament donor2'],
        427:['CNhs11963','Smooth Muscle Cells - Colonic donor2'],
        428:['CNhs11964','Skeletal Muscle Satellite Cells donor2'],
        429:['CNhs11965','Myoblast donor2'],
        430:['CNhs11966','Ciliary Epithelial Cells donor2'],
        431:['CNhs11967','Endothelial Cells - Umbilical vein donor2'],
        432:['CNhs11969','Adipocyte - breast donor2'],
        433:['CNhs11971','Preadipocyte - breast donor2'],
        434:['CNhs11972','Prostate Epithelial Cells donor2'],
        435:['CNhs11973','Prostate Stromal Cells donor2'],
        436:['CNhs11975','Small Airway Epithelial Cells donor2'],
        437:['CNhs11976','Smooth Muscle Cells - Prostate donor2'],
        438:['CNhs11977','Endothelial Cells - Artery donor2'],
        439:['CNhs11978','Endothelial Cells - Thoracic donor2'],
        440:['CNhs11979','Hair Follicle Dermal Papilla Cells donor2'],
        441:['CNhs11980','Osteoblast - differentiated donor2'],
        442:['CNhs11981','Preadipocyte - subcutaneous donor2'],
        443:['CNhs11982','Preadipocyte - visceral donor2'],
        444:['CNhs11987','Smooth Muscle Cells - Coronary Artery donor2'],
        445:['CNhs11988','Smooth Muscle Cells - Internal Thoracic Artery donor2'],
        446:['CNhs11989','Smooth Muscle Cells - Pulmonary Artery donor2'],
        447:['CNhs11990','Smooth Muscle Cells - Subclavian Artery donor2'],
        448:['CNhs11991','Smooth Muscle Cells - Umbilical Artery donor2'],
        449:['CNhs11992','Synoviocyte donor2'],
        450:['CNhs11993','Tracheal Epithelial Cells donor2'],
        451:['CNhs11996','Fibroblast - Periodontal Ligament donor6 (PLH3)'],
        452:['CNhs11997','CD14+ Monocytes donor3'],
        453:['CNhs11998','CD4+ T Cells donor3'],
        454:['CNhs11999','CD8+ T Cells donor3'],
        455:['CNhs12000','Dendritic Cells - monocyte immature derived donor3'],
        456:['CNhs12001','Natural Killer Cells donor3'],
        457:['CNhs12002','Peripheral Blood Mononuclear Cells donor3'],
        458:['CNhs12003','Macrophage - monocyte derived donor3'],
        459:['CNhs12004','Smooth Muscle Cells - Brain Vascular donor3'],
        460:['CNhs12005','Astrocyte - cerebral cortex donor3'],
        461:['CNhs12006','Fibroblast - Gingival donor3'],
        462:['CNhs12007','Smooth Muscle Cells - Colonic donor3'],
        463:['CNhs12008','Skeletal Muscle Satellite Cells donor3'],
        464:['CNhs12009','Ciliary Epithelial Cells donor3'],
        465:['CNhs12010','Endothelial Cells - Umbilical vein donor3'],
        466:['CNhs12011','Fibroblast - Aortic Adventitial donor3'],
        467:['CNhs12012','Mesothelial Cells donor3'],
        468:['CNhs12013','Preadipocyte - omental donor3'],
        469:['CNhs12014','Prostate Epithelial Cells donor3'],
        470:['CNhs12015','Prostate Stromal Cells donor3'],
        471:['CNhs12016','Small Airway Epithelial Cells donor3'],
        472:['CNhs12017','Adipocyte - subcutaneous donor3'],
        473:['CNhs12019','Nucleus Pulposus Cell donor2'],
        474:['CNhs12020','Chondrocyte - de diff donor3'],
        475:['CNhs12021','Chondrocyte - re diff donor3'],
        476:['CNhs12022','Endothelial Cells - Aortic donor3'],
        477:['CNhs12023','Endothelial Cells - Artery donor3'],
        478:['CNhs12024','Endothelial Cells - Microvascular donor3'],
        479:['CNhs12026','Endothelial Cells - Vein donor3'],
        480:['CNhs12027','Fibroblast - Cardiac donor3'],
        481:['CNhs12028','Fibroblast - Dermal donor3'],
        482:['CNhs12030','Hair Follicle Dermal Papilla Cells donor3'],
        483:['CNhs12031','Keratinocyte - epidermal donor3'],
        484:['CNhs12032','Mammary Epithelial Cell donor3'],
        485:['CNhs12033','Melanocyte - light donor3'],
        486:['CNhs12035','Osteoblast - differentiated donor3'],
        487:['CNhs12036','Osteoblast donor3'],
        488:['CNhs12037','Placental Epithelial Cells donor3'],
        489:['CNhs12038','Preadipocyte - subcutaneous donor3'],
        490:['CNhs12039','Preadipocyte - visceral donor3'],
        491:['CNhs12043','Smooth Muscle Cells - Brachiocephalic donor3'],
        492:['CNhs12044','Smooth Muscle Cells - Carotid donor3'],
        493:['CNhs12045','Smooth Muscle Cells - Coronary Artery donor3'],
        494:['CNhs12046','Smooth Muscle Cells - Internal Thoracic Artery donor3'],
        495:['CNhs12048','Smooth Muscle Cells - Subclavian Artery donor3'],
        496:['CNhs12049','Smooth Muscle Cells - Umbilical Artery donor3'],
        497:['CNhs12050','Synoviocyte donor3'],
        498:['CNhs12051','Tracheal Epithelial Cells donor3'],
        499:['CNhs12052','Fibroblast - Dermal donor4'],
        500:['CNhs12053','Skeletal Muscle Cells donor4'],
        501:['CNhs12054','Bronchial Epithelial Cell donor4'],
        502:['CNhs12055','Fibroblast - Dermal donor5'],
        503:['CNhs12056','Skeletal Muscle Cells donor5'],
        504:['CNhs12057','Fibroblast - Cardiac donor5'],
        505:['CNhs12058','Bronchial Epithelial Cell donor5'],
        506:['CNhs12059','Fibroblast - Dermal donor6'],
        507:['CNhs12060','Skeletal Muscle Cells donor6'],
        508:['CNhs12061','Fibroblast - Cardiac donor6'],
        509:['CNhs12062','Bronchial Epithelial Cell donor6'],
        510:['CNhs12063','Nucleus Pulposus Cell donor3'],
        511:['CNhs12064','Anulus Pulposus Cell donor2'],
        512:['CNhs12065','Preadipocyte - perirenal donor1'],
        513:['CNhs12067','Adipocyte - omental donor2'],
        514:['CNhs12068','Adipocyte - omental donor3'],
        515:['CNhs12069','Adipocyte - perirenal donor1'],
        516:['CNhs12074','Renal Glomerular Endothelial Cells donor1'],
        517:['CNhs12075','Hepatic Sinusoidal Endothelial Cells donor1'],
        518:['CNhs12079','Pericytes donor2'],
        519:['CNhs12080','Meningeal Cells donor2'],
        520:['CNhs12081','Astrocyte - cerebellum donor2'],
        521:['CNhs12084','Alveolar Epithelial Cells donor2'],
        522:['CNhs12086','Renal Glomerular Endothelial Cells donor2'],
        523:['CNhs12087','Renal Proximal Tubular Epithelial Cell donor2'],
        524:['CNhs12088','Renal Epithelial Cells donor2'],
        525:['CNhs12091','Urothelial Cells donor2'],
        526:['CNhs12092','Hepatic Sinusoidal Endothelial Cells donor2'],
        527:['CNhs12093','Hepatic Stellate Cells (lipocyte) donor2'],
        528:['CNhs12095','Keratocytes donor2'],
        529:['CNhs12100','Mesenchymal Stem Cells - bone marrow donor2'],
        530:['CNhs12104','Mesenchymal Stem Cells - amniotic membrane donor2'],
        531:['CNhs12105','Multipotent Cord Blood Unrestricted Somatic Stem Cells donor2'],
        532:['CNhs12117','Astrocyte - cerebellum donor3'],
        533:['CNhs12118','Fibroblast - Lymphatic donor3'],
        534:['CNhs12120','Renal Proximal Tubular Epithelial Cell donor3'],
        535:['CNhs12121','Renal Mesangial Cells donor3'],
        536:['CNhs12122','Urothelial Cells donor3'],
        537:['CNhs12123','Corneal Epithelial Cells donor3'],
        538:['CNhs12124','Trabecular Meshwork Cells donor3'],
        539:['CNhs12125','Amniotic Epithelial Cells donor3'],
        540:['CNhs12126','Mesenchymal Stem Cells - bone marrow donor3'],
        541:['CNhs12127','Mesenchymal Stem Cells - umbilical donor3'],
        542:['CNhs12227','spinal cord adult donor10252'],
        543:['CNhs12228','pineal gland adult donor10252'],
        544:['CNhs12229','pituitary gland adult donor10252'],
        545:['CNhs12310','medial temporal gyrus adult donor10252'],
        546:['CNhs12311','amygdala adult donor10252'],
        547:['CNhs12312','hippocampus adult donor10252'],
        548:['CNhs12314','thalamus adult donor10252'],
        549:['CNhs12315','medulla oblongata adult donor10252'],
        550:['CNhs12316','middle temporal gyrus donor10252'],
        551:['CNhs12317','parietal lobe adult donor10252'],
        552:['CNhs12318','substantia nigra adult donor10252'],
        553:['CNhs12319','globus pallidus adult donor10252'],
        554:['CNhs12320','occipital cortex adult donor10252'],
        555:['CNhs12321','caudate nucleus adult donor10252'],
        556:['CNhs12322','locus coeruleus adult donor10252'],
        557:['CNhs12323','cerebellum adult donor10252'],
        558:['CNhs12324','putamen adult donor10196'],
        559:['CNhs12325','epitheloid carcinoma cell line: HelaS3 ENCODE biol_rep1'],
        560:['CNhs12326','epitheloid carcinoma cell line: HelaS3 ENCODE biol_rep2'],
        561:['CNhs12327','epitheloid carcinoma cell line: HelaS3 ENCODE biol_rep3'],
        562:['CNhs12328','hepatocellular carcinoma cell line: HepG2 ENCODE biol_rep1'],
        563:['CNhs12329','hepatocellular carcinoma cell line: HepG2 ENCODE biol_rep2'],
        564:['CNhs12330','hepatocellular carcinoma cell line: HepG2 ENCODE biol_rep3'],
        565:['CNhs12331','B lymphoblastoid cell line: GM12878 ENCODE biol_rep1'],
        566:['CNhs12332','B lymphoblastoid cell line: GM12878 ENCODE biol_rep2'],
        567:['CNhs12333','B lymphoblastoid cell line: GM12878 ENCODE biol_rep3'],
        568:['CNhs12334','chronic myelogenous leukemia cell line:K562 ENCODE biol_rep1'],
        569:['CNhs12335','chronic myelogenous leukemia cell line:K562 ENCODE biol_rep2'],
        570:['CNhs12336','chronic myelogenous leukemia cell line:K562 ENCODE biol_rep3'],
        571:['CNhs12338','Neurons donor1'],
        572:['CNhs12339','Hair Follicle Outer Root Sheath Cells donor1'],
        573:['CNhs12340','Hepatocyte donor1'],
        574:['CNhs12341','Cardiac Myocyte donor1'],
        575:['CNhs12342','Lens Epithelial Cells donor1'],
        576:['CNhs12343','CD19+ B Cells donor1'],
        577:['CNhs12344','Fibroblast - Choroid Plexus donor2'],
        578:['CNhs12347','Hair Follicle Outer Root Sheath Cells donor2'],
        579:['CNhs12348','Smooth Muscle Cells - Bronchial donor2'],
        580:['CNhs12349','Hepatocyte donor2'],
        581:['CNhs12350','Cardiac Myocyte donor2'],
        582:['CNhs12351','testicular germ cell embryonal carcinoma cell line:NEC14'],
        583:['CNhs12352','CD19+ B Cells donor2'],
        584:['CNhs12354','CD19+ B Cells donor3'],
        585:['CNhs12362','testicular germ cell embryonal carcinoma cell line:NEC15'],
        586:['CNhs12363','mesenchymal precursor cell - adipose donor1'],
        587:['CNhs12364','mesenchymal precursor cell - adipose donor2'],
        588:['CNhs12365','mesenchymal precursor cell - adipose donor3'],
        589:['CNhs12366','mesenchymal precursor cell - bone marrow donor1'],
        590:['CNhs12367','mesenchymal precursor cell - bone marrow donor2'],
        591:['CNhs12368','mesenchymal precursor cell - cardiac donor1'],
        592:['CNhs12369','mesenchymal precursor cell - cardiac donor2'],
        593:['CNhs12370','mesenchymal precursor cell - cardiac donor3'],
        594:['CNhs12371','mesenchymal precursor cell - cardiac donor4'],
        595:['CNhs12372','mesenchymal precursor cell - ovarian cancer left ovary donor1'],
        596:['CNhs12373','mesenchymal precursor cell - ovarian cancer right ovary donor1'],
        597:['CNhs12374','mesenchymal precursor cell - ovarian cancer metastasis donor1'],
        598:['CNhs12375','mesenchymal precursor cell - ovarian cancer right ovary donor2'],
        599:['CNhs12376','mesenchymal precursor cell - ovarian cancer left ovary donor3'],
        600:['CNhs12377','mesenchymal precursor cell - ovarian cancer right ovary donor3 (SOC-57-02)'],
        601:['CNhs12378','mesenchymal precursor cell - ovarian cancer metastasis donor3'],
        602:['CNhs12379','amniotic membrane cells donor3'],
        603:['CNhs12380','chorionic membrane cells donor3'],
        604:['CNhs12492','Mesenchymal stem cells - umbilical donor0'],
        605:['CNhs12493','Fibroblast - Periodontal Ligament donor4 (PL29)'],
        606:['CNhs12494','Adipocyte - subcutaneous donor1'],
        607:['CNhs12495','Endothelial Cells - Aortic donor1'],
        608:['CNhs12496','Endothelial Cells - Artery donor1'],
        609:['CNhs12497','Endothelial Cells - Vein donor1'],
        610:['CNhs12498','Fibroblast - Cardiac donor1'],
        611:['CNhs12499','Fibroblast - Cardiac donor1'],
        612:['CNhs12501','Hair Follicle Dermal Papilla Cells donor1'],
        613:['CNhs12502','amniotic membrane cells donor1'],
        614:['CNhs12503','amniotic membrane cells donor2'],
        615:['CNhs12504','chorionic membrane cells donor1'],
        616:['CNhs12506','chorionic membrane cells donor2'],
        617:['CNhs12566','Mast cell donor1'],
        618:['CNhs12568','Lens Epithelial Cells donor2'],
        619:['CNhs12569','Smooth Muscle Cells - Umbilical Vein donor2'],
        620:['CNhs12570','Melanocyte - dark donor3'],
        621:['CNhs12571','Cardiac Myocyte donor3'],
        622:['CNhs12572','Lens Epithelial Cells donor3'],
        623:['CNhs12574','nasal epithelial cells donor2'],
        624:['CNhs12575','Basophils donor3'],
        625:['CNhs12588','CD34+ stem cells - adult bone marrow derived donor1 tech_rep1'],
        626:['CNhs12589','nasal epithelial cells donor1 tech_rep1'],
        627:['CNhs12592','Mast cell donor4'],
        628:['CNhs12593','Mast cell donor3'],
        629:['CNhs12594','Mast cell donor2'],
        630:['CNhs12596','Iris Pigment Epithelial Cells donor1'],
        631:['CNhs12597','Smooth Muscle Cells - Umbilical Vein donor1'],
        632:['CNhs12610','diencephalon adult'],
        633:['CNhs12611','olfactory region adult'],
        634:['CNhs12624','Renal Glomerular Endothelial Cells donor3'],
        635:['CNhs12626','Hepatocyte donor3'],
        636:['CNhs12639','tenocyte donor1'],
        637:['CNhs12640','tenocyte donor2'],
        638:['CNhs12641','tenocyte donor3'],
        639:['CNhs12726','Neurons donor2'],
        640:['CNhs12728','Renal Cortical Epithelial Cells donor2'],
        641:['CNhs12730','Mesenchymal Stem Cells - hepatic donor2'],
        642:['CNhs12731','Meningeal Cells donor3'],
        643:['CNhs12732','Renal Epithelial Cells donor3'],
        644:['CNhs12733','Retinal Pigment Epithelial Cells donor3'],
        645:['CNhs12805','medulloblastoma cell line:D283 Med'],
        646:['CNhs12806','large cell lung carcinoma cell line:NCI-H460'],
        647:['CNhs12807','plasma cell leukemia cell line:ARH-77'],
        648:['CNhs12808','small cell lung carcinoma cell line:DMS 144'],
        649:['CNhs12809','small cell lung carcinoma cell line:NCI-H82'],
        650:['CNhs12810','salivary acinar cells donor1'],
        651:['CNhs12811','salivary acinar cells donor2'],
        652:['CNhs12812','salivary acinar cells donor3'],
        653:['CNhs12838','merkel cell carcinoma cell line:MKL-1'],
        654:['CNhs12839','merkel cell carcinoma cell line:MS-1'],
        655:['CNhs12840','cerebral meninges adult'],
        656:['CNhs12842','appendix adult'],
        657:['CNhs12844','vein adult'],
        658:['CNhs12846','ductus deferens adult'],
        659:['CNhs12847','epididymis adult'],
        660:['CNhs12848','gall bladder adult'],
        661:['CNhs12849','parotid gland adult'],
        662:['CNhs12850','penis adult'],
        663:['CNhs12851','seminal vesicle adult'],
        664:['CNhs12852','submaxillary gland adult'],
        665:['CNhs12853','tongue adult'],
        666:['CNhs12854','vagina adult'],
        667:['CNhs12855','heart - mitral valve adult'],
        668:['CNhs12856','heart - pulmonic valve adult'],
        669:['CNhs12857','heart - tricuspid valve adult'],
        670:['CNhs12858','throat adult'],
        671:['CNhs12894','Smooth Muscle Cells - Tracheal donor3'],
        672:['CNhs12922','Mesenchymal Stem Cells - adipose donor3'],
        673:['CNhs12996','duodenum fetal donor1 tech_rep2'],
        674:['CNhs12997','temporal lobe fetal donor1 tech_rep2'],
        675:['CNhs12998','testis adult pool2'],
        676:['CNhs13049','acute myeloid leukemia (FAB M7) cell line:M-MOK'],
        677:['CNhs13050','acute myeloid leukemia (FAB M5) cell line:NOMO-1'],
        678:['CNhs13051','acute myeloid leukemia (FAB M5) cell line:P31/FUJ'],
        679:['CNhs13052','acute myeloid leukemia (FAB M2) cell line:Kasumi-6'],
        680:['CNhs13053','acute myeloid leukemia (FAB M0) cell line:KG-1'],
        681:['CNhs13054','acute myeloid leukemia (FAB M1) cell line:HYT-1'],
        682:['CNhs13055','acute myeloid leukemia (FAB M3) cell line:HL60'],
        683:['CNhs13056','acute myeloid leukemia (FAB M4eo) cell line:EoL-1'],
        684:['CNhs13057','acute myeloid leukemia (FAB M4eo) cell line:EoL-3'],
        685:['CNhs13058','acute myeloid leukemia (FAB M5) cell line:U-937 DE-4'],
        686:['CNhs13059','acute myeloid leukemia (FAB M6) cell line:EEB'],
        687:['CNhs13060','acute myeloid leukemia (FAB M6) cell line:F-36E'],
        688:['CNhs13061','mesothelioma cell line:NCI-H28'],
        689:['CNhs13062','mesothelioma cell line:NCI-H226'],
        690:['CNhs13063','mesothelioma cell line:NCI-H2052'],
        691:['CNhs13064','mesothelioma cell line:NCI-H2452'],
        692:['CNhs13066','mesothelioma cell line:Mero-25'],
        693:['CNhs13067','mesothelioma cell line:Mero-41'],
        694:['CNhs13068','mesothelioma cell line:Mero-48a'],
        695:['CNhs13069','mesothelioma cell line:Mero-82'],
        696:['CNhs13070','mesothelioma cell line:Mero-83'],
        697:['CNhs13072','mesothelioma cell line:Mero-84'],
        698:['CNhs13073','mesothelioma cell line:Mero-95'],
        699:['CNhs13074','mesothelioma cell line:No36'],
        700:['CNhs13075','mesothelioma cell line:ONE58'],
        701:['CNhs13080','Renal Glomerular Endothelial Cells donor4'],
        702:['CNhs13092','mesenchymal precursor cell - ovarian cancer left ovary donor2'],
        703:['CNhs13093','mesenchymal precursor cell - ovarian cancer metastasis donor2'],
        704:['CNhs13094','mesenchymal precursor cell - ovarian cancer left ovary donor4'],
        705:['CNhs13096','mesenchymal precursor cell - ovarian cancer right ovary donor4'],
        706:['CNhs13097','mesenchymal precursor cell - ovarian cancer metastasis donor4'],
        707:['CNhs13098','mesenchymal precursor cell - bone marrow donor3'],
        708:['CNhs13099','serous adenocarcinoma cell line:SK-OV-3-R biol_rep1'],
        709:['CNhs13195','CD4+CD25+CD45RA- memory regulatory T cells donor1'],
        710:['CNhs13202','CD4+CD25-CD45RA+ naive conventional T cells expanded donor1'],
        711:['CNhs13203','CD4+CD25+CD45RA+ naive regulatory T cells expanded donor1'],
        712:['CNhs13204','CD4+CD25+CD45RA- memory regulatory T cells expanded donor1'],
        713:['CNhs13205','CD4+CD25-CD45RA+ naive conventional T cells donor2'],
        714:['CNhs13206','CD4+CD25+CD45RA- memory regulatory T cells donor2'],
        715:['CNhs13207','CD14-CD16+ Monocytes donor2'],
        716:['CNhs13208','CD14+CD16+ Monocytes donor2'],
        717:['CNhs13215','CD4+CD25-CD45RA- memory conventional T cells expanded donor1'],
        718:['CNhs13216','CD14+CD16- Monocytes donor2'],
        719:['CNhs13223','CD4+CD25-CD45RA+ naive conventional T cells donor1'],
        720:['CNhs13224','CD14+CD16- Monocytes donor1'],
        721:['CNhs13449','optic nerve donor1'],
        722:['CNhs13454','skeletal muscle - soleus muscle donor1'],
        723:['CNhs13465','CD14+ monocytes - treated with BCG donor1'],
        724:['CNhs13466','CD14+ monocytes - treated with IFN + N-hexane donor1'],
        725:['CNhs13467','CD14+ monocytes - treated with Trehalose dimycolate (TDM) donor1'],
        726:['CNhs13468','CD14+ monocytes - mock treated donor1'],
        727:['CNhs13469','CD14+ monocytes - treated with Group A streptococci donor1'],
        728:['CNhs13470','CD14+ monocytes - treated with lipopolysaccharide donor1'],
        729:['CNhs13471','CD14+ monocytes - treated with Salmonella donor1'],
        730:['CNhs13472','CD14+ monocytes - treated with Cryptococcus donor1'],
        731:['CNhs13473','CD14+ monocytes - treated with Candida donor1'],
        732:['CNhs13474','CD14+ monocytes - treated with B-glucan donor1'],
        733:['CNhs13475','CD14+ monocytes - treated with BCG donor2'],
        734:['CNhs13476','CD14+ monocytes - treated with IFN + N-hexane donor2'],
        735:['CNhs13477','Hep-2 cells treated with Streptococci strain 5448 biol_rep1'],
        736:['CNhs13478','Hep-2 cells treated with Streptococci strain JRS4 biol_rep1'],
        737:['CNhs13479','Hep-2 cells mock treated biol_rep1'],
        738:['CNhs13480','immature langerhans cells donor2'],
        739:['CNhs13483','CD14+ monocytes - treated with Trehalose dimycolate (TDM) donor2'],
        740:['CNhs13484','CD14+ monocytes - mock treated donor2'],
        741:['CNhs13485','CD14+ monocytes - treated with Salmonella donor2'],
        742:['CNhs13487','CD14+ monocytes - treated with Cryptococcus donor2'],
        743:['CNhs13488','CD14+ monocytes - treated with Candida donor2'],
        744:['CNhs13489','CD14+ monocytes - treated with B-glucan donor2'],
        745:['CNhs13490','CD14+ monocytes - treated with IFN + N-hexane donor3'],
        746:['CNhs13491','CD14+ monocytes - mock treated donor3'],
        747:['CNhs13492','CD14+ monocytes - treated with Group A streptococci donor3'],
        748:['CNhs13493','CD14+ monocytes - treated with Salmonella donor3'],
        749:['CNhs13494','CD14+ monocytes - treated with Candida donor3'],
        750:['CNhs13495','CD14+ monocytes - treated with B-glucan donor3'],
        751:['CNhs13496','Hep-2 cells treated with Streptococci strain 5448 biol_rep2'],
        752:['CNhs13497','Hep-2 cells treated with Streptococci strain 5448 biol_rep3'],
        753:['CNhs13498','Hep-2 cells treated with Streptococci strain JRS4 biol_rep2'],
        754:['CNhs13499','Hep-2 cells treated with Streptococci strain JRS4 biol_rep3'],
        755:['CNhs13500','Hep-2 cells mock treated biol_rep2'],
        756:['CNhs13501','Hep-2 cells mock treated biol_rep3'],
        757:['CNhs13502','acute myeloid leukemia (FAB M2) cell line:Kasumi-1'],
        758:['CNhs13503','acute myeloid leukemia (FAB M4) cell line:FKH-1'],
        759:['CNhs13504','acute myeloid leukemia (FAB M4) cell line:HNT-34'],
        760:['CNhs13505','acute myeloid leukemia (FAB M6) cell line:F-36P'],
        761:['CNhs13507','mesenchymal precursor cell - ovarian cancer right ovary donor3 (SOC-57-02-G)'],
        762:['CNhs13508','serous adenocarcinoma cell line:SK-OV-3-R after co-culture with SOC-57-02-G biol_rep1'],
        763:['CNhs13512','CD4+CD25-CD45RA+ naive conventional T cells donor3'],
        764:['CNhs13513','CD4+CD25+CD45RA+ naive regulatory T cells donor3'],
        765:['CNhs13532','CD14+ monocytes - treated with Group A streptococci donor2'],
        766:['CNhs13533','CD14+ monocytes - treated with lipopolysaccharide donor2'],
        767:['CNhs13535','migratory langerhans cells donor1'],
        768:['CNhs13536','migratory langerhans cells donor2'],
        769:['CNhs13537','immature langerhans cells donor1'],
        770:['CNhs13538','CD4+CD25+CD45RA- memory regulatory T cells donor3'],
        771:['CNhs13539','CD4+CD25-CD45RA- memory conventional T cells donor3'],
        772:['CNhs13540','CD14+CD16- Monocytes donor3'],
        773:['CNhs13541','CD14+CD16+ Monocytes donor1'],
        774:['CNhs13543','CD14+ monocytes - treated with BCG donor3'],
        775:['CNhs13544','CD14+ monocytes - treated with Trehalose dimycolate (TDM) donor3'],
        776:['CNhs13545','CD14+ monocytes - treated with lipopolysaccharide donor3'],
        777:['CNhs13546','CD14+ monocytes - treated with Cryptococcus donor3'],
        778:['CNhs13547','migratory langerhans cells donor3'],
        779:['CNhs13548','CD14-CD16+ Monocytes donor3'],
        780:['CNhs13549','CD14+CD16+ Monocytes donor3'],
        781:['CNhs13550','Mallassez-derived cells donor2'],
        782:['CNhs13551','Mallassez-derived cells donor3'],
        783:['CNhs13552','Reticulocytes biol_ rep1'],
        784:['CNhs13553','Reticulocytes biol_ rep2'],
        785:['CNhs13793','amygdala - adult donor10196'],
        786:['CNhs13794','thalamus - adult donor10196'],
        787:['CNhs13795','hippocampus - adult donor10196'],
        788:['CNhs13796','medial frontal gyrus - adult donor10196'],
        789:['CNhs13797','parietal lobe - adult donor10196'],
        790:['CNhs13798','occipital cortex - adult donor10196'],
        791:['CNhs13799','cerebellum - adult donor10196'],
        792:['CNhs13800','medulla oblongata - adult donor10196'],
        793:['CNhs13801','globus pallidus - adult donor10196'],
        794:['CNhs13802','caudate nucleus - adult donor10196'],
        795:['CNhs13804','pineal gland - adult donor10196'],
        796:['CNhs13805','pituitary gland - adult donor10196'],
        797:['CNhs13807','spinal cord - adult donor10196'],
        798:['CNhs13808','locus coeruleus - adult donor10196'],
        799:['CNhs13809','medial temporal gyrus - adult donor10196'],
        800:['CNhs13811','CD4+CD25+CD45RA- memory regulatory T cells expanded donor2'],
        801:['CNhs13812','CD4+CD25+CD45RA- memory regulatory T cells expanded donor3'],
        802:['CNhs13813','CD4+CD25-CD45RA+ naive conventional T cells expanded donor2'],
        803:['CNhs13814','CD4+CD25-CD45RA+ naive conventional T cells expanded donor3'],
        804:['CNhs13815','Neurons donor3'],
        805:['CNhs13816','Olfactory epithelial cells donor1'],
        806:['CNhs13817','Olfactory epithelial cells donor2'],
        807:['CNhs13818','Olfactory epithelial cells donor3'],
        808:['CNhs13819','Olfactory epithelial cells donor4'],
    }
    result = []
    for snp in snps:
        for t in tissues:
        ##find tissues based on input
            if t[:1] == 'E':
                for elemt in roadmap:
                    elemt = int(elemt)
                    for i in roadmap[elemt]:
                        if i == t:
                            tissue_selected = elemt
                            tissue_name = roadmap[elemt][1]
                            tissue_id = roadmap[elemt][0]
            elif t[:1] == 'C':
                for elemt in fantom_5:
                    elemt = int(elemt)
                    for i in fantom_5[elemt]:
                        if i == t:
                            tissue_selected = elemt
                            tissue_name = fantom_5[elemt][1]
                            tissue_id = fantom_5[elemt][0]
            result.append( epi_search( snp[1], snp[2], snp[0],
                                       str(tissue_selected),
                                       tissue_name,tissue_id))

    result2 = pd.DataFrame(columns= ['SNP','Location','Gene','Score','Tissue','File_Type'])
    result2 = pd.concat(result)
    result3 = result2.drop_duplicates(subset=['SNP','Location','Gene','Score','Tissue','File_Type'],
                                        keep='first', inplace=False)
    new = result3["Gene"].str.split("$", n = 50, expand = True) 
    print("**---RESULT3---**")
    print(result3)
    print("**---NEW---**")
    print(new)
    if not(result3.empty):
        # making seperate first name column from new data frame 
        result3["Gene ID"]= new[0]
        # making seperate last name column from new data frame 
        result3["Gene"]= new[1]
        result3["Gene Location"]= new[2]+ "-" + new[3].map(str)
        if 4 in new:
            new2 = new[4].str.split("%", n = 50, expand = True)
            gene_id = result3["Gene ID"]
            print(new2)
            if 1 in new2:
                list_new_id = new2[1]
                new_gene_id = []
                for i,j in zip(gene_id,list_new_id):
                    if j is None:
                        new_gene_id.append(i)
                    else:
                        new_gene_id.append(i+'|'+j)
                result3["Gene ID"] = new_gene_id
            result3["Gene Location"] = result3["Gene Location"]+ "" + new2[0].map(str)
            result3 = result3[['SNP', 'Location', 'Gene ID', 'Gene', 'Gene Location', 'Score', 'Tissue', 'File_Type']]
    result3 = result3.reset_index().drop(columns=['index'])
    # result3 = result3.to_json()
    user_id = request.form['user_id']
    file_name = str(user_id)+'_step5.txt'
    save = result3.to_json()
    # save = result3
    # with open('app/users_workflow/'+user_id+'_step5.txt', 'w') as out:  
    #         json.dump(save, out)
    # workflow = Workflow()
    # user_work = Workflow.query.filter_by(user_id_user=user_id).first()
    # user_work.step5 = 'app/users_workflow/'+file_name
    # db.session.commit()
    return result3


@blueprint.route('/index',methods=['GET','POST'])
@cross_origin(origin='*')
@login_required
def index():
    return render_template('index.html')

@blueprint.route('/<template>')
@cross_origin(origin='*')
@login_required
def route_template(template):
    return render_template(template + '.html')

#get only one snp
@blueprint.route('/get_snp_info',methods=['GET','POST'])
@cross_origin(origin='*')
@login_required
def teste():
    #if request method is post execute
    if request.method == 'POST':
        #minor allele list
        minor_allele = []
        #get value from snp_name field from html form
        snp_input_form = str(request.form['snp_name'])
        genome_input = "GRCh37.19"
        snp = snp_input_form[2:]
        a = get_snp_info(snp)

        #dictionary with information based on genome version
        sample_dict = a[a["gnenome_versions"].str.contains(genome_input[:6])]['snp_info_dict'].values[0]

        #if sample dict is bigger than 1(more than one allele in variation)
        if len(sample_dict) > 1:
            for dic in sample_dict:
                minor_allele.append(dic['allele_v'])
            #Stringfy minor allele
            join_symbol = "|"
            str_minor = join_symbol.join(minor_allele)
        else:
            str_minor = sample_dict[0]['allele_v']
        
        return jsonify(sample_dict,snp_input_form,str_minor)
#step 2
@blueprint.route('/verify_snps',methods=['GET','POST'])
@cross_origin(origin='*')
@login_required
def verify_snps():
    dict_snps = []
    #get snp list information from the datatable
    snp_list_rows = request.form['snp_list']
    #transform it in a json file(dictionary?)
    snp_list_rows = json.loads(snp_list_rows)
    if request.method == 'POST':

        user_email = request.form['user_email']

        tissue_list = request.form['tissue_field'].split("|")

        user_id = request.form['user_id']

        # Start threaded verify snps
        # @copy_current_request_context
        # def execute_verify_task(dict_snps, user_id, user_email, tissue_list,snp_list_rows):
        #     verify_snps_threaded(dict_snps, user_id, user_email, tissue_list, snp_list_rows)

        task_res = verify_snps_threaded.apply_async(args=[dict_snps, user_id, user_email, tissue_list, snp_list_rows])
        
        # thread_name = "step2_user_"+user_id
        # threading.Thread(name = thread_name,target=execute_verify_task,
        #                     args=(dict_snps, user_id, user_email, tissue_list, snp_list_rows)).start()
        # funcao_teste_snp(request)
        # return jsonify(resultado = {}, status=202)
        return jsonify({}), 202, {'Location': url_for('home_blueprint.taskstatus2',task_id=task_res.id)}

# step 2 threaded
@celery.task(bind=True)
def verify_snps_threaded(self,dict_snps, user_id, user_email, tissue_list,snp_list_rows):

    # Remove from directory
    user_work = Workflow.query.filter_by(user_id_user=user_id).first()
    if os.path.exists("app/users_workflow/"+user_id+"_step2.txt"):
        os.remove("app/users_workflow/"+user_id+"_step2.txt")

    # # Removing from database
    user_work.step2 = None
    db.session.commit()

    for idx1, tissue in enumerate(tissue_list):
        #iterate through every snp
        for idx2, snp_info in enumerate(snp_list_rows):
            #change chromossome type from number to 'X' or 'Y' string
            if str(snp_info[2]) == '23':
                chrom = 'X'
            elif str(snp_info[2]) == '24':
                chrom = 'Y'
            else:
                chrom = snp_info[2]

            dict_snps.append(apply_state_model(tissue, snp_info, snp_info[0], chrom))
            self.update_state(state='PROGRESS',
                                meta={'current': idx2+1, 'total': len(snp_list_rows),
                                'status': "TISSUE {0}: {1}, id: {2}, is done being processed from a total of {3}".format(tissue,snp_info[0],str(idx2+1),str(len(snp_list_rows)))})
        self.update_state(state='PROGRESS',
                          meta={'current': idx1+1, 
                                'total': len(tissue_list),
                                'status': "TISSUE {0},id: {1}, is done being processed from a total of {2}".format(tissue,str(idx1+1),str(len(tissue_list)))})

    self.update_state(state='DONE',
                          meta={'current': len(tissue_list)*len(snp_list_rows), 
                                'total': len(tissue_list)*len(snp_list_rows),
                                'status': "All SNPs are done. Results should appear shortly"})
    #Add to workflow in database
    # user_id = request.form['user_id']

    file_name = str(user_id)+'_step2.txt'
    with open('app/users_workflow/'+user_id+'_step2.txt', 'w') as out:  
        json.dump(dict_snps, out)
    workflow = Workflow()
    user_work = Workflow.query.filter_by(user_id_user=user_id).first()
    user_work.step2 = 'app/users_workflow/'+file_name

    #delete other steps forward
    user_work.step3 = None
    user_work.step4 = None
    user_work.step5 = None
    if os.path.exists("app/users_workflow/"+user_id+"_step3.txt"):
        os.remove("app/users_workflow/"+user_id+"_step3.txt")
    if os.path.exists("app/users_workflow/"+user_id+"_step3_dictionary.txt"):
        os.remove("app/users_workflow/"+user_id+"_step3_dictionary.txt")
    if os.path.exists("app/users_workflow/"+user_id+"_step4.txt"):
        os.remove("app/users_workflow/"+user_id+"_step4.txt")
    if os.path.exists("app/users_workflow/"+user_id+"_step5.txt"):
        os.remove("app/users_workflow/"+user_id+"_step5.txt")
    
    db.session.commit()

    # user_email = request.form['user_email']

    try:
        server = smtplib.SMTP_SSL('smtp.gmail.com',465)
        server.login("regulomix.temp@gmail.com","regulomix123")
        subject = "Regulomix: Step2"
        body = "Step2 work is done, login to regulomix to continue."
        message = "Subject:{}\n\n{}".format(subject, body)
        server.sendmail("regulomix.temp@gmail.com",user_email,message)
        server.quit()
        print("Email sent")
    except:
        server.quit()
        print("Email failed to send")

    # return jsonify(dict_snps)
    return {'current': 100, 'total': 100, 'status': 'All completed','result': dict_snps}

@blueprint.route('/taskstatus2/<task_id>',methods=['GET','POST'])
@cross_origin(origin='*')
@login_required
def taskstatus2(task_id):
    task = verify_snps_threaded.AsyncResult(task_id)
    if task.state == 'PENDING':
        response = {
            'state': task.state,
            'current': 0,
            'total': 1,
            'status': 'Pending...'
        }
    elif task.state != 'FAILURE':
        response = {
            'state': task.state,
            'current': task.info.get('current', 0),
            'total': task.info.get('total', 1),
            'status': task.info.get('status', '')
        }
        if 'result' in task.info:
            response['result'] = task.info['result']
    else:
        # something went wrong in the background job
        response = {
            'state': task.state,
            'current': 1,
            'total': 1,
            'status': str(task.info),  # this is the exception raised
        }
    return jsonify(response)

#step 3
@blueprint.route('/gen_sequence',methods=['GET','POST'])
@cross_origin(origin='*')
@login_required
def gen_sequence():

    #TODO this method is needed for step 4
    def create_alleles_dictionary(snp_names,info_list):
        dict_snp_allele = {}
        tuple_al = []
        
        index = 0
        #iterate through names and all dictionaies in snp
        for i,snp_name in zip(info_list,snp_names):
            # snp_al = i[index]['allele_wt']
            snp_al = i[3]
            tuple_al.append(snp_al)
            list_minor = []
            minor_allele = i[4].split('|')
            for al in minor_allele:
                tuple_al.append(al)

            list_minor=[]

            dict_snp_allele[snp_name] = tuple_al 
            tuple_al = []

        return dict_snp_allele

    #if request method is post execute
    if request.method == 'POST':
        #get snp list information from the datatable
        filtered_snp = json.loads(request.form['snp_list'])
        #hardcoded genome version
        gnenome_version = 'GRCh37.p13'
        #user selected transcription factor matrix
        tf_selected = request.form['tf']
        #get p-value selected
        p_value = request.form['p_value']
        # checkbox q-value
        q_value = request.form['q_value']
        #get snps ids from filtered list (stage 2)
        snp_ids_list = []
        for x in filtered_snp:
            snp_ids_list.append(x[0])
        
        #only unique values
        snp_ids_list_unique = list(set(snp_ids_list))

        #download snp info
        info_snp_list = []

        snp_names = []
        snp_dicts = []

        for snp_info in (filtered_snp):
        
            snp_name = snp_info[0]
            sample_dict = snp_info

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
        #FIMO
        upload_value = request.form['upload']
        if (int(upload_value) == 0):
            meme= "./TFs/" + tf_selected + ".meme"
            print("not uploaded")
        else:
            meme= "app/home/meme.meme"
            print("uploaded")
        # meme= "./TFs/" + tf_selected + ".meme"
        os.system("export PATH=$HOME/meme/bin:$PATH ")
        if (q_value == "true"):
            os.system("~/meme/bin/fimo --thresh "+ p_value +" --qv-thresh "+meme+" "+fna_filename )
        else:
            os.system("~/meme/bin/fimo --thresh "+ p_value +" "+meme+" "+fna_filename )

        #save in file
        #sequences and dictionary
        user_id = request.form['user_id']
        step_filename = 'app/users_workflow/'+str(user_id)+'_step3.txt'
        with open(step_filename, 'w') as g:
            json.dump(res,g)
        #save sequences
        with open('app/users_workflow/'+user_id+'_step3_dictionary.txt', 'w') as out:  
            json.dump(dictionary_snp_allele, out)
        #save to database
        file_name = str(user_id)+'_step3.txt'
        workflow = Workflow()
        user_work = Workflow.query.filter_by(user_id_user=user_id).first()
        user_work.step3 = 'app/users_workflow/'+file_name

        #delete the other work
        user_work.step4 = None
        user_work.step5 = None

        if os.path.exists("app/users_workflow/"+user_id+"_step4.txt"):
            os.remove("app/users_workflow/"+user_id+"_step4.txt")
        if os.path.exists("app/users_workflow/"+user_id+"_step5.txt"):
            os.remove("app/users_workflow/"+user_id+"_step5.txt")
        db.session.commit()

        user_email = request.form['user_email']

        try:
            server = smtplib.SMTP_SSL('smtp.gmail.com',465)
            server.login("regulomix.temp@gmail.com","regulomix123")
            subject = "Regulomix: Step3"
            body = "Step3 work is done, login to regulomix to continue."
            message = "Subject:{}\n\n{}".format(subject, body)
            server.sendmail("regulomix.temp@gmail.com",user_email,message)
            server.quit()
            print("Email sent")
        except:
            server.quit()
            print("Email failed to send")
    
    return jsonify(res,dictionary_snp_allele)

@blueprint.route('/dif_tf',methods=['GET','POST'])
@cross_origin(origin='*')
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
                count += 1
        
        #starting position of relative positions in seq_name
        start_relative_list = (count*2) + 2
  
        relative_list = seq_name.split('|')[start_relative_list:]
  
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

        #filter by id
        #new_df = df_fimo_output.drop(df_fimo_output.index[index_list],axis=0)
        idx_list_selected = [ i for i in range(len(df_fimo_output)) if not i in index_list ]

        new_df = df_fimo_output.iloc[idx_list_selected]

        return new_df

    def filter_dataframe(table_output,dictionary_snp_allele, log=False):

        new_df = filter_range(table_output)

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
                    condition_csta1 = seq_start == range(int(current_start-31),int(current_start))
                    condition_csto = current_stop == seq_stop
                    condition_csto1 = seq_stop == range(int(current_stop),int(current_stop + 31))

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
                #transfer count to total ?
                total = count
                multi_snps = 1
                count_al = 0
                #iterate through dictionary of snp and alleles
                for snp_name in snp_name_list:
                    for allele in dictionary_snp_allele[snp_name]:
                        count_al +=1
                    multi_snps *= count_al
                    count_al = 0
                # for snp_name,snp_dict in zip(snp_name_list,dictionary_snp_allele):
                #     #if snp in list of combinations is the same in the dictionary 
                #     if snp_dict == snp_name:
                #         #loop through alleles and count it 
                #         for allele in dictionary_snp_allele[snp_dict]:
                #             print ("SNP DICT: ",snp_dict,"--",allele)
                #             count_al +=1
                #         #multiply the amount of alleles to find the amount of combinations
                #         multi_snps *= count_al
                #         count_al = 0
        
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

                    if condition_m and condition_sn2 and condition_csta and condition_csto:
                        # print(motif)
                        # print(count_comb)
                        count_comb+=1
                        #new data frame droping all occurences of specified TF-snps pair
                        #new_df = new_df.drop(new_df[(new_df.motif_alt_id == motif) & (new_df.sequence_name.str.contains(snp_list))].index)
                #print("COMBINATIONS",count_comb)
                if multi_snps <= count_comb:
                    # print("COMBINATIONS",count_comb)
                    # print("MOTIF",current_motif)
                    #only the first snp
                    new_df = new_df.drop(new_df[ (new_df.motif_alt_id == current_motif) & (new_df.sequence_name.str.contains(snp_name_list[0]))].index)
          
        #display (new_df.sort_values(['motif_alt_id','start']))
        print ("DataFrame final size: ",len(new_df))
        return new_df
    
    if request.method == 'POST':
        #insert in a variable fimo path for fimo results
        fimo_res = request.form['fimo_tsv']
        alleles_dict = json.loads(request.form['al_dict'])

        print(alleles_dict)
        df_fimo_output = pd.read_csv(fimo_res, sep='\t') 
        print("DEBUG EMPTY FIMO DATAFRAME")
        print(df_fimo_output)
        print(df_fimo_output.columns)
        if 'start' not in df_fimo_output.columns:
            column_names = ['Motif ID', 'Motif Alt ID', 'Sequence Name', 'Strand', 'Start','End','p-value','q-value','Matched Sequence']
            mock_df = pd.DataFrame(columns = column_names)
            return jsonify(mock_df.to_dict(orient='records'))

        f_dataframe = filter_dataframe(fimo_res,alleles_dict,True)

        f_dataframe = f_dataframe.drop_duplicates(keep="last")

        # filter dataframe for not related snps in combinations
        for i, line in f_dataframe.iterrows():
            #variables of each line
            motif = line['motif_alt_id']
            start = line['start']
            stop = line['stop']
            sequence = line['sequence_name']
            sequence = sequence.split('|')
            if sequence[0] == "sequence_combinations":
                snps_in_each_line = len( list(filter(lambda x: (x.startswith('rs')), sequence)) )
                snps_in_line = list(filter(lambda x: (x.startswith('rs')), sequence))
                # get all the positions
                posits_start = 2*snps_in_each_line + 2 
                positions = sequence[posits_start:]
                chrom = sequence[1+snps_in_each_line:2+snps_in_each_line][0]
                start_alleles = posits_start-snps_in_each_line
                stop_alleles = start_alleles + snps_in_each_line
                alleles = sequence[posits_start-snps_in_each_line:stop_alleles]
                snps_pos_in = []
                pos_elem = []
                # iterate through the position and identify the inbetweener
                for idx,elem in enumerate(positions):
                    if int(elem) >= start and int(elem) <= stop:
                        pos_elem.append(elem)
                        snps_pos_in.append(idx)
                if (len(snps_pos_in) == 1):
                    # change the combinations to wt or variation
                    relevant_snp = snps_in_line[snps_pos_in[0]]
                    relevant_allele = alleles[snps_pos_in[0]]
                    count=0
                    for allele in alleles_dict[relevant_snp]:
                        if allele == relevant_allele:
                            if count == 0:
                                seq_name = "sequence_wild_type"
                            else:
                                seq_name = "sequence_variation"
                        count+=1
                    new_sequence_name = [seq_name, relevant_snp,chrom,relevant_allele,pos_elem[0]]
                    new_sequence_name = str.join('|',new_sequence_name)
                    #edit line
                    f_dataframe.loc[i, 'sequence_name'] = new_sequence_name
        f_dataframe = f_dataframe.drop_duplicates(keep="first")
        dataframe_out = f_dataframe.to_dict(orient='records')
        print("DATA FRAME DICT")
        # print(dataframe_out)
    #save file
    user_id = request.form['user_id']
    file_name = str(user_id)+'_step4.txt'
    with open('app/users_workflow/'+user_id+'_step4.txt', 'w') as out:  
            json.dump(dataframe_out, out)
    workflow = Workflow()
    user_work = Workflow.query.filter_by(user_id_user=user_id).first()
    user_work.step4 = 'app/users_workflow/'+file_name

    #Delete step 5
    user_work.step5 = None

    if os.path.exists("app/users_workflow/"+user_id+"_step5.txt"):
        os.remove("app/users_workflow/"+user_id+"_step5.txt")
    db.session.commit()

    user_email = request.form['user_email']

    try:
        server = smtplib.SMTP_SSL('smtp.gmail.com',465)
        server.login("regulomix.temp@gmail.com","regulomix123")
        subject = "Regulomix: Step4"
        body = "Step4 work is done, login to regulomix to continue."
        message = "Subject:{}\n\n{}".format(subject, body)
        server.sendmail("regulomix.temp@gmail.com",user_email,message)
        server.quit()
        print("Email sent")
    except:
        server.quit()
        print("Email failed to send")
    #save in database
    return jsonify(dataframe_out)

@blueprint.route('/epi',methods=['GET','POST'])
@cross_origin(origin='*')
@login_required
def epi():
    if request.method == 'POST':

        #get snp list information from the datatable
        snp_list_rows = request.form['snp_list']
        #transform it in a json file(dictionary?)
        snp_list_rows = json.loads(snp_list_rows)

        print("EPI")
        print(snp_list_rows)

        tissue_list = request.form['tissue_field_epi'].split("|")

        dataframe_epi = epi_function(snp_list_rows,tissue_list)
        dataframe_new = dataframe_epi.to_dict(orient='records')
        #save to database
        user_id = request.form['user_id']
        with open('app/users_workflow/'+user_id+'_step5.txt', 'w') as out:  
            json.dump(dataframe_new, out)
        workflow = Workflow()
        user_work = Workflow.query.filter_by(user_id_user=user_id).first()
        file_name = str(user_id)+'_step5.txt'
        user_work.step5 = 'app/users_workflow/'+file_name
        db.session.commit()
        print(dataframe_new)

        user_email = request.form['user_email']
        
        try:
            server = smtplib.SMTP_SSL('smtp.gmail.com',465)
            server.login("regulomix.temp@gmail.com","regulomix123")
            subject = "Regulomix: Step5"
            body = "Step5 work is done, login to regulomix to continue."
            message = "Subject:{}\n\n{}".format(subject, body)
            server.sendmail("regulomix.temp@gmail.com",user_email,message)
            server.quit()
            print("Email sent")
        except:
            server.quit()
            print("Email failed to send")

    return jsonify(dataframe_new)

@blueprint.route('/delete_old',methods=['GET','POST'])
@cross_origin(origin='*')
@login_required
def delete_old():
    workflow = Workflow.query.all()
    today = datetime.now()
    #loop entries in db
    for work in workflow:
        #delete if it already expired
        if (work.expire < today):
            #delete files and entries
            work.step1 = None
            work.step2 = None
            work.step3 = None
            work.step4 = None
            work.step5 = None

            work_id = str(work.user_id_user)
            #Delete from directory
            if os.path.exists("app/users_workflow/"+work_id+"_step1.txt"):
                os.remove("app/users_workflow/"+work_id+"_step1.txt")
            if os.path.exists("app/users_workflow/"+work_id+"_step2.txt"):
                os.remove("app/users_workflow/"+work_id+"_step2.txt")
            if os.path.exists("app/users_workflow/"+work_id+"_step3.txt"):
                os.remove("app/users_workflow/"+work_id+"_step3.txt")
            if os.path.exists("app/users_workflow/"+work_id+"_step3_dictionary.txt"):
                os.remove("app/users_workflow/"+work_id+"_step3_dictionary.txt")
            if os.path.exists("app/users_workflow/"+work_id+"_step4.txt"):
                os.remove("app/users_workflow/"+work_id+"_step4.txt")
            if os.path.exists("app/users_workflow/"+work_id+"_step5.txt"):
                os.remove("app/users_workflow/"+work_id+"_step5.txt")
            # commit the record the database
            db.session.commit()
            # Workflow.query.filter(Workflow.id_workflow == work.id_workflow).delete()
    return "Done delete function"

@cross_origin(origin='*')
@login_required
def check_status(task_id):
    pass

@blueprint.route('/data_retrive',methods=['GET','POST'])
@cross_origin(origin='*')
@login_required
def data_retriever():
    id_user = request.form['data_user']
    find_user_work = Workflow.query.filter_by(user_id_user=id_user).first()
    if not find_user_work == None and not find_user_work.step1 == None:
        with open('app/users_workflow/'+str(find_user_work.user_id_user)+'_step1.txt') as json_file:  
            data = json.load(json_file)
            return jsonify(data)
    else:
        return 'No File'

@blueprint.route('/data_retrive1',methods=['GET','POST'])
@cross_origin(origin='*')
@login_required
def data_retriever1():
    id_user = request.form['data_user']
    find_user_work = Workflow.query.filter_by(user_id_user=id_user).first()
    if not find_user_work == None and not find_user_work.step2 == None:
        with open('app/users_workflow/'+str(find_user_work.user_id_user)+'_step2.txt') as json_file:  
            data = json.load(json_file)
            return jsonify(data)
    else:
        return 'No File'

@blueprint.route('/data_retrive2',methods=['GET','POST'])
@cross_origin(origin='*')
@login_required
def data_retriever2():
    id_user = request.form['data_user']
    find_user_work = Workflow.query.filter_by(user_id_user=id_user).first()
    if not find_user_work.step3 == None:
        file1 = 'app/users_workflow/'+str(find_user_work.user_id_user)+'_step3.txt'
        file2 = 'app/users_workflow/'+str(find_user_work.user_id_user)+'_step3_dictionary.txt'
        with open(file1) as json_file, open(file2) as b:  
            data = json.load(json_file)
            dictio = json.load(b)
            return jsonify(data,dictio)
    else:
        return 'No File'

@blueprint.route('/data_retrive3',methods=['GET','POST'])
@cross_origin(origin='*')
@login_required
def data_retriever3():
    id_user = request.form['data_user']
    find_user_work = Workflow.query.filter_by(user_id_user=id_user).first()
    if not find_user_work.step4 == None:
        file1 = 'app/users_workflow/'+str(find_user_work.user_id_user)+'_step4.txt'
        with open(file1) as json_file:  
            data = json.load(json_file)
            return jsonify(data)
    else:
        return 'No File'

@blueprint.route('/data_retrive4',methods=['GET','POST'])
@cross_origin(origin='*')
@login_required
def data_retriever4():
    id_user = request.form['data_user']
    find_user_work = Workflow.query.filter_by(user_id_user=id_user).first()
    if not find_user_work.step5 == None:
        file1 = 'app/users_workflow/'+str(find_user_work.user_id_user)+'_step5.txt'
        with open(file1) as json_file:  
            data = json.load(json_file)
            return jsonify(data)
    else:
        return 'No File'

@blueprint.route('/next_step1',methods=['GET','POST'])
@cross_origin(origin='*')
@login_required
def next_step1():
    if request.method == 'POST':
        #TODO better data retrieve function
        #commit to database the file changed
        user_id = request.form['user']
        data_table = request.form['data']
        data_real = json.loads(data_table)
        #save to file
        file_name = str(user_id)+'_step1.txt'
        with open('app/users_workflow/'+file_name, 'w') as outfile:  
            json.dump(data_real, outfile)
        # set variables
        today = datetime.now()
        expire_date = today + timedelta(days=8)
        workflow = Workflow()
        user_work = Workflow.query.filter_by(user_id_user=user_id).first()
        if(user_work == None):
            # print("ENTROU NA PRIMEIRA ADD")
            workflow.user_id_user = user_id
            workflow.created = datetime.now()
            workflow.step1 = 'app/users_workflow/'+file_name
            workflow.expire = expire_date
            db.session.add(workflow)
        else:
            #update workflow
            if (user_work.step1 == None):
                print("step1 -- empty for user 8")
                user_work.created = datetime.now()
                user_work.expire = expire_date
                user_work.step1 = 'app/users_workflow/'+file_name
                db.session.commit()
        # commit the record the database
        db.session.commit()
        #sqlalchemy insert
        return 'DONE'

@blueprint.route('/next_step2',methods=['GET','POST'])
@cross_origin(origin='*')
@login_required
def next_step2():
    if request.method == 'POST':
        #commit tables one and two
        data_table1 = request.form['data']
        data_table2 = request.form['data2']
        dados1 = json.loads(data_table1)
        dados2 = json.loads(data_table2)
        user_id = request.form['user']
        #update both files
        file_name1 = 'app/users_workflow/'+user_id+'_step1.txt'
        file_name2 = 'app/users_workflow/'+user_id+'_step2.txt'
        with open(file_name1, 'w') as outfile:  
            json.dump(dados1, outfile)
        with open(file_name2, 'w') as f:  
            json.dump(dados2, f)

        return 'DONE'

@blueprint.route('/data_retrive_thread',methods=['GET','POST'])
@cross_origin(origin='*')
@login_required
def data_retriever_thread():
    
    id_user = request.form['data_user']
    # see if user_thread has ended
    # get user thread
    check_thread = False
    print(threading.enumerate())
    for thread in threading.enumerate():
        print(thread.getName())
        if(thread.getName() == "user_"+id_user):
            check_thread = True
    
    file_exists = os.path.exists("app/users_workflow/"+id_user+"_step1.txt")
    if((check_thread == False) and (file_exists == False)):
        return 'Thread Error'

    find_user_work = Workflow.query.filter_by(user_id_user=id_user).first()
    if not find_user_work.step1 == None:
        with open('app/users_workflow/'+str(find_user_work.user_id_user)+'_step1.txt') as json_file:  
            data = json.load(json_file)
            return jsonify(data)
    else:
        return 'No File'

@blueprint.route('/data_retrive2_thread',methods=['GET','POST'])
@cross_origin(origin='*')
@login_required
def data_retriever_thread2():
    
    id_user = request.form['data_user']
    # see if user_thread has ended
    # get user thread
    check_thread = False
    print(threading.enumerate())
    for thread in threading.enumerate():
        print(thread.getName())
        if(thread.getName() == "step2_user_"+id_user):
            check_thread = True
    
    file_exists = os.path.exists("app/users_workflow/"+id_user+"_step2.txt")
    if((check_thread == False) and (file_exists == False)):
        return 'Thread Error'

    find_user_work = Workflow.query.filter_by(user_id_user=id_user).first()
    if not find_user_work.step2 == None:
        with open('app/users_workflow/'+str(find_user_work.user_id_user)+'_step2.txt') as json_file:  
            data = json.load(json_file)
            return jsonify(data)
    else:
        return 'No File'

@blueprint.route('/uploader',methods=['GET','POST'])
@cross_origin(origin='*')
@login_required
def uploader():
    if request.method == 'POST':
        user_id = request.form['user_id']
        user_email = request.form['user_email']
        f = request.files['file']
        file_path = "app/home/drscan.xlsx"
        f.save(file_path)

        user_work = Workflow.query.filter_by(user_id_user=user_id).first()
        if (not user_work is None):
            user_work.step2 = None
            user_work.step3 = None
            user_work.step4 = None
            user_work.step5 = None
            db.session.commit()
        #Delete from directory
        if os.path.exists("app/users_workflow/"+user_id+"_step2.txt"):
            os.remove("app/users_workflow/"+user_id+"_step2.txt")
        if os.path.exists("app/users_workflow/"+user_id+"_step3.txt"):
            os.remove("app/users_workflow/"+user_id+"_step3.txt")
        if os.path.exists("app/users_workflow/"+user_id+"_step3_dictionary.txt"):
            os.remove("app/users_workflow/"+user_id+"_step3_dictionary.txt")
        if os.path.exists("app/users_workflow/"+user_id+"_step4.txt"):
            os.remove("app/users_workflow/"+user_id+"_step4.txt")
        if os.path.exists("app/users_workflow/"+user_id+"_step5.txt"):
            os.remove("app/users_workflow/"+user_id+"_step5.txt")

        # @copy_current_request_context
        # def execute_snp_task(user_id, user_email, file_path):
        #     funcao_teste_snp(user_id, user_email, file_path)
            # snps_from_file(user_id, user_email, file_path)

        #using celery
        task_res = funcao_teste_snp.apply_async(args=[user_id, user_email, file_path])
        
        # thread_name = "user_"+user_id
        # try:
        #     threading.Thread(name = thread_name,target=execute_snp_task,
        #                 args=(user_id, user_email, file_path)).start()
        # except:
        #     print("Excecao Thread")
        # funcao_teste_snp(request)
        return jsonify({}), 202, {'Location': url_for('home_blueprint.taskstatus1',task_id=task_res.id)}
        # return jsonify(resultado = {}, status=202, {'Location': url_for('taskstatus1',task_id=task_res.id)})

def snps_from_file(user_id, user_email, file_path):
    print("Threaded DBSnp File")
    user_work = Workflow.query.filter_by(user_id_user=user_id).first()
    if os.path.exists("app/users_workflow/"+user_id+"_step1.txt"):
        os.remove("app/users_workflow/"+user_id+"_step1.txt")

    # # Removing from database
    user_work.step1 = None
    db.session.commit()
    
    # save xlsx file
    # xlsx file into pandas dataframe
    input = pd.read_excel(file_path)
    # drop duplicates
    input.drop_duplicates(subset=['SNPS'], keep='first',inplace=True)
    # create snp list
    snp_list = input['SNPS']

    #saving json file
    teste = get_all_snps_file(snp_list)
    file_name = str(user_id)+'_step1.txt'
    teste.to_json('app/users_workflow/'+file_name)
    # print(file_name)
    # with open('app/users_workflow/'+file_name, 'w') as outfile:  
    #     json.dump(teste, outfile)
         
    #insert into db (timestamp,id_user,step 1)
    today = datetime.now()
    expire_date = today + timedelta(days=8)
    workflow = Workflow()
 
    #select workflow to see if file slready exists
    user_work = Workflow.query.filter_by(user_id_user=user_id).first()
    if(user_work == None):
        # print("ENTROU NA PRIMEIRA ADD")
        workflow.user_id_user = user_id
        workflow.created = datetime.now()
        workflow.step1 = 'app/users_workflow/'+file_name
        workflow.expire = expire_date
        db.session.add(workflow)
    else:
        #update workflow
        user_work.created = datetime.now()
        user_work.step1 = 'app/users_workflow/'+file_name
        user_work.expire = expire_date
         
    db.session.commit()
         
    #delete previous workflow
    #TODO second delete data from the database
    user_work.step2 = None
    user_work.step3 = None
    user_work.step4 = None
    user_work.step5 = None
    #Delete from directory
    if os.path.exists("app/users_workflow/"+user_id+"_step2.txt"):
        os.remove("app/users_workflow/"+user_id+"_step2.txt")
    if os.path.exists("app/users_workflow/"+user_id+"_step3.txt"):
        os.remove("app/users_workflow/"+user_id+"_step3.txt")
    if os.path.exists("app/users_workflow/"+user_id+"_step3_dictionary.txt"):
        os.remove("app/users_workflow/"+user_id+"_step3_dictionary.txt")
    if os.path.exists("app/users_workflow/"+user_id+"_step4.txt"):
        os.remove("app/users_workflow/"+user_id+"_step4.txt")
    if os.path.exists("app/users_workflow/"+user_id+"_step5.txt"):
        os.remove("app/users_workflow/"+user_id+"_step5.txt")
    # commit the record the database
    db.session.commit()
    #try to send email
    try:
        server = smtplib.SMTP_SSL('smtp.gmail.com',465)
        server.login("regulomix.temp@gmail.com","regulomix123")
        subject = "Regulomix: Step1"
        body = "Step1 work is done, login to regulomix to continue."
        message = "Subject:{}\n\n{}".format(subject, body)
        server.sendmail("regulomix.temp@gmail.com",user_email,message)
        server.quit()
        print("Email sent")
    except:
        server.quit()
        print("Email failed to send")

def get_all_snps_file(snps):
    print("GET ALL SNPS FILE")
    vcf = './GRCh37_vcf/common_all_20180423.vcf'
    df = allel.vcf_to_dataframe(vcf)
    print("DF DONE")
    map(str.strip, snps)
    df = df.loc[df['ID'].isin(snps)]
    df["ALT_AL"] = df["ALT_1"].fillna('') + df["ALT_2"].fillna('') + df["ALT_3"].fillna('')
    # Strip values
    df['ALT_AL'] = df['ALT_AL'].apply(lambda x: x.strip())
    # Split with pipebar |
    df['ALT_AL'] = df['ALT_AL'].apply(lambda x: '|'.join(x))
    df = df[df['REF'].map(len) < 2]
    df = df[df['ALT_1'].map(len) < 2]
    df.drop(['ALT_1','ALT_2', 'ALT_3', 'QUAL','FILTER_PASS'], axis=1, inplace=True)
    df = df.reset_index(drop=True)
    # get the ones from the file
    print("PRINT DF")
    print(df)
    rest = Diff(snps, df['ID'].tolist())
    snps_no_rs = map(lambda x: x[2:], rest)
    result = process(snps_no_rs)
    for i in result:
        chromo = i[0][0]['chrom']
        pos = i[0][0]['location']
        ref = i[0][0]['allele_wt']
        alt = i[2]
        identif = i[1]
        df2 = pd.DataFrame([(chromo, pos, identif,ref,alt)],columns=['CHROM','POS','ID','REF','ALT_AL'])
        df = df.append(df2,ignore_index = True)
    return df

def Diff(li1, li2): 
    return (list(set(li1) - set(li2)))

#Will not be used
@celery.task(bind=True)
def funcao_teste_snp(self, user_id, user_email, file_path):
    print("Funcao Teste SNP Thread")
    # delete archive from users_workflow and db if exists
    # # Removing from directory
    user_work = Workflow.query.filter_by(user_id_user=user_id).first()
    if os.path.exists("app/users_workflow/"+user_id+"_step1.txt"):
        os.remove("app/users_workflow/"+user_id+"_step1.txt")

    # # Removing from database
    if (not user_work is None):
        user_work.step1 = None
        db.session.commit()
    
    # save xlsx file
    # xlsx file into pandas dataframe
    input = pd.read_excel(file_path)
    # drop duplicates
    input.drop_duplicates(subset=['SNPS'], keep='first',inplace=True)
    # create snp list
    snp_list = input['SNPS']
    # only snp ids/ without rs
    snp_ids = []
    for snp in snp_list:
        snp_ids.append(snp.strip()[2:])
 
    #saving json file
    #maybe take out the process function
    # teste = process(snp_ids)
    list_all_snps = []
    for idx, snp in enumerate(snp_ids):

        minor_allele = []
        a = get_snp_info(snp)
        #verify if its an empty dataframe
        if a.empty:
            # print("NOT SNV")
            # print("rs"+snp)
            continue
        #dictionary with information based on genome version
        sample_dict = a[a["gnenome_versions"].str.contains("GRCh37")]['snp_info_dict'].values[0]
        #if sample dict is bigger than 1(more than one allele in variation)
        if len(sample_dict) > 1:
            for dic in sample_dict:
                minor_allele.append(dic['allele_v'])
            #Stringfy minor allele
            join_symbol = "|"
            str_minor = join_symbol.join(minor_allele)
        else:
            str_minor = sample_dict[0]['allele_v']
        snp_id_form = 'rs'+snp
        snp_in = [sample_dict,snp_id_form,str_minor]
        list_all_snps.append(snp_in)

        self.update_state(state='PROGRESS',
                          meta={'current': idx+1, 
                                'total': len(snp_ids),
                                'status': "SNP {0}:{1} is done being processed from a total of {2}".format(str(idx+1),snp,str(len(snp_ids)))})
    
    self.update_state(state='DONE',meta={'current': idx+1, 
                                'total': len(snp_ids),
                                'status': "All SNPs are done. Results should appear shortly"})
    # return list_all_snps
    
    file_name = str(user_id)+'_step1.txt'
    # print(file_name)
    with open('app/users_workflow/'+file_name, 'w') as outfile:  
        json.dump(list_all_snps, outfile)
         
    #insert into db (timestamp,id_user,step 1)
    today = datetime.now()
    expire_date = today + timedelta(days=8)
    workflow = Workflow()
 
    #select workflow to see if file slready exists
    user_work = Workflow.query.filter_by(user_id_user=user_id).first()
    if(user_work == None):
        # print("ENTROU NA PRIMEIRA ADD")
        workflow.user_id_user = user_id
        workflow.created = datetime.now()
        workflow.step1 = 'app/users_workflow/'+file_name
        workflow.expire = expire_date
        db.session.add(workflow)
    else:
        #update workflow
        user_work.created = datetime.now()
        user_work.step1 = 'app/users_workflow/'+file_name
        user_work.expire = expire_date
         
    db.session.commit()
    #Added recently
    user_work = Workflow.query.filter_by(user_id_user=user_id).first()
    #delete previous workflow
    #TODO second delete data from the database
    user_work.step2 = None
    user_work.step3 = None
    user_work.step4 = None
    user_work.step5 = None
    #Delete from directory
    if os.path.exists("app/users_workflow/"+user_id+"_step2.txt"):
        os.remove("app/users_workflow/"+user_id+"_step2.txt")
    if os.path.exists("app/users_workflow/"+user_id+"_step3.txt"):
        os.remove("app/users_workflow/"+user_id+"_step3.txt")
    if os.path.exists("app/users_workflow/"+user_id+"_step3_dictionary.txt"):
        os.remove("app/users_workflow/"+user_id+"_step3_dictionary.txt")
    if os.path.exists("app/users_workflow/"+user_id+"_step4.txt"):
        os.remove("app/users_workflow/"+user_id+"_step4.txt")
    if os.path.exists("app/users_workflow/"+user_id+"_step5.txt"):
        os.remove("app/users_workflow/"+user_id+"_step5.txt")
    # commit the record the database
    db.session.commit()
    #try to send email
    try:
        server = smtplib.SMTP_SSL('smtp.gmail.com',465)
        server.login("regulomix.temp@gmail.com","regulomix123")
        subject = "Regulomix: Step1"
        body = "Step1 work is done, login to regulomix to continue."
        message = "Subject:{}\n\n{}".format(subject, body)
        server.sendmail("regulomix.temp@gmail.com",user_email,message)
        server.quit()
        print("Email sent")
    except:
        server.quit()
        print("Email failed to send")
    
    # Pass all parameters in result
    snp_info_list = []
    for info in list_all_snps:
        print(info)
        dict_snp = {}
        dict_snp['snp'] = info[1]
        dict_snp['al_w'] = info[0][0]['allele_wt']
        dict_snp['al_a'] = info[0][0]['allele_v']
        dict_snp['local'] = info[0][0]['location']
        dict_snp['chrom'] = info[0][0]['chrom']
        snp_info_list.append(dict_snp)
    return {'current': 100, 'total': 100, 'status': 'All completed','result': snp_info_list}
    #sqlalchemy inserts
    # return jsonify(teste)


@blueprint.route('/taskstatus1/<task_id>',methods=['GET','POST'])
@cross_origin(origin='*')
@login_required
def taskstatus1(task_id):
    task = funcao_teste_snp.AsyncResult(task_id)
    if task.state == 'PENDING':
        response = {
            'state': task.state,
            'current': 0,
            'total': 1,
            'status': 'Pending...'
        }
    elif task.state != 'FAILURE':
        response = {
            'state': task.state,
            'current': task.info.get('current', 0),
            'total': task.info.get('total', 1),
            'status': task.info.get('status', '')
        }
        if 'result' in task.info:
            response['result'] = task.info['result']
    else:
        # something went wrong in the background job
        response = {
            'state': task.state,
            'current': 1,
            'total': 1,
            'status': str(task.info),  # this is the exception raised
        }
    return jsonify(response)

#Asynchronous function
async def async_process_snp_upload(snp_ids):
    data_frame_snps = process(snp_ids)
    return data_frame_snps

@blueprint.route('/upload_matrix',methods=['GET','POST'])
@cross_origin(origin='*')
@login_required
def upload_matrix():
    if request.method == 'POST':
        # print("REQUEST")
        # print(request.files) 
        f = request.files['file_matrix']
        file_path = "app/home/meme.meme"
        #save xlsx file
        f.save(file_path)
        
        return "OK"

def process(snp_ids):
    list_all_snps = []
    for snp in snp_ids:

        minor_allele = []
        a = get_snp_info(snp)
        #verify if its an empty dataframe
        if a.empty:
            # print("NOT SNV")
            # print("rs"+snp)
            continue
        #dictionary with information based on genome version
        sample_dict = a[a["gnenome_versions"].str.contains("GRCh37")]['snp_info_dict'].values[0]
        #if sample dict is bigger than 1(more than one allele in variation)
        if len(sample_dict) > 1:
            for dic in sample_dict:
                minor_allele.append(dic['allele_v'])
            #Stringfy minor allele
            join_symbol = "|"
            str_minor = join_symbol.join(minor_allele)
        else:
            str_minor = sample_dict[0]['allele_v']
        snp_id_form = 'rs'+snp
        snp_in = [sample_dict,snp_id_form,str_minor]
        list_all_snps.append(snp_in)
    return list_all_snps
        #return dataframe

@blueprint.route('index2/get_table',methods=['GET','POST'])
@cross_origin(origin='*')
@login_required
def get_info_dictionary():
    ##dictionary roadmap epigenomics
    roadmap = {
    
        1:['E001','ES-I3 Cells'],
        2:['E002','ES-WA7 Cells'],
        3:['E003','H1 Cells'],
        4:['E004','H1 BMP4 Derived Mesendoderm Cultured Cells'],
        5:['E005','H1 BMP4 Derived Trophoblast Cultured Cells'],
        6:['E006','H1 Derived Mesenchymal Stem Cells'],
        7:['E007','H1 Derived Neuronal Progenitor Cultured Cells'],
        8:['E008','H9 Cells'],
        9:['E009','H9 Derived Neuronal Progenitor Cultured Cells'],
        10:['E010','H9 Derived Neuron Cultured Cells'],
        11:['E011','hESC Derived CD184+ Endoderm Cultured Cells'],
        12:['E012','hESC Derived CD56+ Ectoderm Cultured Cells'],
        13:['E013','hESC Derived CD56+ Mesoderm Cultured Cells'],
        14:['E014','HUES48 Cells'],
        15:['E015','HUES6 Cells'],
        16:['E016','HUES64 Cells'],
        17:['E017','IMR90 fetal lung fibroblasts Cell Line'],
        18:['E018','iPS-15b Cells'],
        19:['E019','iPS-18 Cells'],
        20:['E020','iPS-20b Cells'],
        21:['E021','iPS DF 6.9 Cells'],
        22:['E022','iPS DF 19.11 Cells'],
        23:['E023','Mesenchymal Stem Cell Derived Adipocyte Cultured Cells'],
        24:['E024','ES-UCSF4 Cells'],
        25:['E025','Adipose Derived Mesenchymal Stem Cell Cultured Cells'],
        26:['E026','Bone Marrow Derived Cultured Mesenchymal Stem Cells'],
        27:['E027','Breast Myoepithelial Primary Cells'],
        28:['E028','Breast variant Human Mammary Epithelial Cells (vHMEC)'],
        29:['E029','Primary monocytes from peripheral blood'],
        30:['E030','Primary neutrophils from peripheral blood'],
        31:['E031','Primary B cells from cord blood'],
        32:['E032','Primary B cells from peripheral blood'],
        33:['E033','Primary T cells from cord blood'],
        34:['E034','Primary T cells from peripheral blood'],
        35:['E035','Primary hematopoietic stem cells'],
        36:['E036','Primary hematopoietic stem cells short term culture'],
        37:['E037','Primary T helper memory cells from peripheral blood 2'],
        38:['E038','Primary T helper naive cells from peripheral blood'],
        39:['E039','Primary T helper naive cells from peripheral blood'],
        40:['E040','Primary T helper memory cells from peripheral blood 1'],
        41:['E041','Primary T helper cells PMA-I stimulated'],
        42:['E042','Primary T helper 17 cells PMA-I stimulated'],
        43:['E043','Primary T helper cells from peripheral blood'],
        44:['E044','Primary T regulatory cells from peripheral blood'],
        45:['E045','Primary T cells effector/memory enriched from peripheral blood'],
        46:['E046','Primary Natural Killer cells from peripheral blood'],
        47:['E047','Primary T CD8+ naive cells from peripheral blood'],
        48:['E048','Primary T CD8+ memory cells from peripheral blood'],
        49:['E049','Mesenchymal Stem Cell Derived Chondrocyte Cultured Cells'],
        50:['E050','Primary hematopoietic stem cells G-CSF-mobilized Female'],
        51:['E051','Primary hematopoietic stem cells G-CSF-mobilized Male'],
        52:['E052','Muscle Satellite Cultured Cells'],
        53:['E053','Cortex derived primary cultured neurospheres'],
        54:['E054','Ganglion Eminence derived primary cultured neurospheres'],
        55:['E055','Foreskin Fibroblast Primary Cells skin01'],
        56:['E056','Foreskin Fibroblast Primary Cells skin02'],
        57:['E057','Foreskin Keratinocyte Primary Cells skin02'],
        58:['E058','Foreskin Keratinocyte Primary Cells skin03'],
        59:['E059','Foreskin Melanocyte Primary Cells skin01'],
        60:['E061','Foreskin Melanocyte Primary Cells skin03'],
        61:['E062','Primary mononuclear cells from peripheral blood'],
        62:['E063','Adipose Nuclei'],
        63:['E065','Aorta'],
        64:['E066','Liver'],
        65:['E067','Brain Angular Gyrus'],
        66:['E068','Brain Anterior Caudate'],
        67:['E069','Brain Cingulate Gyrus'],
        68:['E070','Brain Germinal Matrix'],
        69:['E071','Brain Hippocampus Middle'],
        70:['E072','Brain Inferior Temporal Lobe'],
        71:['E073','Brain_Dorsolateral_Prefrontal_Cortex'],
        72:['E074','Brain Substantia Nigra'],
        73:['E075','Colonic Mucosa'],
        74:['E076','Colon Smooth Muscle'],
        75:['E077','Duodenum Mucosa'],
        76:['E078','Duodenum Smooth Muscle'],
        77:['E079','Esophagus'],
        78:['E080','Fetal Adrenal Gland'],
        79:['E081','Fetal Brain Male'],
        80:['E082','Fetal Brain Female'],
        81:['E083','Fetal Heart'],
        82:['E084','Fetal Intestine Large'],
        83:['E085','Fetal Intestine Small'],
        84:['E086','Fetal Kidney'],
        85:['E087','Pancreatic Islets'],
        86:['E088','Fetal Lung'],
        87:['E089','Fetal Muscle Trunk'],
        88:['E090','Fetal Muscle Leg'],
        89:['E091','Placenta'],
        90:['E092','Fetal Stomach'],
        91:['E093','Fetal Thymus'],
        92:['E094','Gastric'],
        93:['E095','Left Ventricle'],
        94:['E096','Lung'],
        95:['E097','Ovary'],
        96:['E098','Pancreas'],
        97:['E099','Placenta Amnion'],
        98:['E100','Psoas Muscle'],
        99:['E101','Rectal Mucosa Donor 29'],
        100:['E102','Rectal Mucosa Donor 31'],
        101:['E103','Rectal Smooth Muscle'],
        102:['E104','Right Atrium'],
        103:['E105','Right Ventricle'],
        104:['E106','Sigmoid Colon'],
        105:['E107','Skeletal Muscle Male'],
        106:['E108','Skeletal Muscle Female'],
        107:['E109','Small Intestine'],
        108:['E110','Stomach Mucosa'],
        109:['E111','Stomach Smooth Muscle'],
        110:['E112','Thymus'],
        111:['E113','Spleen'],
        112:['E114','A549 EtOH 0.02pct Lung Carcinoma Cell Line'],
        113:['E115','Dnd41 TCell Leukemia Cell Line'],
        114:['E116','GM12878 Lymphoblastoid Cells'],
        115:['E117','HeLa-S3 Cervical Carcinoma Cell Line'],
        116:['E118','HepG2 Hepatocellular Carcinoma Cell Line'],
        117:['E119','HMEC Mammary Epithelial Primary Cells'],
        118:['E120','HSMM Skeletal Muscle Myoblasts Cells'],
        119:['E121','HSMM cell derived Skeletal Muscle Myotubes Cells'],
        120:['E122','HUVEC Umbilical Vein Endothelial Primary Cells'],
        121:['E123','K562 Leukemia Cells'],
        122:['E124','Monocytes-CD14+ RO01746 Primary Cells'],
        123:['E125','NH-A Astrocytes Primary Cells'],
        124:['E126','NHDF-Ad Adult Dermal Fibroblast Primary Cells'],
        125:['E127','NHEK-Epidermal Keratinocyte Primary Cells'],
        126:['E128','NHLF Lung Fibroblast Primary Cells'],
        127:['E129','Osteoblast Primary Cells'],
    }
    ##dictionary fantom 5
    fantom_5 = {
    
        1:['CNhs10608','Clontech Human Universal Reference Total RNA pool1'],
        2:['CNhs10610','SABiosciences XpressRef Human Universal Total RNA pool1'],
        3:['CNhs10612','Universal RNA - Human Normal Tissues Biochain pool1'],
        4:['CNhs10615','adipose tissue adult pool1'],
        5:['CNhs10616','bladder adult pool1'],
        6:['CNhs10617','brain adult pool1'],
        7:['CNhs10618','cervix adult pool1'],
        8:['CNhs10619','colon adult pool1'],
        9:['CNhs10620','esophagus adult pool1'],
        10:['CNhs10621','heart adult pool1'],
        11:['CNhs10622','kidney adult pool1'],
        12:['CNhs10624','liver adult pool1'],
        13:['CNhs10625','lung adult pool1'],
        14:['CNhs10626','ovary adult pool1'],
        15:['CNhs10627','placenta adult pool1'],
        16:['CNhs10628','prostate adult pool1'],
        17:['CNhs10629','skeletal muscle adult pool1'],
        18:['CNhs10630','small intestine adult pool1'],
        19:['CNhs10631','spleen adult pool1'],
        20:['CNhs10632','testis adult pool1'],
        21:['CNhs10633','thymus adult pool1'],
        22:['CNhs10634','thyroid adult pool1'],
        23:['CNhs10635','trachea adult pool1'],
        24:['CNhs10636','retina adult pool1'],
        25:['CNhs10637','temporal lobe adult pool1'],
        26:['CNhs10638','postcentral gyrus adult pool1'],
        27:['CNhs10640','pons adult pool1'],
        28:['CNhs10641','parietal lobe adult pool1'],
        29:['CNhs10642','paracentral gyrus adult pool1'],
        30:['CNhs10643','occipital pole adult pool1'],
        31:['CNhs10644','nucleus accumbens adult pool1'],
        32:['CNhs10645','medulla oblongata adult pool1'],
        33:['CNhs10646','insula adult pool1'],
        34:['CNhs10647','frontal lobe adult pool1'],
        35:['CNhs10648','dura mater adult donor1'],
        36:['CNhs10649','corpus callosum adult pool1'],
        37:['CNhs10650','thymus fetal pool1'],
        38:['CNhs10651','spleen fetal pool1'],
        39:['CNhs10652','kidney fetal pool1'],
        40:['CNhs10653','heart fetal pool1'],
        41:['CNhs10654','tonsil adult pool1'],
        42:['CNhs10722','acute myeloid leukemia (FAB M5) cell line:THP-1 (fresh)'],
        43:['CNhs10723','acute myeloid leukemia (FAB M5) cell line:THP-1 (revived)'],
        44:['CNhs10724','acute myeloid leukemia (FAB M5) cell line:THP-1 (thawed)'],
        45:['CNhs10726','lung adenocarcinoma cell line:PC-14'],
        46:['CNhs10727','chronic myelogenous leukemia cell line:KU812'],
        47:['CNhs10728','extraskeletal myxoid chondrosarcoma cell line:H-EMC-SS'],
        48:['CNhs10729','renal cell carcinoma cell line:OS-RC-2'],
        49:['CNhs10730','malignant trichilemmal cyst cell line:DJM-1'],
        50:['CNhs10731','glioma cell line:GI-1'],
        51:['CNhs10732','maxillary sinus tumor cell line:HSQ-89'],
        52:['CNhs10733','gall bladder carcinoma cell line:TGBC2TKB'],
        53:['CNhs10734','papillotubular adenocarcinoma cell line:TGBC18TKB'],
        54:['CNhs10735','transitional-cell carcinoma cell line:5637'],
        55:['CNhs10736','breast carcinoma cell line:MDA-MB-453'],
        56:['CNhs10737','colon carcinoma cell line:COLO-320'],
        57:['CNhs10738','adult T-cell leukemia cell line:ATN-1'],
        58:['CNhs10739','Burkitt\'s lymphoma cell line:DAUDI'],
        59:['CNhs10740','choriocarcinoma cell line:BeWo'],
        60:['CNhs10741','splenic lymphoma with villous lymphocytes cell line:SLVL'],
        61:['CNhs10742','astrocytoma cell line:TM-31'],
        62:['CNhs10743','epidermoid carcinoma cell line:A431'],
        63:['CNhs10744','b cell line:RPMI1788'],
        64:['CNhs10745','anaplastic carcinoma cell line:8305C'],
        65:['CNhs10746','acute lymphoblastic leukemia (T-ALL) cell line:HPB-ALL'],
        66:['CNhs10747','non T non B acute lymphoblastic leukemia (ALL) cell line:P30/OHK'],
        67:['CNhs10748','epidermoid carcinoma cell line:Ca Ski'],
        68:['CNhs10750','bile duct carcinoma cell line:HuCCT1'],
        69:['CNhs10751','giant cell carcinoma cell line:Lu99B'],
        70:['CNhs10752','oral squamous cell carcinoma cell line:Ca9-22'],
        71:['CNhs10753','signet ring carcinoma cell line:Kato III'],
        72:['CNhs10837','Endothelial Cells - Aortic donor0'],
        73:['CNhs10838','Smooth Muscle Cells - Aortic donor0'],
        74:['CNhs10839','Smooth Muscle Cells - Umbilical artery donor0'],
        75:['CNhs10842','Retinal Pigment Epithelial Cells donor0'],
        76:['CNhs10843','Urothelial cells donor0'],
        77:['CNhs10844','Mesenchymal stem cells - adipose donor0'],
        78:['CNhs10845','Mesenchymal stem cells - hepatic donor0'],
        79:['CNhs10847','Sebocyte donor1'],
        80:['CNhs10848','Fibroblast - Gingival donor4 (GFH2)'],
        81:['CNhs10850','Mesothelial Cells donor1'],
        82:['CNhs10851','Sertoli Cells donor1'],
        83:['CNhs10852','CD14+ Monocytes donor1'],
        84:['CNhs10853','CD4+ T Cells donor1'],
        85:['CNhs10854','CD8+ T Cells donor1'],
        86:['CNhs10855','Dendritic Cells - monocyte immature derived donor1 tech_rep1'],
        87:['CNhs10857','Dendritic Cells - plasmacytoid donor1'],
        88:['CNhs10858','CD14+ monocyte derived endothelial progenitor cells donor1'],
        89:['CNhs10859','Natural Killer Cells donor1'],
        90:['CNhs10860','Peripheral Blood Mononuclear Cells donor1'],
        91:['CNhs10861','Macrophage - monocyte derived donor1'],
        92:['CNhs10862','Neutrophils donor1'],
        93:['CNhs10863','Smooth Muscle Cells - Brain Vascular donor1'],
        94:['CNhs10864','Astrocyte - cerebral cortex donor1'],
        95:['CNhs10865','Endothelial Cells - Lymphatic donor1'],
        96:['CNhs10866','Fibroblast - Gingival donor1'],
        97:['CNhs10867','Fibroblast - Periodontal Ligament donor1'],
        98:['CNhs10868','Smooth Muscle Cells - Colonic donor1'],
        99:['CNhs10869','Skeletal Muscle Satellite Cells donor1'],
        100:['CNhs10870','Myoblast donor1'],
        101:['CNhs10871','Ciliary Epithelial Cells donor1'],
        102:['CNhs10872','Endothelial Cells - Umbilical vein donor1'],
        103:['CNhs10874','Fibroblast - Aortic Adventitial donor1'],
        104:['CNhs10875','Intestinal epithelial cells (polarized) donor1'],
        105:['CNhs10876','Anulus Pulposus Cell donor1'],
        106:['CNhs10877','Pancreatic stromal cells donor1'],
        107:['CNhs10878','Fibroblast - Pulmonary Artery donor1'],
        108:['CNhs10879','Keratinocyte - oral donor1'],
        109:['CNhs10881	','Nucleus Pulposus Cell donor1'],
        110:['CNhs10882','Prostate Epithelial Cells (polarized) donor1'],
        111:['CNhs10883','Prostate Stromal Cells donor1'],
        112:['CNhs10884','Small Airway Epithelial Cells donor1'],
        113:['CNhs11045','cord blood derived cell line:COBL-a untreated'],
        114:['CNhs11046','embryonic kidney cell line: HEK293/SLAM untreated'],
        115:['CNhs11047','embryonic kidney cell line: HEK293/SLAM infection 24hr'],
        116:['CNhs11049','cord blood derived cell line:COBL-a 24h infection(-C)'],
        117:['CNhs11050','cord blood derived cell line:COBL-a 24h infection'],
        118:['CNhs11051','Adipocyte - breast donor1'],
        119:['CNhs11052','Preadipocyte - breast donor1'],
        120:['CNhs11054','Adipocyte - omental donor1'],
        121:['CNhs11057','Mesenchymal Stem Cells - Wharton\'s Jelly donor1'],
        122:['CNhs11061','Gingival epithelial cells donor1 (GEA11)'],
        123:['CNhs11062','Dendritic Cells - monocyte immature derived donor1 tech_rep2'],
        124:['CNhs11063','Neural stem cells donor1'],
        125:['CNhs11064','Keratinocyte - epidermal donor1'],
        126:['CNhs11065','Preadipocyte - omental donor1'],
        127:['CNhs11073','Mast cell - stimulated donor1'],
        128:['CNhs11074','Fibroblast - skin spinal muscular atrophy donor1'],
        129:['CNhs11075','Whole blood (ribopure) donor090325 donation1'],
        130:['CNhs11076','Whole blood (ribopure) donor090325 donation2'],
        131:['CNhs11077','Mammary Epithelial Cell donor1'],
        132:['CNhs11079','Placental Epithelial Cells donor1'],
        133:['CNhs11082','Preadipocyte - visceral donor1'],
        134:['CNhs11083','Skeletal Muscle Cells donor1'],
        135:['CNhs11084','Skeletal muscle cells differentiated into Myotubes - multinucleated donor1'],
        136:['CNhs11085','Smooth Muscle Cells - Aortic donor1'],
        137:['CNhs11086','Smooth Muscle Cells - Brachiocephalic donor1'],
        138:['CNhs11087','Smooth Muscle Cells - Carotid donor1'],
        139:['CNhs11088','Smooth Muscle Cells - Coronary Artery donor1'],
        140:['CNhs11090','Smooth Muscle Cells - Subclavian Artery donor1'],
        141:['CNhs11091','Smooth Muscle Cells - Umbilical Artery donor1'],
        142:['CNhs11092','Tracheal Epithelial Cells donor1'],
        143:['CNhs11100','ductal cell carcinoma cell line:KLM-1'],
        144:['CNhs11183','schwannoma cell line:HS-PSS'],
        145:['CNhs11185','glioblastoma cell line:A172'],
        146:['CNhs11243','prostate cancer cell line:PC-3'],
        147:['CNhs11244','synovial sarcoma cell line:HS-SY-II'],
        148:['CNhs11245','schwannoma cell line:HS-PSS tech_rep2'],
        149:['CNhs11247','epithelioid sarcoma cell line:HS-ES-1'],
        150:['CNhs11248','glioblastoma cell line:A172 tech_rep2'],
        151:['CNhs11249','endometrial stromal sarcoma cell line:OMC-9'],
        152:['CNhs11250','chronic myelogenous leukemia cell line:K562'],
        153:['CNhs11251','acute lymphoblastic leukemia (B-ALL) cell line:BALL-1'],
        154:['CNhs11252','squamous cell carcinoma cell line:EC-GI-10'],
        155:['CNhs11253','acute lymphoblastic leukemia (T-ALL) cell line:Jurkat'],
        156:['CNhs11254','melanoma cell line:G-361'],
        157:['CNhs11255','rectal cancer cell line:TT1TKB'],
        158:['CNhs11256','gall bladder carcinoma cell line:TGBC14TKB'],
        159:['CNhs11257','renal cell carcinoma cell line:TUHR10TKB'],
        160:['CNhs11258','myeloma cell line:PCM6'],
        161:['CNhs11259','ductal cell carcinoma cell line:MIA Paca2'],
        162:['CNhs11260','prostate cancer cell line:DU145'],
        163:['CNhs11261','transitional-cell carcinoma cell line:JMSU1'],
        164:['CNhs11263','mesothelioma cell line:ACC-MESO-1'],
        165:['CNhs11264','mesothelioma cell line:ACC-MESO-4'],
        166:['CNhs11265','bile duct carcinoma cell line:TFK-1'],
        167:['CNhs11266','endometrial carcinoma cell line:OMC-2'],
        168:['CNhs11267','retinoblastoma cell line:Y79'],
        169:['CNhs11268','Burkitt\'s lymphoma cell line:RAJI'],
        170:['CNhs11269','rhabdomyosarcoma cell line:RMS-YM'],
        171:['CNhs11270','signet ring carcinoma cell line:NUGC-4'],
        172:['CNhs11271','hepatoma cell line:Li-7'],
        173:['CNhs11272','glioblastoma cell line:T98G'],
        174:['CNhs11273','squamous cell lung carcinoma cell line:EBC-1'],
        175:['CNhs11274','giant cell carcinoma cell line:LU65'],
        176:['CNhs11275','lung adenocarcinoma cell line:A549'],
        177:['CNhs11276','neuroblastoma cell line:CHP-134'],
        178:['CNhs11277','large cell lung carcinoma cell line:IA-LM'],
        179:['CNhs11279','osteosarcoma cell line:143B/TK^(-)neo^(R)'],
        180:['CNhs11280','colon carcinoma cell line:CACO-2'],
        181:['CNhs11281','melanoma cell line:COLO 679'],
        182:['CNhs11282','acute lymphoblastic leukemia (B-ALL) cell line:NALM-6'],
        183:['CNhs11283','cholangiocellular carcinoma cell line:HuH-28'],
        184:['CNhs11284','neuroblastoma cell line:NB-1'],
        185:['CNhs11285','small cell lung carcinoma cell line:LK-2'],
        186:['CNhs11286','gastric cancer cell line:AZ521'],
        187:['CNhs11287','oral squamous cell carcinoma cell line:HO-1-u-1'],
        188:['CNhs11288','cervical cancer cell line:D98-AH2'],
        189:['CNhs11289','cervical cancer cell line:ME-180'],
        190:['CNhs11290','osteosarcoma cell line:HS-Os-1'],
        191:['CNhs11303','Melanocyte - light donor1'],
        192:['CNhs11305','Smooth Muscle Cells - Aortic donor2'],
        193:['CNhs11309','Smooth Muscle Cells - Aortic donor3'],
        194:['CNhs11311','Osteoblast - differentiated donor1'],
        195:['CNhs11317','Pericytes donor1'],
        196:['CNhs11319','Fibroblast - Choroid Plexus donor1'],
        197:['CNhs11320','Meningeal Cells donor1'],
        198:['CNhs11321','Astrocyte - cerebellum donor1'],
        199:['CNhs11322','Fibroblast - Lymphatic donor1'],
        200:['CNhs11323','Esophageal Epithelial Cells donor1'],
        201:['CNhs11324','Smooth Muscle Cells - Esophageal donor1'],
        202:['CNhs11325','Alveolar Epithelial Cells donor1'],
        203:['CNhs11328','Smooth Muscle Cells - Bronchial donor1'],
        204:['CNhs11329','Smooth Muscle Cells - Tracheal donor1'],
        205:['CNhs11330','Renal Proximal Tubular Epithelial Cell donor1'],
        206:['CNhs11331','Renal Cortical Epithelial Cells donor1'],
        207:['CNhs11332','Renal Epithelial Cells donor1'],
        208:['CNhs11333','Renal Mesangial Cells donor1'],
        209:['CNhs11334','Urothelial Cells donor1'],
        210:['CNhs11335','Hepatic Stellate Cells (lipocyte) donor1'],
        211:['CNhs11336','Corneal Epithelial Cells donor1'],
        212:['CNhs11337','Keratocytes donor1'],
        213:['CNhs11338','Retinal Pigment Epithelial Cells donor1'],
        214:['CNhs11339','Fibroblast - Conjunctival donor1'],
        215:['CNhs11340','Trabecular Meshwork Cells donor1'],
        216:['CNhs11341','Amniotic Epithelial Cells donor1'],
        217:['CNhs11344','Mesenchymal Stem Cells - bone marrow donor1'],
        218:['CNhs11345','Mesenchymal Stem Cells - adipose donor1'],
        219:['CNhs11347','Mesenchymal Stem Cells - umbilical donor1'],
        220:['CNhs11349','Mesenchymal Stem Cells - amniotic membrane donor1'],
        221:['CNhs11350','Multipotent Cord Blood Unrestricted Somatic Stem Cells donor1'],
        222:['CNhs11351','Fibroblast - skin normal donor1'],
        223:['CNhs11352','Fibroblast - skin walker warburg donor1'],
        224:['CNhs11353','Fibroblast - skin dystrophia myotonica donor1'],
        225:['CNhs11354','Fibroblast - skin dystrophia myotonica donor2'],
        226:['CNhs11371','Adipocyte - subcutaneous donor2'],
        227:['CNhs11372','Chondrocyte - de diff donor2'],
        228:['CNhs11373','Chondrocyte - re diff donor2'],
        229:['CNhs11375','Endothelial Cells - Aortic donor2'],
        230:['CNhs11376','Endothelial Cells - Microvascular donor2'],
        231:['CNhs11377','Endothelial Cells - Vein donor2'],
        232:['CNhs11378','Fibroblast - Cardiac donor2'],
        233:['CNhs11379','Fibroblast - Dermal donor2'],
        234:['CNhs11381','Keratinocyte - epidermal donor2'],
        235:['CNhs11382','Mammary Epithelial Cell donor2'],
        236:['CNhs11383','Melanocyte - light donor2'],
        237:['CNhs11384','Neural stem cells donor2'],
        238:['CNhs11385','Osteoblast donor2'],
        239:['CNhs11386','Placental Epithelial Cells donor2'],
        240:['CNhs11671','Whole blood (ribopure) donor090309 donation2'],
        241:['CNhs11672','Whole blood (ribopure) donor090612 donation1'],
        242:['CNhs11673','Whole blood (ribopure) donor090612 donation2'],
        243:['CNhs11675','Whole blood (ribopure) donor090309 donation1'],
        244:['CNhs11676','uterus adult pool1'],
        245:['CNhs11677','salivary gland adult pool1'],
        246:['CNhs11680','lung fetal donor1'],
        247:['CNhs11714','chronic lymphocytic leukemia (T-CLL) cell line:SKW-3'],
        248:['CNhs11715','Hodgkin\'s lymphoma cell line:HD-Mar2'],
        249:['CNhs11716','papillary adenocarcinoma cell line:8505C'],
        250:['CNhs11717','oral squamous cell carcinoma cell line:HSC-3'],
        251:['CNhs11718','mesenchymal stem cell line:Hu5/E18'],
        252:['CNhs11722','leiomyoma cell line:10964C'],
        253:['CNhs11723','leiomyoma cell line:15242A'],
        254:['CNhs11724','leiomyoma cell line:15425'],
        255:['CNhs11725','argyrophil small cell carcinoma cell line:TC-YIK'],
        256:['CNhs11726','testicular germ cell embryonal carcinoma cell line:NEC8'],
        257:['CNhs11728','Wilms\' tumor cell line:HFWT'],
        258:['CNhs11729','myxofibrosarcoma cell line:MFH-ino'],
        259:['CNhs11731','embryonic pancreas cell line:1B2C6'],
        260:['CNhs11732','embryonic pancreas cell line:1C3D3'],
        261:['CNhs11733','embryonic pancreas cell line:1C3IKEI'],
        262:['CNhs11734','small-cell gastrointestinal carcinoma cell line:ECC4'],
        263:['CNhs11736','small cell gastrointestinal carcinoma cell line:ECC10'],
        264:['CNhs11737','gastric adenocarcinoma cell line:MKN1'],
        265:['CNhs11738','gastrointestinal carcinoma cell line:ECC12'],
        266:['CNhs11739','squamous cell carcinoma cell line:T3M-5'],
        267:['CNhs11740','granulosa cell tumor cell line:KGN'],
        268:['CNhs11741','diffuse large B-cell lymphoma cell line:CTB-1'],
        269:['CNhs11742','hepatoblastoma cell line:HuH-6'],
        270:['CNhs11744','neuroectodermal tumor cell line:FU-RPNT-1'],
        271:['CNhs11745','clear cell carcinoma cell line:JHOC-5'],
        272:['CNhs11746','serous adenocarcinoma cell line:JHOS-2'],
        273:['CNhs11747','carcinosarcoma cell line:JHUCS-1'],
        274:['CNhs11748','endometrioid adenocarcinoma cell line:JHUEM-1'],
        275:['CNhs11750','lens epithelial cell line:SRA 01/04'],
        276:['CNhs11752','mucinous adenocarcinoma cell line:JHOM-1'],
        277:['CNhs11753','neuroectodermal tumor cell line:FU-RPNT-2'],
        278:['CNhs11755','smooth muscle adult pool1'],
        279:['CNhs11756','pancreas adult donor1'],
        280:['CNhs11757','heart adult diseased post-infarction donor1'],
        281:['CNhs11758','heart adult diseased donor1'],
        282:['CNhs11760','aorta adult pool1'],
        283:['CNhs11761','blood adult pool1'],
        284:['CNhs11762','eye fetal donor1'],
        285:['CNhs11763','uterus fetal donor1'],
        286:['CNhs11764','spinal cord fetal donor1'],
        287:['CNhs11765','umbilical cord fetal donor1'],
        288:['CNhs11766','trachea fetal donor1'],
        289:['CNhs11768','tongue fetal donor1'],
        290:['CNhs11769','thyroid fetal donor1'],
        291:['CNhs11770','throat fetal donor1'],
        292:['CNhs11771','stomach fetal donor1'],
        293:['CNhs11772','temporal lobe fetal donor1 tech_rep1'],
        294:['CNhs11773','small intestine fetal donor1'],
        295:['CNhs11774','skin fetal donor1'],
        296:['CNhs11776','skeletal muscle fetal donor1'],
        297:['CNhs11777','rectum fetal donor1'],
        298:['CNhs11779','diaphragm fetal donor1'],
        299:['CNhs11780','colon fetal donor1'],
        300:['CNhs11781','duodenum fetal donor1 tech_rep1'],
        301:['CNhs11782','parietal lobe fetal donor1'],
        302:['CNhs11784','occipital lobe fetal donor1'],
        303:['CNhs11786','lung right lower lobe adult donor1'],
        304:['CNhs11787','occipital lobe adult donor1'],
        305:['CNhs11788','lymph node adult donor1'],
        306:['CNhs11789','left ventricle adult donor1'],
        307:['CNhs11790','left atrium adult donor1'],
        308:['CNhs11792','breast adult donor1'],
        309:['CNhs11794','colon adult donor1'],
        310:['CNhs11795','cerebellum adult pool1'],
        311:['CNhs11796','brain adult donor1'],
        312:['CNhs11797','brain fetal pool1'],
        313:['CNhs11798','liver fetal pool1'],
        314:['CNhs11810','oral squamous cell carcinoma cell line:SAS'],
        315:['CNhs11811','neuroblastoma cell line:NH-12'],
        316:['CNhs11812','small cell lung carcinoma cell line:WA-hT'],
        317:['CNhs11813','xeroderma pigentosum b cell line:XPL 17'],
        318:['CNhs11814','embryonic pancreas cell line:2C6'],
        319:['CNhs11818','neuroblastoma cell line:NBsusSR'],
        320:['CNhs11819','gastric adenocarcinoma cell line:MKN45'],
        321:['CNhs11820','choriocarcinoma cell line:T3M-3'],
        322:['CNhs11821','myxofibrosarcoma cell line:NMFH-1'],
        323:['CNhs11824','glassy cell carcinoma cell line:HOKUG'],
        324:['CNhs11825','large cell non-keratinizing squamous carcinoma cell line:SKG-II-SF'],
        325:['CNhs11827','serous cystadenocarcinoma cell line:HTOA'],
        326:['CNhs11828','tridermal teratoma cell line:HGRT'],
        327:['CNhs11829','sacrococcigeal teratoma cell line:HTST'],
        328:['CNhs11830','peripheral neuroectodermal tumor cell line:KU-SN'],
        329:['CNhs11832','pancreatic carcinoma cell line:NOR-P1'],
        330:['CNhs11834','carcinoid cell line:NCI-H1770'],
        331:['CNhs11835','osteoclastoma cell line:Hs 706.T'],
        332:['CNhs11836','Ewing\'s sarcoma cell line:Hs 863.T'],
        333:['CNhs11838','alveolar cell carcinoma cell line:SW 1573'],
        334:['CNhs11840','bronchioalveolar carcinoma cell line:NCI-H358'],
        335:['CNhs11841','bronchogenic carcinoma cell line:ChaGo-K-1'],
        336:['CNhs11842','fibrous histiocytoma cell line:GCT TIB-223'],
        337:['CNhs11843','hairy cell leukemia cell line:Mo'],
        338:['CNhs11844','acantholytic squamous carcinoma cell line:HCC1806'],
        339:['CNhs11845','biphenotypic B myelomonocytic leukemia cell line:MV-4-11'],
        340:['CNhs11846','carcinoid cell line:SK-PN-DW'],
        341:['CNhs11848','leiomyoblastoma cell line:G-402'],
        342:['CNhs11849','pharyngeal carcinoma cell line:Detroit 562'],
        343:['CNhs11851','liposarcoma cell line:SW 872'],
        344:['CNhs11852','lymphangiectasia cell line:DS-1'],
        345:['CNhs11853','neuroepithelioma cell line:SK-N-MC'],
        346:['CNhs11854','neurofibroma cell line:Hs 53.T'],
        347:['CNhs11856','pagetoid sarcoma cell line:Hs 925.T'],
        348:['CNhs11857','spindle cell sarcoma cell line:Hs 132.T'],
        349:['CNhs11858','mycosis fungoides T cell lymphoma cell line:HuT 102 TIB-162'],
        350:['CNhs11859','leukemia chronic megakaryoblastic cell line:MEG-01'],
        351:['CNhs11860','fibrosarcoma cell line:HT-1080'],
        352:['CNhs11861','medulloblastoma cell line:ONS-76'],
        353:['CNhs11862','bronchial squamous cell carcinoma cell line:KNS-62'],
        354:['CNhs11864','acute myeloid leukemia (FAB M2) cell line:NKM-1'],
        355:['CNhs11865','chronic myelogenous leukemia (CML) cell line:MEG-A2'],
        356:['CNhs11866','neuroectodermal tumor cell line:TASK1'],
        357:['CNhs11867','NK T cell leukemia cell line:KHYG-1'],
        358:['CNhs11868','hepatic mesenchymal tumor cell line:LI90'],
        359:['CNhs11869','somatostatinoma cell line:QGP-1'],
        360:['CNhs11870','liposarcoma cell line:KMLS-1'],
        361:['CNhs11872','thyroid carcinoma cell line:TCO-1'],
        362:['CNhs11873','mucinous cystadenocarcinoma cell line:MCAS'],
        363:['CNhs11875','choriocarcinoma cell line:SCH'],
        364:['CNhs11876','testicular germ cell embryonal carcinoma cell line:ITO-II'],
        365:['CNhs11877','rhabdomyosarcoma cell line:KYM-1'],
        366:['CNhs11878','teratocarcinoma cell line:NCC-IT-A3'],
        367:['CNhs11880','keratoacanthoma cell line:HKA-1'],
        368:['CNhs11881','anaplastic large cell lymphoma cell line:Ki-JK'],
        369:['CNhs11882','adenocarcinoma cell line:IM95m'],
        370:['CNhs11883','tubular adenocarcinoma cell line:SUIT-2'],
        371:['CNhs11884','teratocarcinoma cell line:NCR-G1'],
        372:['CNhs11885','small cell cervical cancer cell line:HCSC-1'],
        373:['CNhs11886','chronic myeloblastic leukemia (CML) cell line:KCL-22'],
        374:['CNhs11888','acute myeloid leukemia (FAB M7) cell line:MKPL-1'],
        375:['CNhs11889','anaplastic squamous cell carcinoma cell line:RPMI 2650'],
        376:['CNhs11890','teratocarcinoma cell line:PA-1'],
        377:['CNhs11891','hereditary spherocytic anemia cell line:WIL2-NS'],
        378:['CNhs11892','Wilms\' tumor cell line:G-401'],
        379:['CNhs11893','adrenal cortex adenocarcinoma cell line:SW-13'],
        380:['CNhs11894','normal embryonic palatal mesenchymal cell line:HEPM'],
        381:['CNhs11896','Gingival epithelial cells donor2 (GEA14)'],
        382:['CNhs11897','CD14+ monocyte derived endothelial progenitor cells donor2'],
        383:['CNhs11899','Macrophage - monocyte derived donor2'],
        384:['CNhs11900','Smooth Muscle Cells - Brain Vascular donor2'],
        385:['CNhs11901','Endothelial Cells - Lymphatic donor2'],
        386:['CNhs11902','Preadipocyte - omental donor2'],
        387:['CNhs11903','Gingival epithelial cells donor3 (GEA15)'],
        388:['CNhs11904','CD14+ monocyte derived endothelial progenitor cells donor3'],
        389:['CNhs11905','Neutrophils donor3'],
        390:['CNhs11906','Endothelial Cells - Lymphatic donor3'],
        391:['CNhs11907','Fibroblast - Periodontal Ligament donor3'],
        392:['CNhs11908','Myoblast donor3'],
        393:['CNhs11909','Fibroblast - Cardiac donor4'],
        394:['CNhs11911','Fibroblast - skin spinal muscular atrophy donor2'],
        395:['CNhs11912','Fibroblast - skin spinal muscular atrophy donor3'],
        396:['CNhs11913','Fibroblast - skin dystrophia myotonica donor3'],
        397:['CNhs11914','Fibroblast - skin normal donor2'],
        398:['CNhs11920','Smooth Muscle Cells - Prostate donor1'],
        399:['CNhs11923','Chondrocyte - de diff donor1'],
        400:['CNhs11925','Endothelial Cells - Microvascular donor1'],
        401:['CNhs11926','Endothelial Cells - Thoracic donor1'],
        402:['CNhs11927','Smooth Muscle Cells - Uterine donor3'],
        403:['CNhs11930','clear cell carcinoma cell line:TEN'],
        404:['CNhs11931','bone marrow stromal cell line:StromaNKtert'],
        405:['CNhs11932','basal cell carcinoma cell line:TE 354.T'],
        406:['CNhs11933','pleomorphic hepatocellular carcinoma cell line:SNU-387'],
        407:['CNhs11934','myelodysplastic syndrome cell line:SKM-1'],
        408:['CNhs11935','lymphoma malignant hairy B-cell cell line:MLMA'],
        409:['CNhs11943','breast carcinoma cell line:MCF7'],
        410:['CNhs11944','mixed mullerian tumor cell line:HTMMT'],
        411:['CNhs11945','meningioma cell line:HKBMM'],
        412:['CNhs11948','Whole blood (ribopure) donor090309 donation3'],
        413:['CNhs11949','Whole blood (ribopure) donor090612 donation3'],
        414:['CNhs11950','normal intestinal epithelial cell line:FHs 74 Int'],
        415:['CNhs11951','Sebocyte donor2'],
        416:['CNhs11952','Fibroblast - Gingival donor5 (GFH3)'],
        417:['CNhs11953','Fibroblast - Periodontal Ligament donor5 (PL30)'],
        418:['CNhs11954','CD14+ Monocytes donor2'],
        419:['CNhs11955','CD4+ T Cells donor2'],
        420:['CNhs11956','CD8+ T Cells donor2'],
        421:['CNhs11957','Natural Killer Cells donor2'],
        422:['CNhs11958','Peripheral Blood Mononuclear Cells donor2'],
        423:['CNhs11959','Neutrophils donor2'],
        424:['CNhs11960','Astrocyte - cerebral cortex donor2'],
        425:['CNhs11961','Fibroblast - Gingival donor2'],
        426:['CNhs11962','Fibroblast - Periodontal Ligament donor2'],
        427:['CNhs11963','Smooth Muscle Cells - Colonic donor2'],
        428:['CNhs11964','Skeletal Muscle Satellite Cells donor2'],
        429:['CNhs11965','Myoblast donor2'],
        430:['CNhs11966','Ciliary Epithelial Cells donor2'],
        431:['CNhs11967','Endothelial Cells - Umbilical vein donor2'],
        432:['CNhs11969','Adipocyte - breast donor2'],
        433:['CNhs11971','Preadipocyte - breast donor2'],
        434:['CNhs11972','Prostate Epithelial Cells donor2'],
        435:['CNhs11973','Prostate Stromal Cells donor2'],
        436:['CNhs11975','Small Airway Epithelial Cells donor2'],
        437:['CNhs11976','Smooth Muscle Cells - Prostate donor2'],
        438:['CNhs11977','Endothelial Cells - Artery donor2'],
        439:['CNhs11978','Endothelial Cells - Thoracic donor2'],
        440:['CNhs11979','Hair Follicle Dermal Papilla Cells donor2'],
        441:['CNhs11980','Osteoblast - differentiated donor2'],
        442:['CNhs11981','Preadipocyte - subcutaneous donor2'],
        443:['CNhs11982','Preadipocyte - visceral donor2'],
        444:['CNhs11987','Smooth Muscle Cells - Coronary Artery donor2'],
        445:['CNhs11988','Smooth Muscle Cells - Internal Thoracic Artery donor2'],
        446:['CNhs11989','Smooth Muscle Cells - Pulmonary Artery donor2'],
        447:['CNhs11990','Smooth Muscle Cells - Subclavian Artery donor2'],
        448:['CNhs11991','Smooth Muscle Cells - Umbilical Artery donor2'],
        449:['CNhs11992','Synoviocyte donor2'],
        450:['CNhs11993','Tracheal Epithelial Cells donor2'],
        451:['CNhs11996','Fibroblast - Periodontal Ligament donor6 (PLH3)'],
        452:['CNhs11997','CD14+ Monocytes donor3'],
        453:['CNhs11998','CD4+ T Cells donor3'],
        454:['CNhs11999','CD8+ T Cells donor3'],
        455:['CNhs12000','Dendritic Cells - monocyte immature derived donor3'],
        456:['CNhs12001','Natural Killer Cells donor3'],
        457:['CNhs12002','Peripheral Blood Mononuclear Cells donor3'],
        458:['CNhs12003','Macrophage - monocyte derived donor3'],
        459:['CNhs12004','Smooth Muscle Cells - Brain Vascular donor3'],
        460:['CNhs12005','Astrocyte - cerebral cortex donor3'],
        461:['CNhs12006','Fibroblast - Gingival donor3'],
        462:['CNhs12007','Smooth Muscle Cells - Colonic donor3'],
        463:['CNhs12008','Skeletal Muscle Satellite Cells donor3'],
        464:['CNhs12009','Ciliary Epithelial Cells donor3'],
        465:['CNhs12010','Endothelial Cells - Umbilical vein donor3'],
        466:['CNhs12011','Fibroblast - Aortic Adventitial donor3'],
        467:['CNhs12012','Mesothelial Cells donor3'],
        468:['CNhs12013','Preadipocyte - omental donor3'],
        469:['CNhs12014','Prostate Epithelial Cells donor3'],
        470:['CNhs12015','Prostate Stromal Cells donor3'],
        471:['CNhs12016','Small Airway Epithelial Cells donor3'],
        472:['CNhs12017','Adipocyte - subcutaneous donor3'],
        473:['CNhs12019','Nucleus Pulposus Cell donor2'],
        474:['CNhs12020','Chondrocyte - de diff donor3'],
        475:['CNhs12021','Chondrocyte - re diff donor3'],
        476:['CNhs12022','Endothelial Cells - Aortic donor3'],
        477:['CNhs12023','Endothelial Cells - Artery donor3'],
        478:['CNhs12024','Endothelial Cells - Microvascular donor3'],
        479:['CNhs12026','Endothelial Cells - Vein donor3'],
        480:['CNhs12027','Fibroblast - Cardiac donor3'],
        481:['CNhs12028','Fibroblast - Dermal donor3'],
        482:['CNhs12030','Hair Follicle Dermal Papilla Cells donor3'],
        483:['CNhs12031','Keratinocyte - epidermal donor3'],
        484:['CNhs12032','Mammary Epithelial Cell donor3'],
        485:['CNhs12033','Melanocyte - light donor3'],
        486:['CNhs12035','Osteoblast - differentiated donor3'],
        487:['CNhs12036','Osteoblast donor3'],
        488:['CNhs12037','Placental Epithelial Cells donor3'],
        489:['CNhs12038','Preadipocyte - subcutaneous donor3'],
        490:['CNhs12039','Preadipocyte - visceral donor3'],
        491:['CNhs12043','Smooth Muscle Cells - Brachiocephalic donor3'],
        492:['CNhs12044','Smooth Muscle Cells - Carotid donor3'],
        493:['CNhs12045','Smooth Muscle Cells - Coronary Artery donor3'],
        494:['CNhs12046','Smooth Muscle Cells - Internal Thoracic Artery donor3'],
        495:['CNhs12048','Smooth Muscle Cells - Subclavian Artery donor3'],
        496:['CNhs12049','Smooth Muscle Cells - Umbilical Artery donor3'],
        497:['CNhs12050','Synoviocyte donor3'],
        498:['CNhs12051','Tracheal Epithelial Cells donor3'],
        499:['CNhs12052','Fibroblast - Dermal donor4'],
        500:['CNhs12053','Skeletal Muscle Cells donor4'],
        501:['CNhs12054','Bronchial Epithelial Cell donor4'],
        502:['CNhs12055','Fibroblast - Dermal donor5'],
        503:['CNhs12056','Skeletal Muscle Cells donor5'],
        504:['CNhs12057','Fibroblast - Cardiac donor5'],
        505:['CNhs12058','Bronchial Epithelial Cell donor5'],
        506:['CNhs12059','Fibroblast - Dermal donor6'],
        507:['CNhs12060','Skeletal Muscle Cells donor6'],
        508:['CNhs12061','Fibroblast - Cardiac donor6'],
        509:['CNhs12062','Bronchial Epithelial Cell donor6'],
        510:['CNhs12063','Nucleus Pulposus Cell donor3'],
        511:['CNhs12064','Anulus Pulposus Cell donor2'],
        512:['CNhs12065','Preadipocyte - perirenal donor1'],
        513:['CNhs12067','Adipocyte - omental donor2'],
        514:['CNhs12068','Adipocyte - omental donor3'],
        515:['CNhs12069','Adipocyte - perirenal donor1'],
        516:['CNhs12074','Renal Glomerular Endothelial Cells donor1'],
        517:['CNhs12075','Hepatic Sinusoidal Endothelial Cells donor1'],
        518:['CNhs12079','Pericytes donor2'],
        519:['CNhs12080','Meningeal Cells donor2'],
        520:['CNhs12081','Astrocyte - cerebellum donor2'],
        521:['CNhs12084','Alveolar Epithelial Cells donor2'],
        522:['CNhs12086','Renal Glomerular Endothelial Cells donor2'],
        523:['CNhs12087','Renal Proximal Tubular Epithelial Cell donor2'],
        524:['CNhs12088','Renal Epithelial Cells donor2'],
        525:['CNhs12091','Urothelial Cells donor2'],
        526:['CNhs12092','Hepatic Sinusoidal Endothelial Cells donor2'],
        527:['CNhs12093','Hepatic Stellate Cells (lipocyte) donor2'],
        528:['CNhs12095','Keratocytes donor2'],
        529:['CNhs12100','Mesenchymal Stem Cells - bone marrow donor2'],
        530:['CNhs12104','Mesenchymal Stem Cells - amniotic membrane donor2'],
        531:['CNhs12105','Multipotent Cord Blood Unrestricted Somatic Stem Cells donor2'],
        532:['CNhs12117','Astrocyte - cerebellum donor3'],
        533:['CNhs12118','Fibroblast - Lymphatic donor3'],
        534:['CNhs12120','Renal Proximal Tubular Epithelial Cell donor3'],
        535:['CNhs12121','Renal Mesangial Cells donor3'],
        536:['CNhs12122','Urothelial Cells donor3'],
        537:['CNhs12123','Corneal Epithelial Cells donor3'],
        538:['CNhs12124','Trabecular Meshwork Cells donor3'],
        539:['CNhs12125','Amniotic Epithelial Cells donor3'],
        540:['CNhs12126','Mesenchymal Stem Cells - bone marrow donor3'],
        541:['CNhs12127','Mesenchymal Stem Cells - umbilical donor3'],
        542:['CNhs12227','spinal cord adult donor10252'],
        543:['CNhs12228','pineal gland adult donor10252'],
        544:['CNhs12229','pituitary gland adult donor10252'],
        545:['CNhs12310','medial temporal gyrus adult donor10252'],
        546:['CNhs12311','amygdala adult donor10252'],
        547:['CNhs12312','hippocampus adult donor10252'],
        548:['CNhs12314','thalamus adult donor10252'],
        549:['CNhs12315','medulla oblongata adult donor10252'],
        550:['CNhs12316','middle temporal gyrus donor10252'],
        551:['CNhs12317','parietal lobe adult donor10252'],
        552:['CNhs12318','substantia nigra adult donor10252'],
        553:['CNhs12319','globus pallidus adult donor10252'],
        554:['CNhs12320','occipital cortex adult donor10252'],
        555:['CNhs12321','caudate nucleus adult donor10252'],
        556:['CNhs12322','locus coeruleus adult donor10252'],
        557:['CNhs12323','cerebellum adult donor10252'],
        558:['CNhs12324','putamen adult donor10196'],
        559:['CNhs12325','epitheloid carcinoma cell line: HelaS3 ENCODE biol_rep1'],
        560:['CNhs12326','epitheloid carcinoma cell line: HelaS3 ENCODE biol_rep2'],
        561:['CNhs12327','epitheloid carcinoma cell line: HelaS3 ENCODE biol_rep3'],
        562:['CNhs12328','hepatocellular carcinoma cell line: HepG2 ENCODE biol_rep1'],
        563:['CNhs12329','hepatocellular carcinoma cell line: HepG2 ENCODE biol_rep2'],
        564:['CNhs12330','hepatocellular carcinoma cell line: HepG2 ENCODE biol_rep3'],
        565:['CNhs12331','B lymphoblastoid cell line: GM12878 ENCODE biol_rep1'],
        566:['CNhs12332','B lymphoblastoid cell line: GM12878 ENCODE biol_rep2'],
        567:['CNhs12333','B lymphoblastoid cell line: GM12878 ENCODE biol_rep3'],
        568:['CNhs12334','chronic myelogenous leukemia cell line:K562 ENCODE biol_rep1'],
        569:['CNhs12335','chronic myelogenous leukemia cell line:K562 ENCODE biol_rep2'],
        570:['CNhs12336','chronic myelogenous leukemia cell line:K562 ENCODE biol_rep3'],
        571:['CNhs12338','Neurons donor1'],
        572:['CNhs12339','Hair Follicle Outer Root Sheath Cells donor1'],
        573:['CNhs12340','Hepatocyte donor1'],
        574:['CNhs12341','Cardiac Myocyte donor1'],
        575:['CNhs12342','Lens Epithelial Cells donor1'],
        576:['CNhs12343','CD19+ B Cells donor1'],
        577:['CNhs12344','Fibroblast - Choroid Plexus donor2'],
        578:['CNhs12347','Hair Follicle Outer Root Sheath Cells donor2'],
        579:['CNhs12348','Smooth Muscle Cells - Bronchial donor2'],
        580:['CNhs12349','Hepatocyte donor2'],
        581:['CNhs12350','Cardiac Myocyte donor2'],
        582:['CNhs12351','testicular germ cell embryonal carcinoma cell line:NEC14'],
        583:['CNhs12352','CD19+ B Cells donor2'],
        584:['CNhs12354','CD19+ B Cells donor3'],
        585:['CNhs12362','testicular germ cell embryonal carcinoma cell line:NEC15'],
        586:['CNhs12363','mesenchymal precursor cell - adipose donor1'],
        587:['CNhs12364','mesenchymal precursor cell - adipose donor2'],
        588:['CNhs12365','mesenchymal precursor cell - adipose donor3'],
        589:['CNhs12366','mesenchymal precursor cell - bone marrow donor1'],
        590:['CNhs12367','mesenchymal precursor cell - bone marrow donor2'],
        591:['CNhs12368','mesenchymal precursor cell - cardiac donor1'],
        592:['CNhs12369','mesenchymal precursor cell - cardiac donor2'],
        593:['CNhs12370','mesenchymal precursor cell - cardiac donor3'],
        594:['CNhs12371','mesenchymal precursor cell - cardiac donor4'],
        595:['CNhs12372','mesenchymal precursor cell - ovarian cancer left ovary donor1'],
        596:['CNhs12373','mesenchymal precursor cell - ovarian cancer right ovary donor1'],
        597:['CNhs12374','mesenchymal precursor cell - ovarian cancer metastasis donor1'],
        598:['CNhs12375','mesenchymal precursor cell - ovarian cancer right ovary donor2'],
        599:['CNhs12376','mesenchymal precursor cell - ovarian cancer left ovary donor3'],
        600:['CNhs12377','mesenchymal precursor cell - ovarian cancer right ovary donor3 (SOC-57-02)'],
        601:['CNhs12378','mesenchymal precursor cell - ovarian cancer metastasis donor3'],
        602:['CNhs12379','amniotic membrane cells donor3'],
        603:['CNhs12380','chorionic membrane cells donor3'],
        604:['CNhs12492','Mesenchymal stem cells - umbilical donor0'],
        605:['CNhs12493','Fibroblast - Periodontal Ligament donor4 (PL29)'],
        606:['CNhs12494','Adipocyte - subcutaneous donor1'],
        607:['CNhs12495','Endothelial Cells - Aortic donor1'],
        608:['CNhs12496','Endothelial Cells - Artery donor1'],
        609:['CNhs12497','Endothelial Cells - Vein donor1'],
        610:['CNhs12498','Fibroblast - Cardiac donor1'],
        611:['CNhs12499','Fibroblast - Cardiac donor1'],
        612:['CNhs12501','Hair Follicle Dermal Papilla Cells donor1'],
        613:['CNhs12502','amniotic membrane cells donor1'],
        614:['CNhs12503','amniotic membrane cells donor2'],
        615:['CNhs12504','chorionic membrane cells donor1'],
        616:['CNhs12506','chorionic membrane cells donor2'],
        617:['CNhs12566','Mast cell donor1'],
        618:['CNhs12568','Lens Epithelial Cells donor2'],
        619:['CNhs12569','Smooth Muscle Cells - Umbilical Vein donor2'],
        620:['CNhs12570','Melanocyte - dark donor3'],
        621:['CNhs12571','Cardiac Myocyte donor3'],
        622:['CNhs12572','Lens Epithelial Cells donor3'],
        623:['CNhs12574','nasal epithelial cells donor2'],
        624:['CNhs12575','Basophils donor3'],
        625:['CNhs12588','CD34+ stem cells - adult bone marrow derived donor1 tech_rep1'],
        626:['CNhs12589','nasal epithelial cells donor1 tech_rep1'],
        627:['CNhs12592','Mast cell donor4'],
        628:['CNhs12593','Mast cell donor3'],
        629:['CNhs12594','Mast cell donor2'],
        630:['CNhs12596','Iris Pigment Epithelial Cells donor1'],
        631:['CNhs12597','Smooth Muscle Cells - Umbilical Vein donor1'],
        632:['CNhs12610','diencephalon adult'],
        633:['CNhs12611','olfactory region adult'],
        634:['CNhs12624','Renal Glomerular Endothelial Cells donor3'],
        635:['CNhs12626','Hepatocyte donor3'],
        636:['CNhs12639','tenocyte donor1'],
        637:['CNhs12640','tenocyte donor2'],
        638:['CNhs12641','tenocyte donor3'],
        639:['CNhs12726','Neurons donor2'],
        640:['CNhs12728','Renal Cortical Epithelial Cells donor2'],
        641:['CNhs12730','Mesenchymal Stem Cells - hepatic donor2'],
        642:['CNhs12731','Meningeal Cells donor3'],
        643:['CNhs12732','Renal Epithelial Cells donor3'],
        644:['CNhs12733','Retinal Pigment Epithelial Cells donor3'],
        645:['CNhs12805','medulloblastoma cell line:D283 Med'],
        646:['CNhs12806','large cell lung carcinoma cell line:NCI-H460'],
        647:['CNhs12807','plasma cell leukemia cell line:ARH-77'],
        648:['CNhs12808','small cell lung carcinoma cell line:DMS 144'],
        649:['CNhs12809','small cell lung carcinoma cell line:NCI-H82'],
        650:['CNhs12810','salivary acinar cells donor1'],
        651:['CNhs12811','salivary acinar cells donor2'],
        652:['CNhs12812','salivary acinar cells donor3'],
        653:['CNhs12838','merkel cell carcinoma cell line:MKL-1'],
        654:['CNhs12839','merkel cell carcinoma cell line:MS-1'],
        655:['CNhs12840','cerebral meninges adult'],
        656:['CNhs12842','appendix adult'],
        657:['CNhs12844','vein adult'],
        658:['CNhs12846','ductus deferens adult'],
        659:['CNhs12847','epididymis adult'],
        660:['CNhs12848','gall bladder adult'],
        661:['CNhs12849','parotid gland adult'],
        662:['CNhs12850','penis adult'],
        663:['CNhs12851','seminal vesicle adult'],
        664:['CNhs12852','submaxillary gland adult'],
        665:['CNhs12853','tongue adult'],
        666:['CNhs12854','vagina adult'],
        667:['CNhs12855','heart - mitral valve adult'],
        668:['CNhs12856','heart - pulmonic valve adult'],
        669:['CNhs12857','heart - tricuspid valve adult'],
        670:['CNhs12858','throat adult'],
        671:['CNhs12894','Smooth Muscle Cells - Tracheal donor3'],
        672:['CNhs12922','Mesenchymal Stem Cells - adipose donor3'],
        673:['CNhs12996','duodenum fetal donor1 tech_rep2'],
        674:['CNhs12997','temporal lobe fetal donor1 tech_rep2'],
        675:['CNhs12998','testis adult pool2'],
        676:['CNhs13049','acute myeloid leukemia (FAB M7) cell line:M-MOK'],
        677:['CNhs13050','acute myeloid leukemia (FAB M5) cell line:NOMO-1'],
        678:['CNhs13051','acute myeloid leukemia (FAB M5) cell line:P31/FUJ'],
        679:['CNhs13052','acute myeloid leukemia (FAB M2) cell line:Kasumi-6'],
        680:['CNhs13053','acute myeloid leukemia (FAB M0) cell line:KG-1'],
        681:['CNhs13054','acute myeloid leukemia (FAB M1) cell line:HYT-1'],
        682:['CNhs13055','acute myeloid leukemia (FAB M3) cell line:HL60'],
        683:['CNhs13056','acute myeloid leukemia (FAB M4eo) cell line:EoL-1'],
        684:['CNhs13057','acute myeloid leukemia (FAB M4eo) cell line:EoL-3'],
        685:['CNhs13058','acute myeloid leukemia (FAB M5) cell line:U-937 DE-4'],
        686:['CNhs13059','acute myeloid leukemia (FAB M6) cell line:EEB'],
        687:['CNhs13060','acute myeloid leukemia (FAB M6) cell line:F-36E'],
        688:['CNhs13061','mesothelioma cell line:NCI-H28'],
        689:['CNhs13062','mesothelioma cell line:NCI-H226'],
        690:['CNhs13063','mesothelioma cell line:NCI-H2052'],
        691:['CNhs13064','mesothelioma cell line:NCI-H2452'],
        692:['CNhs13066','mesothelioma cell line:Mero-25'],
        693:['CNhs13067','mesothelioma cell line:Mero-41'],
        694:['CNhs13068','mesothelioma cell line:Mero-48a'],
        695:['CNhs13069','mesothelioma cell line:Mero-82'],
        696:['CNhs13070','mesothelioma cell line:Mero-83'],
        697:['CNhs13072','mesothelioma cell line:Mero-84'],
        698:['CNhs13073','mesothelioma cell line:Mero-95'],
        699:['CNhs13074','mesothelioma cell line:No36'],
        700:['CNhs13075','mesothelioma cell line:ONE58'],
        701:['CNhs13080','Renal Glomerular Endothelial Cells donor4'],
        702:['CNhs13092','mesenchymal precursor cell - ovarian cancer left ovary donor2'],
        703:['CNhs13093','mesenchymal precursor cell - ovarian cancer metastasis donor2'],
        704:['CNhs13094','mesenchymal precursor cell - ovarian cancer left ovary donor4'],
        705:['CNhs13096','mesenchymal precursor cell - ovarian cancer right ovary donor4'],
        706:['CNhs13097','mesenchymal precursor cell - ovarian cancer metastasis donor4'],
        707:['CNhs13098','mesenchymal precursor cell - bone marrow donor3'],
        708:['CNhs13099','serous adenocarcinoma cell line:SK-OV-3-R biol_rep1'],
        709:['CNhs13195','CD4+CD25+CD45RA- memory regulatory T cells donor1'],
        710:['CNhs13202','CD4+CD25-CD45RA+ naive conventional T cells expanded donor1'],
        711:['CNhs13203','CD4+CD25+CD45RA+ naive regulatory T cells expanded donor1'],
        712:['CNhs13204','CD4+CD25+CD45RA- memory regulatory T cells expanded donor1'],
        713:['CNhs13205','CD4+CD25-CD45RA+ naive conventional T cells donor2'],
        714:['CNhs13206','CD4+CD25+CD45RA- memory regulatory T cells donor2'],
        715:['CNhs13207','CD14-CD16+ Monocytes donor2'],
        716:['CNhs13208','CD14+CD16+ Monocytes donor2'],
        717:['CNhs13215','CD4+CD25-CD45RA- memory conventional T cells expanded donor1'],
        718:['CNhs13216','CD14+CD16- Monocytes donor2'],
        719:['CNhs13223','CD4+CD25-CD45RA+ naive conventional T cells donor1'],
        720:['CNhs13224','CD14+CD16- Monocytes donor1'],
        721:['CNhs13449','optic nerve donor1'],
        722:['CNhs13454','skeletal muscle - soleus muscle donor1'],
        723:['CNhs13465','CD14+ monocytes - treated with BCG donor1'],
        724:['CNhs13466','CD14+ monocytes - treated with IFN + N-hexane donor1'],
        725:['CNhs13467','CD14+ monocytes - treated with Trehalose dimycolate (TDM) donor1'],
        726:['CNhs13468','CD14+ monocytes - mock treated donor1'],
        727:['CNhs13469','CD14+ monocytes - treated with Group A streptococci donor1'],
        728:['CNhs13470','CD14+ monocytes - treated with lipopolysaccharide donor1'],
        729:['CNhs13471','CD14+ monocytes - treated with Salmonella donor1'],
        730:['CNhs13472','CD14+ monocytes - treated with Cryptococcus donor1'],
        731:['CNhs13473','CD14+ monocytes - treated with Candida donor1'],
        732:['CNhs13474','CD14+ monocytes - treated with B-glucan donor1'],
        733:['CNhs13475','CD14+ monocytes - treated with BCG donor2'],
        734:['CNhs13476','CD14+ monocytes - treated with IFN + N-hexane donor2'],
        735:['CNhs13477','Hep-2 cells treated with Streptococci strain 5448 biol_rep1'],
        736:['CNhs13478','Hep-2 cells treated with Streptococci strain JRS4 biol_rep1'],
        737:['CNhs13479','Hep-2 cells mock treated biol_rep1'],
        738:['CNhs13480','immature langerhans cells donor2'],
        739:['CNhs13483','CD14+ monocytes - treated with Trehalose dimycolate (TDM) donor2'],
        740:['CNhs13484','CD14+ monocytes - mock treated donor2'],
        741:['CNhs13485','CD14+ monocytes - treated with Salmonella donor2'],
        742:['CNhs13487','CD14+ monocytes - treated with Cryptococcus donor2'],
        743:['CNhs13488','CD14+ monocytes - treated with Candida donor2'],
        744:['CNhs13489','CD14+ monocytes - treated with B-glucan donor2'],
        745:['CNhs13490','CD14+ monocytes - treated with IFN + N-hexane donor3'],
        746:['CNhs13491','CD14+ monocytes - mock treated donor3'],
        747:['CNhs13492','CD14+ monocytes - treated with Group A streptococci donor3'],
        748:['CNhs13493','CD14+ monocytes - treated with Salmonella donor3'],
        749:['CNhs13494','CD14+ monocytes - treated with Candida donor3'],
        750:['CNhs13495','CD14+ monocytes - treated with B-glucan donor3'],
        751:['CNhs13496','Hep-2 cells treated with Streptococci strain 5448 biol_rep2'],
        752:['CNhs13497','Hep-2 cells treated with Streptococci strain 5448 biol_rep3'],
        753:['CNhs13498','Hep-2 cells treated with Streptococci strain JRS4 biol_rep2'],
        754:['CNhs13499','Hep-2 cells treated with Streptococci strain JRS4 biol_rep3'],
        755:['CNhs13500','Hep-2 cells mock treated biol_rep2'],
        756:['CNhs13501','Hep-2 cells mock treated biol_rep3'],
        757:['CNhs13502','acute myeloid leukemia (FAB M2) cell line:Kasumi-1'],
        758:['CNhs13503','acute myeloid leukemia (FAB M4) cell line:FKH-1'],
        759:['CNhs13504','acute myeloid leukemia (FAB M4) cell line:HNT-34'],
        760:['CNhs13505','acute myeloid leukemia (FAB M6) cell line:F-36P'],
        761:['CNhs13507','mesenchymal precursor cell - ovarian cancer right ovary donor3 (SOC-57-02-G)'],
        762:['CNhs13508','serous adenocarcinoma cell line:SK-OV-3-R after co-culture with SOC-57-02-G biol_rep1'],
        763:['CNhs13512','CD4+CD25-CD45RA+ naive conventional T cells donor3'],
        764:['CNhs13513','CD4+CD25+CD45RA+ naive regulatory T cells donor3'],
        765:['CNhs13532','CD14+ monocytes - treated with Group A streptococci donor2'],
        766:['CNhs13533','CD14+ monocytes - treated with lipopolysaccharide donor2'],
        767:['CNhs13535','migratory langerhans cells donor1'],
        768:['CNhs13536','migratory langerhans cells donor2'],
        769:['CNhs13537','immature langerhans cells donor1'],
        770:['CNhs13538','CD4+CD25+CD45RA- memory regulatory T cells donor3'],
        771:['CNhs13539','CD4+CD25-CD45RA- memory conventional T cells donor3'],
        772:['CNhs13540','CD14+CD16- Monocytes donor3'],
        773:['CNhs13541','CD14+CD16+ Monocytes donor1'],
        774:['CNhs13543','CD14+ monocytes - treated with BCG donor3'],
        775:['CNhs13544','CD14+ monocytes - treated with Trehalose dimycolate (TDM) donor3'],
        776:['CNhs13545','CD14+ monocytes - treated with lipopolysaccharide donor3'],
        777:['CNhs13546','CD14+ monocytes - treated with Cryptococcus donor3'],
        778:['CNhs13547','migratory langerhans cells donor3'],
        779:['CNhs13548','CD14-CD16+ Monocytes donor3'],
        780:['CNhs13549','CD14+CD16+ Monocytes donor3'],
        781:['CNhs13550','Mallassez-derived cells donor2'],
        782:['CNhs13551','Mallassez-derived cells donor3'],
        783:['CNhs13552','Reticulocytes biol_ rep1'],
        784:['CNhs13553','Reticulocytes biol_ rep2'],
        785:['CNhs13793','amygdala - adult donor10196'],
        786:['CNhs13794','thalamus - adult donor10196'],
        787:['CNhs13795','hippocampus - adult donor10196'],
        788:['CNhs13796','medial frontal gyrus - adult donor10196'],
        789:['CNhs13797','parietal lobe - adult donor10196'],
        790:['CNhs13798','occipital cortex - adult donor10196'],
        791:['CNhs13799','cerebellum - adult donor10196'],
        792:['CNhs13800','medulla oblongata - adult donor10196'],
        793:['CNhs13801','globus pallidus - adult donor10196'],
        794:['CNhs13802','caudate nucleus - adult donor10196'],
        795:['CNhs13804','pineal gland - adult donor10196'],
        796:['CNhs13805','pituitary gland - adult donor10196'],
        797:['CNhs13807','spinal cord - adult donor10196'],
        798:['CNhs13808','locus coeruleus - adult donor10196'],
        799:['CNhs13809','medial temporal gyrus - adult donor10196'],
        800:['CNhs13811','CD4+CD25+CD45RA- memory regulatory T cells expanded donor2'],
        801:['CNhs13812','CD4+CD25+CD45RA- memory regulatory T cells expanded donor3'],
        802:['CNhs13813','CD4+CD25-CD45RA+ naive conventional T cells expanded donor2'],
        803:['CNhs13814','CD4+CD25-CD45RA+ naive conventional T cells expanded donor3'],
        804:['CNhs13815','Neurons donor3'],
        805:['CNhs13816','Olfactory epithelial cells donor1'],
        806:['CNhs13817','Olfactory epithelial cells donor2'],
        807:['CNhs13818','Olfactory epithelial cells donor3'],
        808:['CNhs13819','Olfactory epithelial cells donor4'],
    }
    if request.method == 'POST':
        return jsonify(roadmap,fantom_5)