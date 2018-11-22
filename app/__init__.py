from flask import Flask, url_for
from flask_login import LoginManager
from flask_sqlalchemy import SQLAlchemy
from importlib import import_module
from logging import basicConfig, DEBUG, getLogger, StreamHandler
from os import path

import requests, sys
import numpy as np
import pandas as pd
import re
from flask import jsonify
from flask import request

db = SQLAlchemy()
login_manager = LoginManager()


def register_extensions(app):
    db.init_app(app)
    login_manager.init_app(app)


def register_blueprints(app):
    for module_name in ('base', 'forms', 'ui', 'home', 'tables', 'data', 'additional', 'base'):
        module = import_module('app.{}.routes'.format(module_name))
        app.register_blueprint(module.blueprint)


def configure_database(app):

    @app.before_first_request
    def initialize_database():
        db.create_all()

    @app.teardown_request
    def shutdown_session(exception=None):
        db.session.remove()


def configure_logs(app):
    basicConfig(filename='error.log', level=DEBUG)
    logger = getLogger()
    logger.addHandler(StreamHandler())


def apply_themes(app):
    """
    Add support for themes.

    If DEFAULT_THEME is set then all calls to
      url_for('static', filename='')
      will modfify the url to include the theme name

    The theme parameter can be set directly in url_for as well:
      ex. url_for('static', filename='', theme='')

    If the file cannot be found in the /static/<theme>/ lcation then
      the url will not be modified and the file is expected to be
      in the default /static/ location
    """
    @app.context_processor
    def override_url_for():
        return dict(url_for=_generate_url_for_theme)

    def _generate_url_for_theme(endpoint, **values):
        if endpoint.endswith('static'):
            themename = values.get('theme', None) or \
                app.config.get('DEFAULT_THEME', None)
            if themename:
                theme_file = "{}/{}".format(themename, values.get('filename', ''))
                if path.isfile(path.join(app.static_folder, theme_file)):
                    values['filename'] = theme_file
        return url_for(endpoint, **values)

#DR_SCAN_get_snp_info
def drscan_features(app):
    @app.route('/background_process_test')
    def background_process_test():
        return 'teste'	

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

    @app.route('/get_snp_info',methods=['GET','POST'])
    def get_snp_info_web():
        snp = request.form['snp_input']
        snp1 = request.arg.get('snp_input')
        snp2 = request.form['snp_name']
        snp3 = request.arg.get('snp_name')
        print("SNP_INPUT----: ",snp)
        print("SNP_INPUT1: ",snp1)
        print("SNP_INPUT2: ",snp2)
        print("SNP_INPUT3: ",snp3)
        #snp_id = request.args.get('snp_name')
        #snp = get_snp_info(snp_id)
        #print("OK")
        a = get_snp_info('11897559')
        #TODO: return more than one dna version
        print(a[0][0])
        return jsonify(a[0][0])


def create_app(config, selenium=False):
    app = Flask(__name__, static_folder='base/static')
    app.config.from_object(config)
    if selenium:
        app.config['LOGIN_DISABLED'] = True
    register_extensions(app)
    register_blueprints(app)
    configure_database(app)
    configure_logs(app)
    apply_themes(app)
    return app
