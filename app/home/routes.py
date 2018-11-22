from app.home import blueprint
from flask import render_template
from flask_login import login_required
from flask import request

@blueprint.route('/index')
@login_required
def index():
    ###working method
    ###this is when you do a get request
    #if request.method == 'GET':
        #get parameter snp_name
        #snp_id = request.args.get('snp_name')
        #print snp id
        #print("SNP ID:",snp_id)
    
    return render_template('index.html')

@blueprint.route('/<template>')
@login_required
def route_template(template):
    return render_template(template + '.html')

