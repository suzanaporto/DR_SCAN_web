from bcrypt import checkpw
from flask import jsonify, render_template, redirect, request, url_for
from flask_login import (
    current_user,
    login_required,
    login_user,
    logout_user
)

from app import db, login_manager
from app.base import blueprint
from app.base.forms import LoginForm, CreateAccountForm
from app.base.models import User


@blueprint.route('/')
def route_default():
    # return render_template('home_1.html')
    return redirect(url_for('base_blueprint.login'))

# @blueprint.route('/login')
# def login_default():
#     # return redirect(url_for('base_blueprint.login'))


@blueprint.route('/<template>')
@login_required
def route_template(template):
    return render_template(template + '.html')


@blueprint.route('/fixed_<template>')
@login_required
def route_fixed_template(template):
    return render_template('fixed/fixed_{}.html'.format(template))


@blueprint.route('/page_<error>')
def route_errors(error):
    return render_template('errors/page_{}.html'.format(error))

## Login & Registration


@blueprint.route('/login', methods=['GET', 'POST'])
def login():
    login_form = LoginForm(request.form)
    create_account_form = CreateAccountForm(request.form)
    if 'login' in request.form:
        username = request.form['username']
        password = request.form['password']
        user = User.query.filter_by(username=username).first()
        if user and checkpw(password.encode('utf8'), user.password):
            login_user(user)
            return redirect(url_for('home_blueprint.index'))
            # return redirect(url_for('base_blueprint.route_default'))
        return render_template('errors/page_403.html')
    if not current_user.is_authenticated:
        return render_template(
            'login/login.html',
            login_form=login_form,
            create_account_form=create_account_form
        )
    return redirect(url_for('home_blueprint.index'))


@blueprint.route('/create_user', methods=['POST'])
def create_user():
    user = User(**request.form)
    COUNTRY = [('1','Afghanistan'), 
    ('2','Albania'),
    ('3','Algeria'),
    ('4','American Samoa'),
    ('5','Andorra'),
    ('6','Angola'),
    ('7','Anguilla'),
    ('8','Antigua & Barbuda'),
    ('9','Argentina'),
    ('10','Armenia'),
    ('11','Aruba'),
    ('12','Australia'),
    ('13','Austria'),
    ('14','Azerbaijan'),
    ('15','Bahamas, The'),
    ('16','Bahrain'),
    ('17','Bangladesh'),
    ('18','Barbados '),
    ('19','Belarus '),
    ('20','Belgium '),
    ('21','Belize '),
    ('22','Benin '),
    ('23','Bermuda '),
    ('24','Bhutan '),
    ('25','Bolivia '),
    ('26','Bosnia & Herzegovina '),
    ('27','Botswana '),
    ('28','Brazil '),
    ('29','British Virgin Is. '),
    ('30','Brunei '),
    ('31','Bulgaria '),
    ('32','Burkina Faso '),
    ('33','Burma '),
    ('34','Burundi '),
    ('35','Cambodia '),
    ('36','Cameroon '),
    ('37','Canada '),
    ('38','Cape Verde '),
    ('39','Cayman Islands '),
    ('40','Central African Rep. '),
    ('40','Chad '),
    ('41','Chile '),
    ('42','China '),
    ('43','Colombia '),
    ('44','Comoros '),
    ('45','Congo, Dem. Rep. '),
    ('46','Congo, Repub. of the '),
    ('47','Cook Islands'),
    ('48','Costa Rica'),
    ('49','Cote dIvoire '),
    ('50','Croatia '),
    ('51','Cuba '),
    ('52','Cyprus '),
    ('53','Czech Republic '),
    ('54','Denmark '),
    ('55','Djibouti '),
    ('56','Dominica '),
    ('57','Dominican Republic '),
    ('58','East Timor '),
    ('59','Ecuador '),
    ('60','Egypt '),
    ('61','El Salvador '),
    ('62','Equatorial Guinea '),
    ('63','Eritrea '),
    ('64','Estonia '),
    ('65','Ethiopia '),
    ('66','Faroe Islands '),
    ('67','Fiji '),
    ('68','Finland '),
    ('69','France'),
    ('70','French Guiana'),
    ('71','French Polynesia '),
    ('72','Gabon '),
    ('73','Gambia, The '),
    ('74','Gaza Strip '),
    ('75','Georgia '),
    ('76','Germany '),
    ('77','Ghana '),
    ('78','Gibraltar '),
    ('79','Greece'),
    ('80','Greenland '),
    ('81','Grenada '),
    ('82','Guadeloupe '),
    ('83','Guam '),
    ('84','Guatemala '),
    ('85','Guernsey '),
    ('86','Guinea '),
    ('87','Guinea-Bissau '),
    ('88','Guyana '),
    ('89','Haiti '),
    ('90','Honduras '),
    ('91','Hong Kong '),
    ('92','Hungary '),
    ('93','Iceland '),
    ('94','India '),
    ('95','Indonesia '),
    ('96','Iran '),
    ('97','Iraq '),
    ('98','Ireland '),
    ('99','Isle of Man '),
    ('100','Israel '),
    ('101','Italy '),
    ('102','Jamaica '),
    ('103','Japan '),
    ('104','Jersey '),
    ('105','Jordan '),
    ('106','Kazakhstan '),
    ('107','Kenya '),
    ('108','Kiribati '),
    ('109','Korea, North '),
    ('110','Korea, South '),
    ('111','Kuwait '),
    ('112','Kyrgyzstan '),
    ('113','Laos '),
    ('114','Latvia '),
    ('115','Lebanon '),
    ('116','Lesotho '),
    ('117','Liberia '),
    ('118','Libya '),
    ('119','Liechtenstein '),
    ('120','Lithuania '),
    ('121','Luxembourg '),
    ('122','Macau '),
    ('123','Macedonia '),
    ('124','Madagascar '),
    ('125','Malawi '),
    ('126','Malaysia '),
    ('127','Maldives '),
    ('128','Mali '),
    ('129','Malta '),
    ('130','Marshall Islands '),
    ('131','Martinique '),
    ('132','Mauritania '),
    ('133','Mauritius '),
    ('134','Mayotte '),
    ('135','Mexico '),
    ('136','Micronesia, Fed. St.'),
    ('137','Moldova '),
    ('138','Monaco '),
    ('139','Mongolia '),
    ('140','Montserrat '),
    ('141','Morocco '),
    ('142','Mozambique '),
    ('143','Namibia '),
    ('144','Nauru '),
    ('145','Nepal '),
    ('146','Netherlands '),
    ('147','Netherlands Antilles '),
    ('148','New Caledonia',),
    ('149','New Zealand '),
    ('150','Nicaragua '),
    ('151','Niger '),
    ('152','Nigeria '),
    ('153','N. Mariana Islands '),
    ('154','Norway '),
    ('155','Oman '),
    ('156','Pakistan '),
    ('157','Palau '),
    ('158','Panama '),
    ('159','Papua New Guinea '),
    ('160','Paraguay '),
    ('161','Peru '),
    ('162','Philippines '),
    ('163','Poland '),
    ('164','Portugal '),
    ('165','Puerto Rico '),
    ('166','Qatar '),
    ('167','Reunion '),
    ('168','Romania '),
    ('169','Russia '),
    ('170','Rwanda '),
    ('171','Saint Helena '),
    ('172','Saint Kitts & Nevis '),
    ('173','Saint Lucia '),
    ('174','St Pierre & Miquelon '),
    ('175','Saint Vincent and the Grenadines '),
    ('176','Samoa '),
    ('177','San Marino '),
    ('178','Sao Tome & Principe '),
    ('179','Saudi Arabia '),
    ('180','Senegal '),
    ('181','Serbia '),
    ('182','Seychelles '),
    ('183','Sierra Leone '),
    ('184','Singapore '),
    ('185','Slovakia '),
    ('186','Slovenia '),
    ('187','Solomon Islands '),
    ('188','Somalia '),
    ('189','South Africa '),
    ('190','Spain '),
    ('191','Sri Lanka '),
    ('192','Sudan '),
    ('193','Suriname '),
    ('194','Swaziland '),
    ('195','Sweden '),
    ('196','Switzerland '),
    ('197','Syria '),
    ('198','Taiwan '),
    ('199','Tajikistan '),
    ('200','Tanzania '),
    ('201','Thailand '),
    ('202','Togo '),
    ('203','Tonga '),
    ('204','Trinidad & Tobago '),
    ('205','Tunisia '),
    ('206','Turkey '),
    ('207','Turkmenistan '),
    ('208','Turks & Caicos Is '),
    ('209','Tuvalu '),
    ('210','Uganda '),
    ('211','Ukraine '),
    ('212','United Arab Emirates '),
    ('213','United Kingdom '),
    ('214','United States '),
    ('215','Uruguay '),
    ('216','Uzbekistan '),
    ('217','Vanuatu '),
    ('218','Venezuela '),
    ('219','Vietnam '),
    ('220','Virgin Islands '),
    ('221','Wallis and Futuna '),
    ('222','West Bank '),
    ('223','Western Sahara '),
    ('224','Yemen'),
    ('225','Zambia'),
    ('226','Zimbabwe ')]
    countries = dict(COUNTRY)
    user.country = countries[str(user.country)]
    db.session.add(user)
    db.session.commit()
    return jsonify('success')


@blueprint.route('/logout')
@login_required
def logout():
    logout_user()
    return redirect(url_for('base_blueprint.login'))


@blueprint.route('/about')
def about():
    return render_template('login/about.html'),200

@blueprint.route('/privacy')
def privacy():
    return render_template('login/privacy.html'),200


@blueprint.route('/shutdown')
def shutdown():
    func = request.environ.get('werkzeug.server.shutdown')
    if func is None:
        raise RuntimeError('Not running with the Werkzeug Server')
    func()
    return 'Server shutting down...'

## Errors


@login_manager.unauthorized_handler
def unauthorized_handler():
    return render_template('errors/page_403.html'), 403


@blueprint.errorhandler(403)
def access_forbidden(error):
    return render_template('errors/page_403.html'), 403


@blueprint.errorhandler(404)
def not_found_error(error):
    return render_template('errors/page_404.html'), 404


@blueprint.errorhandler(500)
def internal_error(error):
    return render_template('errors/page_500.html'), 500
