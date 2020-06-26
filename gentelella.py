from flask_migrate import Migrate
from os import environ
from sys import exit

from config import config_dict
from app import create_app, db, celery
import app

from app.celery_utils import init_celery

get_config_mode = environ.get('GENTELELLA_CONFIG_MODE', 'Debug')

try:
    config_mode = config_dict[get_config_mode.capitalize()]
except KeyError:
    exit('Error: Invalid GENTELELLA_CONFIG_MODE environment variable entry.')

app = create_app(config_mode,False, celery=app.celery)
init_celery(celery, app)

Migrate(app, db)
