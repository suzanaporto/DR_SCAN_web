from os import environ


class Config(object):

    # Template deafult database
    # SQLALCHEMY_DATABASE_URI = 'sqlite:///database.db'
    SECRET_KEY = 'key'
    SQLALCHEMY_TRACK_MODIFICATIONS = False

    # Created database
    SQLALCHEMY_DATABASE_URI = 'postgresql://{}:{}@{}:{}/{}'.format(
        environ.get('GENTELELLA_DATABASE_USER', 'postgres'),
        environ.get('GENTELELLA_DATABASE_PASSWORD', 'unifor'),
        environ.get('GENTELELLA_DATABASE_HOST', '127.0.0.1'),
        environ.get('GENTELELLA_DATABASE_PORT', 5432),
        environ.get('GENTELELLA_DATABASE_NAME', 'regulomix')
    )

    # THEME SUPPORT
    #  if set then url_for('static', filename='', theme='')
    #  will add the theme name to the static URL:
    #    /static/<DEFAULT_THEME>/filename
    # DEFAULT_THEME = "themes/dark"
    DEFAULT_THEME = None


class ProductionConfig(Config):
    DEBUG = False

    # PostgreSQL database
    # SQLALCHEMY_DATABASE_URI = 'postgresql://{}:{}@{}:{}/{}'.format(
    #     environ.get('GENTELELLA_DATABASE_USER', 'gentelella'),
    #     environ.get('GENTELELLA_DATABASE_PASSWORD', 'gentelella'),
    #     environ.get('GENTELELLA_DATABASE_HOST', 'db'),
    #     environ.get('GENTELELLA_DATABASE_PORT', 5432),
    #     environ.get('GENTELELLA_DATABASE_NAME', 'gentelella')
    # )
    SQLALCHEMY_DATABASE_URI = 'postgresql://{}:{}@{}:{}/{}'.format(
        environ.get('GENTELELLA_DATABASE_USER', 'postgres'),
        environ.get('GENTELELLA_DATABASE_PASSWORD', 'unifor'),
        environ.get('GENTELELLA_DATABASE_HOST', '127.0.0.1'),
        environ.get('GENTELELLA_DATABASE_PORT', 5432),
        environ.get('GENTELELLA_DATABASE_NAME', 'regulomix')
    )

class DebugConfig(Config):
    DEBUG = True


config_dict = {
    'Production': ProductionConfig,
    'Debug': DebugConfig
}
