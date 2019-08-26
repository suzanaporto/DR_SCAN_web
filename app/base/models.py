from bcrypt import gensalt, hashpw
from flask_login import UserMixin
from sqlalchemy import Binary, Column, Integer, String, DateTime, ForeignKey
from sqlalchemy.orm import relationship
import datetime

from app import db, login_manager


class User(db.Model, UserMixin):

    __tablename__ = 'user'

    id = Column(Integer, primary_key=True)
    username = Column(String, unique=True)
    email = Column(String, unique=True)
    password = Column(Binary)
    country = Column(String, unique=False)
    workflow = relationship('Workflow', backref='user', lazy='dynamic')

    def __init__(self, **kwargs):
        for property, value in kwargs.items():
            # depending on whether value is an iterable or not, we must
            # unpack it's value (when **kwargs is request.form, some values
            # will be a 1-element list)
            if hasattr(value, '__iter__') and not isinstance(value, str):
                # the ,= unpack of a singleton fails PEP8 (travis flake8 test)
                value = value[0]
            if property == 'password':
                value = hashpw(value.encode('utf8'), gensalt())
            setattr(self, property, value)

    def __repr__(self):
        return str(self.username)


@login_manager.user_loader
def user_loader(id):
    return User.query.filter_by(id=id).first()


@login_manager.request_loader
def request_loader(request):
    username = request.form.get('username')
    user = User.query.filter_by(username=username).first()
    return user if user else None

class Workflow(db.Model):

    __tablename__ = 'workflow'
    # __table_args__ = {'extend_existing': True} 

    id_workflow = Column(Integer, primary_key=True)
    created = Column(DateTime,index=True, default=datetime.datetime.utcnow)
    user_id_user = Column(Integer, ForeignKey('user.id'))
    step1 = Column(String(128))
    step2 = Column(String(128))
    step3 = Column(String(128))
    step4 = Column(String(128))
    step5 = Column(String(128))
    expire = Column(DateTime)

    def __repr__(self):
        return '<Workflow {}>'.format(self.id_workflow)