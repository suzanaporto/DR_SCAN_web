FROM python:3.6

ENV FLASK_APP gentelella.py

COPY gentelella.py gunicorn.py requirements.txt config.py .env ./
COPY app app
COPY migrations migrations

RUN pip install -r requirements.txt

RUN wget -c http://meme-suite.org/meme-software/5.0.5/meme-5.0.5.tar.gz
RUN tar zxf meme-5.0.5.tar.gz
RUN cd meme-5.0.5
RUN ./configure --prefix=$HOME/meme --with-url=http://meme-suite.org/ --enable-build-libxml2 --enable-build-libxslt
RUN make
RUN make test
RUN make install

EXPOSE 5000
CMD ["gunicorn", "--config", "gunicorn.py", "gentelella:app"]