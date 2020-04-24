FROM python:3.6

ENV FLASK_APP gentelella.py

COPY gentelella.py gunicorn.py requirements.txt config.py .env ./
COPY app app
COPY migrations migrations
COPY TFs TFs
COPY Genome_GRCh37 Genome_GRCh37

RUN wget -c http://meme-suite.org/meme-software/5.1.0/meme-5.1.0.tar.gz
RUN tar -xvf meme-5.1.0.tar.gz
RUN ./meme-5.1.0/configure --prefix=$HOME/meme --with-url=http://meme-suite.org/ --enable-build-libxml2 --enable-build-libxslt

RUN make
RUN make install

RUN pip install --default-timeout=1000 future

RUN pip install -r requirements.txt


EXPOSE 5000
CMD ["gunicorn", "--config", "gunicorn.py", "gentelella:app"]