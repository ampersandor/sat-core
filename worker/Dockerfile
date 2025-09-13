FROM python:3.13

# Install MAFFT
COPY ./lib/mafft-7.487-without-extensions-src.tgz .
RUN tar xvzf mafft-7.487-without-extensions-src.tgz && \
    cd mafft-7.487-without-extensions/core && \
    make clean && make && make install && \
    cd ../.. && rm -rf mafft-7.487-without-extensions*

# Copy UCLUST binary
COPY ./lib/uclust /usr/local/bin/
RUN chmod +x /usr/local/bin/uclust

WORKDIR /app

RUN pip install celery redis requests python-dotenv wonderwords biopython