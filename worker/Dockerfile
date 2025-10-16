FROM python:3.13

# Install MAFFT
COPY ./lib/mafft-7.487-without-extensions-src.tgz .
RUN tar xvzf mafft-7.487-without-extensions-src.tgz && \
    cd mafft-7.487-without-extensions/core && \
    make clean && make && make install && \
    cd ../.. && rm -rf mafft-7.487-without-extensions*

RUN mkdir -p /app/tmp && chmod 777 /app/tmp
ENV MAFFT_TMPDIR=/app/tmp
ENV TMPDIR=/app/tmp

# Install UCLUST and VSEARCH
COPY ./lib/uclust /usr/local/bin/
RUN chmod +x /usr/local/bin/uclust

COPY ./lib/vsearch-2.19.0.tar.gz /opt/vsearch-2.19.0.tar.gz
RUN tar xvzf /opt/vsearch-2.19.0.tar.gz -C /opt/ && \
    cd /opt/vsearch-2.19.0 && \
    ./autogen.sh && \
    ./configure CFLAGS="-O2" CXXFLAGS="-O2" && \
    make ARFLAGS="cr" && make install && \
    rm -rf /opt/vsearch-2.19.0

WORKDIR /app


RUN pip install celery redis requests python-dotenv wonderwords biopython psutil