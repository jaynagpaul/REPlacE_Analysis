FROM conda/miniconda2:latest

RUN apt-get update
RUN apt-get install -y git
RUN apt-get install -y wget

# Create Conda Env
RUN git clone https://github.com/jaynagpaul/REPlacE_Analysis/
RUN conda env create -f /REPlacE_Analysis/uditas_env.yml


# REFERENCE DATA
# RUN apt-get install -y tar
# RUN mkdir /reference

# RUN mkdir /reference/genome_2bit
# RUN wget -q -o /reference/genome_2bit/hg38.2bit https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit

# ENV BOWTIE2_INDEXES=/reference/bowtie2_index/
# RUN mkdir /reference/bowtie2_index/
# WORKDIR /reference/bowtie2_index/
# RUN wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz
# RUN tar zxvf GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz
# RUN mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.1.bt2 hg38.1.bt2
# RUN mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.2.bt2 hg38.2.bt2
# RUN mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.3.bt2 hg38.3.bt2
# RUN mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.4.bt2 hg38.4.bt2
# RUN mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.rev.1.bt2 hg38.rev.1.bt2
# RUN mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.rev.2.bt2 hg38.rev.2.bt2
# RUN chmod a+rx *
# WORKDIR /


# Activate conda environment workaround
SHELL ["/bin/bash", "-c"]
RUN echo "source activate uditas_env" > ~/.bashrc
ENV PATH /usr/local/envs/uditas_env/bin/:$PATH

WORKDIR /REPlacE_Analysis
RUN conda run -n uditas_env /bin/bash -c "python setup.py install"

WORKDIR /


