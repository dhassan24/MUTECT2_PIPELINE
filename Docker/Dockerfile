# Base image: Ubuntu 18.04 LTS - provides a stable, well-supported Linux environment
FROM ubuntu:18.04
########################## DO NO EDIT ABOVE THIS LINE ##########################

# MAINTAINER: identifies the author of this Dockerfile for collaborators
MAINTAINER Danya Hassan <dhassan24@kgi.edu>

###############################################################################
###### Environment settings and dependency installation - DO NOT EDIT #########
###############################################################################

# Set terminal type for proper text rendering in interactive sessions
ENV TERM=vt100

# Prevent interactive prompts during apt package installation (non-interactive mode)
ARG DEBIAN_FRONTEND=noninteractive

# Update apt package lists and install core system dependencies:
# - build-essential: gcc, make, and other compilation tools
# - dialog: UI for shell scripts
# - wget/curl: file downloads
# - libx11/x11 libs: required by some genomics tools with GUI components
# - R base and dev libraries: needed for FACETS (R package)
# - python3.7: scripting and GATK dependencies
# - java (default-jdk): required to run GATK4 (JVM-based)
# - zlib, libbz2, liblzma: compression libraries for BAM/VCF handling
# --no-install-recommends keeps the image lean
RUN apt-get update && \
  apt-get -y --no-install-recommends install apt-utils build-essential dialog \
  wget curl git unzip \
  libx11-dev x11-common \
  r-base r-base-dev \
  python3.7 python3-pip \
  default-jdk \
  zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev \
  libncurses5-dev libncursesw5-dev

# Install LaTeX packages required for R to generate PDF reports (e.g., from FACETS)
RUN apt-get -y install texinfo texlive-latex-base texlive-fonts-recommended \
  texlive-fonts-extra texlive-latex-extra

# Create a symlink so that 'python' calls python3.7 by default
# (Ubuntu 18.04 defaults to python2; this ensures modern python is used)
RUN ln -sf /usr/bin/python3.7 /usr/bin/python

###############################################################################
#################### Take home assignment starts here #########################
###############################################################################

# ── GATK4 Installation ───────────────────────────────────────────────────────
# Download the latest GATK4 release from the Broad Institute's GitHub releases.
# We use wget to fetch the zip archive, unzip it into /opt, and create a symlink
# so that 'gatk' is accessible system-wide without specifying the full path.
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip \
    -O /tmp/gatk.zip && \
    unzip /tmp/gatk.zip -d /opt/ && \
    rm /tmp/gatk.zip

# Add GATK to PATH so it can be called by just typing 'gatk' in the shell
ENV PATH="/opt/gatk-4.5.0.0:${PATH}"

# ── Samtools Installation ─────────────────────────────────────────────────────
# Samtools is required for BAM indexing and manipulation (used alongside GATK).
# We download the source, compile it, and install to /usr/local/bin.
RUN wget https://github.com/samtools/samtools/releases/download/1.19/samtools-1.19.tar.bz2 \
    -O /tmp/samtools.tar.bz2 && \
    tar -xjf /tmp/samtools.tar.bz2 -C /tmp/ && \
    cd /tmp/samtools-1.19 && \
    ./configure --prefix=/usr/local && \
    make && make install && \
    rm -rf /tmp/samtools*

# ── FACETS R Package Installation ─────────────────────────────────────────────
# FACETS is an R package for allele-specific copy number analysis.
# It requires pctGCdata as a dependency.
# We install both from GitHub using devtools inside R.
# - pctGCdata: GC content reference data required by FACETS
# - facets: the main copy number analysis package
# RUN Rscript -e "install.packages('devtools', repos='http://cran.rstudio.com/')" && \
   # Rscript -e "install.packages('BiocManager', repos='http://cran.rstudio.com/')" && \
   # Rscript -e "BiocManager::install('IRanges')" && \
   # Rscript -e "devtools::install_github('mskcc/pctGCdata')" && \
   # Rscript -e "devtools::install_github('mskcc/facets')"

# Install Python dependencies for filtering pipeline
####### RUN pip3 install pandas cyvcf2 numpy

# Install pip, setuptools, and Python dependencies
# setuptools is required for building packages like cyvcf2 from source

# Install pip, setuptools, and Python dependencies
# - python3-pip: provides pip3 package manager
# - python3-dev: Python header files needed to compile C extensions (cyvcf2)
# - python3-setuptools: required by cyvcf2 to compile from source
# - pandas: dataframe handling for FACETS seg/qc files
# - numpy: numerical operations
# - cyvcf2: fast VCF parsing (compiled C extension, needs python3-dev headers)

RUN apt-get -y install python3-pip python3-dev python3-setuptools && \
    pip3 install setuptools --upgrade && \
    pip3 install pandas numpy && \
    pip3 install cyvcf2

# ── facets-suite Installation ──────────────────────────────────────────────────
# facets-suite provides command-line scripts (runFacets.R, annotate, etc.)
# that wrap the FACETS R package for easier pipeline use.
# Cloning from GitHub and adding to PATH makes scripts directly executable.
# Clone facets-suite and make any executable scripts accessible
# Using find to chmod only files that actually exist
RUN git clone https://github.com/mskcc/facets-suite.git /opt/facets-suite && \
    find /opt/facets-suite -name "*.py" -exec chmod +x {} \; && \
    find /opt/facets-suite -name "*.R"  -exec chmod +x {} \;


# Add facets-suite scripts to PATH so they're callable by name
ENV PATH="/opt/facets-suite:${PATH}"

# Install Rhtslib R package first — required to compile snp-pileup
# Rhtslib provides the htslib C headers needed for the g++ compilation step
# RUN Rscript -e "install.packages('BiocManager', repos='http://cran.rstudio.com/')" && \
   # Rscript -e "BiocManager::install('Rhtslib', ask=FALSE)"


# Compile snp-pileup from FACETS source using system htslib headers
# - libboost-all-dev: C++ boost libraries required by snp-pileup
# - libhts-dev: provides htslib headers (bgzf.h, sam.h) AND the compiled library
# snp-pileup generates allele counts at SNP positions from BAM files
# which is the required input format for FACETS copy number analysis
RUN apt-get -y install libboost-all-dev libhts-dev && \
    git clone https://github.com/mskcc/facets.git /tmp/facets-src && \
    cd /tmp/facets-src/inst/extcode && \
    g++ -std=c++11 -o /usr/local/bin/snp-pileup snp-pileup.cpp \
        -lhts && \
    rm -rf /tmp/facets-src

###############################################################################
#################### Take home assignment ends here ###########################
###############################################################################

######################### DO NO EDIT BELOW THIS LINE ##########################

# Set default Rscript flags: --vanilla suppresses loading of .Rprofile,
# .Renviron, history, and workspace — ensures reproducible, clean R execution
ENV Rscript="Rscript --vanilla"

# Default command: drop into bash shell when container is run interactively
CMD ["/bin/bash"]
