# Use posit base image on Ubuntu Jammy (22)
FROM posit/r-base:4.3-jammy

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive




# Update and install necessary packages
RUN apt-get update && apt-get install -y \
    g++ \
    wget \
    curl \
    git \
    software-properties-common \
    build-essential \
    cmake \
    && apt-get clean

RUN apt-get update && apt-get install -y \
    libboost-all-dev \
    && apt-get clean

RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-dev \
    python3-numpy\
    && apt-get clean




# Install R packages
# RUN apt-get update && apt-get install -y \  
#     libssl-dev \
#     r-base\
#     r-base-dev \
#     gdebi-core \
#     && apt-get clean
# Verify R installation
RUN R --version
# Install RDKit
# RUN apt-get update && apt-get install -y \
#     librdkit-dev \
#     python3-rdkit \
#     && apt-get clean
RUN pip install rdkit

# Install R Shiny and Shiny Server
RUN R -e "install.packages('shiny', repos='https://cran.rstudio.com/')"
RUN apt-get update && apt-get install -y \
    wget \ 
    gdebi \
    && apt-get clean
RUN wget "https://download3.rstudio.org/ubuntu-20.04/x86_64/shiny-server-1.5.23.1030-amd64.deb"
RUN gdebi -n shiny-server-1.5.23.1030-amd64.deb
RUN rm shiny-server-1.5.23.1030-amd64.deb


# Copy your custom Shiny Server configuration
COPY shiny-server.conf /etc/shiny-server/shiny-server.conf

# Copy your Shiny app into the container
COPY cyto_db_shiny_app/ /srv/shiny-server/cyto_db_shiny_app/
COPY cyto_db_data/CytoLabs_Database.sdf /srv/shiny-server/cyto_db_shiny_app/
RUN Rscript /srv/shiny-server/cyto_db_shiny_app/requirements.R
RUN R -e "install.packages('dplyr',repos = 'http://cran.rstudio.com/')"
RUN R -e "install.packages('here',repos = 'http://cran.rstudio.com/')"
RUN R -e "install.packages(c('reticulate', 'base64enc', 'png'),repos = 'http://cran.rstudio.com/' )"


# Expose Shiny Server port
EXPOSE 3838

# Expose RDKit Python library path
ENV PYTHONPATH="/usr/lib/python3/dist-packages"

# Start Shiny Server
CMD ["/usr/bin/shiny-server"]
