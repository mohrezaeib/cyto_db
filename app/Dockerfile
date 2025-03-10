# Use posit base image on Ubuntu Jammy (22)
FROM rocker/shiny

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
    python3-venv\
    && apt-get clean

RUN python3 -m venv /venv
RUN /venv/bin/pip install rdkit

RUN R --version

# Copy your custom Shiny Server configuration
COPY shiny-server.conf /etc/shiny-server/shiny-server.conf

# Copy your Shiny app into the container
COPY cyto_db_shiny_app/ /srv/shiny-server/cyto_db_shiny_app/
COPY cyto_db_data/CytoLabs_Database.sdf /srv/shiny-server/cyto_db_shiny_app/
RUN R -e "install.packages('dplyr',repos = 'http://cran.rstudio.com/')"
RUN R -e "install.packages('here',repos = 'http://cran.rstudio.com/')"
RUN R -e "install.packages(c('reticulate', 'base64enc', 'png'),repos = 'http://cran.rstudio.com/' )"


# Expose Shiny Server port
EXPOSE 3838

# Expose RDKit Python library path
ENV PYTHONPATH="/usr/lib/python3/dist-packages"
ENV PATH="/venv/bin:$PATH"

# Start Shiny Server
CMD ["/usr/bin/shiny-server"]
