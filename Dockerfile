# Use Rocker R base image
FROM rocker/r-ver:4.3.1

# Install necessary system libraries
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libxt-dev \
    && apt-get clean

# Install required R packages
RUN Rscript -e "install.packages(c('data.table', 'edgeR', 'limma', 'ggplot2', 'ggrepel', 'clusterProfiler', 'org.Hs.eg.db', 'ReactomePA', 'enrichplot', 'ggnewscale', 'biomaRt'), repos='https://cloud.r-project.org/')"

# Set working directory
WORKDIR /app

# Copy the R script and data file into the container
COPY STML_TXNIP_Script_DEA_GSE_analysis.md /app/STML_TXNIP_Script_DEA_GSE_analysis.md
COPY TXNIP_raw_counts.csv /app/TXNIP_raw_counts.csv

# Set the default command to run the R script
CMD ["Rscript", "/app/STML_TXNIP_Script_DEA_GSE_analysis.md"]

