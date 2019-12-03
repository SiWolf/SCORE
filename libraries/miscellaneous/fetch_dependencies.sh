# Script for downloading SCORE dependencies

# Download dependencies
wget ftp://ftp.jcvi.org/pub/data/TIGRFAMs/TIGRFAMs_15.0_HMM.LIB.gz -P references/tigrfam/
wget ftp://ftp.jcvi.org/pub/data/TIGRFAMs/TIGRFAMS_ROLE_LINK -P references/tigrfam/
wget ftp://ftp.jcvi.org/pub/data/TIGRFAMs/TIGR_ROLE_NAMES -P references/tigrfam/
wget https://genescf.kandurilab.org/ftp/geneSCF-master-v1.1-p2.tar.gz -P libraries/
wget http://purl.obolibrary.org/obo/go.obo -P references/

# Unpack and install them
gunzip references/tigrfam/TIGRFAMs_15.0_HMM.LIB.gz 

tar -xvzf libraries/geneSCF-master-v1.1-p2.tar.gz

mv geneSCF-master-source-v1.1-p2 libraries/geneSCF/

# Final clean up
rm -r libraries/geneSCF-master-v1.1-p2.tar.gz