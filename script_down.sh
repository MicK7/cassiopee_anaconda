#!/bin/bash
version="2.4"
wget http://elsa.onera.fr/Cassiopee/Download/KCore-${version}F.tar.gz
wget http://elsa.onera.fr/Cassiopee/Download/Connector-${version}F.tar.gz
wget http://elsa.onera.fr/Cassiopee/Download/Converter-${version}F.tar.gz
wget http://elsa.onera.fr/Cassiopee/Download/Geom-${version}F.tar.gz
wget http://elsa.onera.fr/Cassiopee/Download/Transform-${version}F.tar.gz
wget http://elsa.onera.fr/Cassiopee/Download/Generator-${version}F.tar.gz
wget http://elsa.onera.fr/Cassiopee/Download/Post-${version}F.tar.gz
wget http://elsa.onera.fr/Cassiopee/Download/Initiator-${version}F.tar.gz
wget http://elsa.onera.fr/Cassiopee/Download/Distributor2-${version}F.tar.gz
wget http://elsa.onera.fr/Cassiopee/Download/Dist2Walls-${version}F.tar.gz

#wget http://elsa.onera.fr/Cassiopee/Download/Compressor-${version}F.tar.gz
#wget http://elsa.onera.fr/Cassiopee/Download/Intersector-${version}F.tar.gz
#
tar -xvf KCore-${version}F.tar.gz
tar -xvf Converter-${version}F.tar.gz
tar -xvf Geom-${version}F.tar.gz
tar -xvf Transform-${version}F.tar.gz
tar -xvf Generator-${version}F.tar.gz
tar -xvf Post-${version}F.tar.gz
tar -xvf Initiator-${version}F.tar.gz
tar -xvf Connector-${version}F.tar.gz
tar -xvf Distributor2-${version}F.tar.gz
tar -xvf Dist2Walls-${version}F.tar.gz
#tar -xvf Compressor-${version}F.tar.gz
#tar -xvf Intersector-${version}F.tar.gz

