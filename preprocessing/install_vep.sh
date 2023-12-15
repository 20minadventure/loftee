#!/usr/bin/env bash

set -xe

mkdir -p $HOME/.vep/plugins
chmod a+rwx $HOME/.vep

docker pull ensemblorg/ensembl-vep:release_110.1
docker run -v $HOME/.vep:/data ensemblorg/ensembl-vep:release_110.1 INSTALL.pl --CACHEDIR /data --CACHE_VERSION 110 --AUTO cf --SPECIES homo_sapiens --ASSEMBLY GRCh38

git clone --depth 1 --branch v1.0.4_GRCh38 https://github.com/konradjk/loftee.git $HOME/.vep/plugins/loftee

wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz -P $HOME/.vep/plugins/loftee
wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.fai -P $HOME/.vep/plugins/loftee
wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.gzi -P $HOME/.vep/plugins/loftee
wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw -P $HOME/.vep/plugins/loftee

wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/loftee.sql.gz -P $HOME/.vep/plugins/loftee
gzip -dc $HOME/.vep/plugins/loftee/loftee.sql.gz > $HOME/.vep/plugins/loftee/loftee.sql
