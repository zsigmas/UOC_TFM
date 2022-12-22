#!/usr/bin/bash
# Fixes: https://github.com/Bioconductor/bioconductor_docker/issues/22#issue-767695068
cd /
git clone https://github.com/bmbolstad/preprocessCore.git
cd preprocessCore
R CMD INSTALL --configure-args="--disable-threading"  .

