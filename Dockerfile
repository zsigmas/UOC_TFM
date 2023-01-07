# 4.2.1 As it is the first one with a frozen date in Posit Package Manager
FROM ghcr.io/rocker-org/tidyverse:4.2.1

ADD docker_scripts/install_system_dep.sh docker_scripts/install_system_dep.sh
RUN chmod +x docker_scripts/install_system_dep.sh
RUN docker_scripts/install_system_dep.sh
ADD docker_scripts/install_script.R docker_scripts/install_script.R
RUN chmod +x docker_scripts/install_script.R
RUN docker_scripts/install_script.R
ADD docker_scripts/install_bioc_dep.R docker_scripts/install_bioc_dep.R
RUN chmod +x docker_scripts/install_bioc_dep.R
RUN docker_scripts/install_bioc_dep.R

# Fixes bug
ADD docker_scripts/install_fix.sh docker_scripts/install_fix.sh
RUN chmod +x docker_scripts/install_fix.sh
RUN ls docker_scripts/*
RUN docker_scripts/install_fix.sh

CMD ["sh"]