# 4.2.1 As it is the first one with a frozen date
FROM ghcr.io/rocker-org/tidyverse:4.2.1

ADD docker_scripts/install_script.R docker_scripts/install_script.R
RUN chmod +x docker_scripts/install_script.R
RUN docker_scripts/install_script.R
ADD docker_scripts/install_bioc_dep.sh docker_scripts/install_bioc_dep.sh
RUN chmod +x docker_scripts/install_bioc_dep.sh
RUN docker_scripts/install_bioc_dep.sh
ADD docker_scripts/install_clean.R docker_scripts/install_clean.R
RUN chmod +x docker_scripts/install_clean.R
RUN docker_scripts/install_clean.R
ADD docker_scripts/install_bioc_dep.R docker_scripts/install_bioc_dep.R
RUN chmod +x docker_scripts/install_bioc_dep.R
RUN docker_scripts/install_bioc_dep.R

RUN apt-get install -y libxtst6 libxt6 && ln -s /usr/local/lib/R/lib/libR.so /lib/x86_64-linux-gnu/libR.so

CMD ["sh"]