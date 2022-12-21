docker run --name target_trial -v $(pwd)/:/workspace -w /workspace --rm uoc:latest Rscript -e "targets::tar_make()"
sudo chown -R zsigmas ./_targets

