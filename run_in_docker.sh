docker run --name target_trial -v $(pwd)/:/workspace -w /workspace --rm uoc:final Rscript -e "targets::tar_make()"
echo "Remember the property of _targets is transferred to root due to running the pipeline in a docker!"
echo "Return property to the current user running: sudo chown -R ${USER} ./_targets"

