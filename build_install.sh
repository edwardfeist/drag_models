conda remove drag_models
conda build .
conda install --offline --use-local drag_models
conda build purge
