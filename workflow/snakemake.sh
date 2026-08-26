export PYTHONNOUSERSITE=1
unset PYTHONPATH
snakemake --profile ../config/profile/ $@ 2> smk-"`date +"%d-%m-%Y"`".log
