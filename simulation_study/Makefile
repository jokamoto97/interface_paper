main: setup run
setup:

	if [ ! -d "sim_rst/" ]; then mkdir sim_rst; fi

run:
	Rscript scripts/power_fdr_rmse.R
	Rscript scripts/effect_scatterplots.R 
	Rscript scripts/pve_histogram.R
