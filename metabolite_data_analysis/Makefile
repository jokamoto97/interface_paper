main: setup run
setup:
	
	if [ ! -d "results/" ]; then mkdir results; fi

run:
	Rscript scripts/compare_pcg_implication_rst.R
	Rscript scripts/single_variant_pred_model_sensitivty_analysis.R
	Rscript scripts/h2_histogram.R
	Rscript scripts/analyze_isusie_effect_estimates.R

