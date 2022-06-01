"Dockerfile" builds the cryanking/temperature_container on docker hub (digest:)

Data (if not included) is assumed to be in "ob_temp_working_data.csv" in the folder mounted as /research

The bland-altman style calculation is done by ba_analysis.R
- This largely follows https://pubmed.ncbi.nlm.nih.gov/32532218/
- outputs include 
    - agreement_output.csv with LOA calculations
    - oral_ir_ba.jpg, drager_oral_ba.jpg, drager_ir_ba.jpg BA type plots
    - pairwise_correlations.csv with raw correlations

The Hypothermia incidence calculation is done by hypothermia_incidence.R
- outputs include
    - smoothed_temp.csv with individual smoothed data using random slope spline models
    - curve_with_ci.csv with the overall temperature trends
    - main_outs.csv with the target quantity estimates and ci
    - scatter_smoothering.jpg figure showing smoothing effect
    - overlay_variation.jpeg figure showing sample of temporal changes
    - hypothermia_indicators.csv individual level occurance of smoothed severe hypothermia
    -  table1.csv patient characteristics stratified by severe hypothermia outcome
