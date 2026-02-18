# ROI-RFA_MedSea
This repo contains the Matlab scripts to run the analysis presented in:

De Leo, F., & Solari, S. (2025). Regional frequency analysis of extreme waves based on regions of influence in the Mediterranean Sea. Ocean Engineering, 342, 123088.

func_find_peaks.m --> is an example function for extraction of a series of peaks from a time series of significant wave heights

run_MeanShift_Clustering.m --> once you have the series of peaks for each node, this script defines the events (clusters) using the MeanShift++ algorithm, and finds the annual maxima at each node

run_H_Test.m --> estimates the H statstic for a given node, for a series of posible p_ROI values. It calls func_cmp_H.m, that requieres R to be installed.

run_RFA_GEV.m --> performs EVA following both RFA and at-site analysis for a target node. It uses func_LMOM_GEV.m and func_LMOM_GEV_ROI.m

The archive containing the AM series of SWH in the Mediterranean Sea is freely available for research and can be requested by contacting francesco.deleo@unige.it
