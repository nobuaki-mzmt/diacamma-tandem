# DiacammaTandem
R codes for the analysis of Diacamma Tandem project

The directory is structured:
	Diacamma-tandem-DataAnalysis-R.Rproj:		Rproject file
	scripts:					All script for the analysis
		NetworkProcessing.R			Codes for processing ant recruitment datasets
		TrajectoryProcessing.R			Codes for processing tandem trajectory datasets
		NetworkOutput.R				Codes for the results of rectuitment analysis / network analysis
		TrajectoryOutputs.R			Codes for the results of trajectory analysis / trnsfer entropy analysis
		Phylogeny.R				Codes for phylogeentic comparative analysis
	data						All the data used for the analysis
		network
			processed			All .Rdata file (created during analysis) will be storaged here
			raw				Raw data for each species tandem 
				diacamma_indicum	Data from Kolay and Annagiri 2015, 10.1098/rsos.150104
				diacamma_sp		Diacamma cf indicum from Japan
				diacamma_albipennis	Data from Richardson et al.  2018, 10.1098/rspb.2017.2726
				diacamma_nylanderi	Data from Richardson et al.  2021, 10.1038/s42003-021-02048-7
				diacamma_rugatulus	Data from Valentini et al.   2020, 10.1098/rspb.2019.2950
		trajectory
			processed			All .Rdata file (created during analysis) will be storaged here
			raw				Raw data for each species tandem 
				Diacamma-csv		Tracking results of D. cf indicum in .csv file
		phylogeny				Include tree data and csv data for ant recruitment from Reeves and Moreau 2019
	img						All output will be storaged here
