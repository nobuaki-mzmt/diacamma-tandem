# README

This repository provides access to the data and source code used for the manuscript

#### Cryptic functional diversity in ant tandem runs

by Nobuaki Mizumoto<sup>1</sup>, Yasunari Tanaka<sup>2</sup>, Gabriele Valentini<sup>3</sup>, Thomas O. Richardson<sup>4</sup>, Sumana Annagiri<sup>5</sup>, Stephen C Pratt<sup>6</sup>, Hiroyuki Shimoji<sup>2</sup>

<sup>1</sup> Okinawa Institute of Science & Technology Graduate University, Onna-son, Okinawa, Japan<br />
<sup>2</sup> School of Biological and Environmental Sciences, Kwansei Gakuin University, Sanda, Hyogo 669-1337, Japan<br />
<sup>3</sup> IRIDIA, Université Libre de Bruxelles, Brussels, Belgium<br />
<sup>4</sup> School of Biological Sciences, University of Bristol, UK <br />
<sup>5</sup> Behaviour and Ecology Lab, Department of Biological Sciences, Indian Institute of Science Education and Research, Mohanpur, India<br />
<sup>6</sup> School of Life Sciences, Arizona State University, Tempe, AZ 85287, USA<br />

## Table of Contents
* [README](./README.md) - this file
* [Rproject](./Rproject) - folder containing data and script of data analysis
	* [data](./data) - folder containing raw data and created rda data during analysis
	 	* [network](./data/network) - datasets for network / recruitment analysis
	 		* [raw](./data/network/raw) - Raw data for each species recruitment 
	 			* [diacamma_indicum](./data/network/raw/diacamma_indicum) - Data from [Kolay and Annagiri 2015](https://doi.org/10.1098/rsos.150104)
	 			* [diacamma_sp](./data/network/raw/diacamma_sp) - Data for Diacamma cf indicum from Japan
	 			* [temnothorax_albipennis](./data/network/raw/temnothorax_albipennis) - Data from [Richardson et al.  2018](https://doi.org/10.1098/rspb.2017.2726)
	 			* [temnothorax_nylanderi](./data/network/raw/temnothorax_nylanderi) - Data from [Richardson et al.  2021](https://doi.org/10.1038/s42003-021-02048-7)
	 			* [temnothorax_rugatulus](./data/network/raw/temnothorax_rugatulus) - Data from [Valentini et al.   2020](https://doi.org/10.1098/rspb.2019.2950)
	 		* [processed](./data/network/processed) - all .Rdata file (created during analysis) will be storaged here
	 	* [trajectory](./data/trajectory) - datasets for network / recruitment analysis
	 		* [raw](./data/trajectory/raw) - Raw data for each species tandem. Data from [Valentini et al. 2020](https://doi.org/10.7554/eLife.55395) and [Mizumoto and Dobata 2019](https://doi.org/10.1126/sciadv.aau6108), or this study [Diacamma-csv](./data/trajectory/raw/Diacamma-csv)
	 		* [processed](./data/trajectory/processed) - all .Rdata file (created during analysis) will be storaged here
	 		* [Scale.xlsx](./data/trajectory/Scale.xlsx) - Scale information for D. cf indicum experiments
	 		* [Tracking_failed.xlsx](./data/trajectory/Tracking_failed.xlsx) - Information for detection faulure during observation for D. cf indicum experiments
	 	* [phylogeny](./data/phylogeny) -  Include tree data and csv data for ant recruitment from [Reeves and Moreau 2019](https://doi.org/10.26049/ASP77-2-2019-10)
	* [scripts](./scripts) - folder containing all script for the analysis
		* [NetworkProcessing.R](./scripts/Output.R) - Codes for processing ant recruitment datasets
		* [TrajectoryProcessing.R](./scripts/Phylogeny.R) - Codes for processing tandem trajectory datasets
		* [NetworkOutput.R](./scripts/Preprocess.R) - Codes for the results of rectuitment analysis / network analysis
		* [TrajectoryOutputs.R](./scripts/Sources.R) - Codes for the results of trajectory analysis / trnsfer entropy analysis
		* [Phylogeny.R](./scripts/Sources.R) - Codes for phylogeentic comparative analysis
	* [img](./img) - All output will be storaged here
	* [Diacamma-tandem-DataAnalysis-R.Rproj](./Diacamma-tandem-DataAnalysis-R.Rproj) - R project file
* [TrackingAssist](./TrackingAssist) - folder containing source codes for tracking
	* [TrackingConverter.R](./TrackingAssist/TrackingConverter.R) - Scripts to convert the results of FastTrack into UMATracker
	* [TrimFocusedTandem.html](./TrackingAssist/TrimFocusedTandem.html) - python scripts to creat masked videos

## Usage notes
To reproduce our analysis, open [Diacamma-tandem-DataAnalysis-R.Rproj](./Diacamma-tandem-DataAnalysis-R.Rproj) in [RStudio](https://www.rstudio.com/). 
[NetworkProcessing.R](./scripts/NetworkProcessing.R) needs to run before [NetworkOutput.R](./scripts/NetworkOutput.R), and
[TrajectoryProcessing.R](./scripts/TrajectoryProcessing.R) needs to run before [TrajectoryOutputs.R](./scripts/TrajectoryOutputs.R).
For each script, load all of the function first. Results will be produced in [img](./img) or printed in the Console.

## Recruitment analysis of D. cf indicum
We recorded individual IDs and role of tandem runs during nest emigration of D. cf indidum. Data are manually created from the observation of videos.

## Tandem analysis of D. cf indicum
We recorded individual movement patterns of tandem running pairs in D. cf indicum, combining two different tracking softwares, [UMATracker](https://ymnk13.github.io/UMATracker/) and [FastTrack](https://www.fasttrack.sh/docs/interactiveTracking/), and our python coding in [TrackingAssist](./TrackingAssist).
Coordinates files were storaged in [Diacamma-csv](./data/trajectory/raw/Diacamma-csv).

See the Method and Supplementary Materials for more details.
