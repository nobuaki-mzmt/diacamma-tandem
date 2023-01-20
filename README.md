# README

This repository provides access to the data and source code used for the manuscript  
**Cryptic functional diversity in ant tandem runs**  
Nobuaki Mizumoto<sup>1</sup>, Yasunari Tanaka<sup>2</sup>, Gabriele Valentini<sup>3</sup>, Thomas O. Richardson<sup>4</sup>, Sumana Annagiri<sup>5</sup>, Stephen C Pratt<sup>6</sup>, Hiroyuki Shimoji<sup>2</sup>  

<sup>1</sup> Okinawa Institute of Science & Technology Graduate University, Onna-son, Okinawa, Japan  
<sup>2</sup> School of Biological and Environmental Sciences, Kwansei Gakuin University, Sanda, Hyogo 669-1337, Japan  
<sup>3</sup> IRIDIA, Université Libre de Bruxelles, Brussels, Belgium  
<sup>4</sup> School of Biological Sciences, University of Bristol, UK  
<sup>5</sup> Behaviour and Ecology Lab, Department of Biological Sciences, Indian Institute of Science Education and Research, Mohanpur, India  
<sup>6</sup> School of Life Sciences, Arizona State University, Tempe, AZ 85287, USA  

## Usage notes
### Tracking
Trajectories of tandem running pairs in D. cf indicum were obtained by combining two different tracking softwares, [UMATracker](https://ymnk13.github.io/UMATracker/) and [FastTrack](https://www.fasttrack.sh/docs/interactiveTracking/), with python coding [TrimFocusedTandem.py](./tracking/TrimFocusedTandem.py) and R coding [TrackingConverter.R](./tracking/test/TrackingConverter.R).
Produced files for coordinates were stored in [Diacamma-csv](./data/trajectory/raw/Diacamma-csv).

### Data analysis
#### Trajectory analysis
[TrajectoryProcessing.R](./scripts/TrajectoryProcessing.R) produces Rdata in [data/trajectory/processed](./data/trajectory/processed) based on raw data in [data/trajectory/raw](./data/trajectory/raw).  
Then, [TrajectoryOutputs.R](./scripts/TrajectoryOutputs.R) outputs results in [img/trajectory](./img/trajectory).  
  
#### Network analysis
[NetworkProcessing.R](./scripts/NetworkProcessing.R) produces Rdata in [data/network/processed](./data/network/processed) based on raw data in [data/network/raw](./data/network/raw).  
Then, [NetworkOutput.R](./scripts/NetworkOutput.R) outputs results in [img/network](./img/network).  

#### Phylogentic analysis
run [Phylogeny.R](./scripts/Phylogeny.R). results are stored in [img/phylogeny](./img/phylogeny). 
  
See the Method and Supplementary Materials for more details.


## Table of Contents
* [README](./README.md) - this file
* [data](./data)
 	* [network](./data/network)
 		* [raw](./data/network/raw)
 			* [diacamma_indicum](./data/network/raw/diacamma_indicum) - Data from [Kolay and Annagiri 2015](https://doi.org/10.1098/rsos.150104)
 			* [diacamma_sp](./data/network/raw/diacamma_sp) - Data for Diacamma cf indicum from Japan
 			* [temnothorax_albipennis](./data/network/raw/temnothorax_albipennis) - Data from [Richardson et al.  2018](https://doi.org/10.1098/rspb.2017.2726)
 			* [temnothorax_nylanderi](./data/network/raw/temnothorax_nylanderi) - Data from [Richardson et al.  2021](https://doi.org/10.1038/s42003-021-02048-7)
 			* [temnothorax_rugatulus](./data/network/raw/temnothorax_rugatulus) - Data from [Valentini et al.   2020](https://doi.org/10.1098/rspb.2019.2950)
 		* [processed](./data/network/processed)
 	* [trajectory](./data/trajectory)
 		* [raw](./data/trajectory/raw) - Raw data for each species tandem. Data from [Valentini et al. 2020](https://doi.org/10.7554/eLife.55395) and [Mizumoto and Dobata 2019](https://doi.org/10.1126/sciadv.aau6108), or this study [Diacamma-csv](./data/trajectory/raw/Diacamma-csv)
 		* [processed](./data/trajectory/processed)
 		* [Scale.xlsx](./data/trajectory/Scale.xlsx) - Scale information for D. cf indicum experiments
 		* [Tracking_failed.xlsx](./data/trajectory/Tracking_failed.xlsx) - Information for detection faulure during observation for D. cf indicum experiments
 	* [phylogeny](./data/phylogeny) -  Include tree data and csv data for ant recruitment from [Reeves and Moreau 2019](https://doi.org/10.26049/ASP77-2-2019-10)
* [scripts](./scripts)
	* [NetworkProcessing.R](./scripts/Output.R)
	* [TrajectoryProcessing.R](./scripts/Phylogeny.R)
	* [NetworkOutput.R](./scripts/Preprocess.R)
	* [TrajectoryOutputs.R](./scripts/Sources.R)
	* [Phylogeny.R](./scripts/Sources.R)
* [img](./img)
* [tracking](./tracking)
  * [TrackingConverter.R](./tracking/TrackingConverter.R) - Scripts to convert the results of FastTrack into UMATracker
  * [TrimFocusedTandem.py](./tracking/TrimFocusedTandem.py) - python scripts to create masked videos
  * [test](./tracking/test) - example of tracking
    * [TandemTrackingMethods.md](./tracking/test/TandemTrackingMethods.md) - a notebook for tracking
    * 1. [FastTrack](./tracking/test/FastTrack)
    * 2. [TrackingConverter.R](./tracking/TrackingConverter.R) with [data](./tracking/test/data)
    * 3. [VideoPrepPython](./tracking/test/VideoPrepPython)
    * 4. [UMATracker](./tracking/test/UMATracker)
