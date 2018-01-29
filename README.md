# MixtureAnalysis

This repo will be a collection of tools created to assist with the mixture analysis from WGBS bam files generated by bismark

Documentation will be added once the package is more complete. This repo is a WORK IN PROGRESS.

### Requirements
* Python3
* All dependencies listed in ```requirements.txt```
* All of [pysam's dependencies](http://pysam.readthedocs.io/en/latest/#)

## Temporary usage directions
### Install
* Clone or download the master branch to your local machine
* __(Optional, but HIGHLY recommended)__ Create a new python virtual environment and activate that virualenv
* Install all python requirements using ```pip3 install -r /path/to/requirements.txt```

### Plot read-level matrices
* Create a text file containing the bins you wish to plot (will be refered to below as __bins.txt__ below). One bin per line in this format: __chr19_55555__. This will generate a plot for the genomic region chr19:55455-55555
* To view all program options run ```python3 /path/to/PlotBinLevelReads.py -h```
* To generate a clustered heatmaps of individual bins in the file you created run:
```bash
python3 /path/to/PlotBinLevelReads.py -c -o /path/to/desired/output /path/to/bins.txt /path/to/alignment.bam
```
* If you are inputting a very large list of bins to plot you can speed up the processing by using multiple CPU cores. Specificy the number of CPU cores you wish to use using the -n flag. For example to use 24 processors:
```bash
python3 /path/to/PlotBinLevelReads.py -c -n 24 -o /path/to/desired/output /path/to/bins.txt /path/to/alignment.bam
```
