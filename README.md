# CluBCpG: Cluster-Based analysis of CpG methylation

[![Build Status](https://travis-ci.com/waterlandlab/CluBCpG.svg?branch=master)](https://travis-ci.com/waterlandlab/CluBCpG)
[![Documentation Status](https://readthedocs.org/projects/clubcpg/badge/?version=latest)](https://clubcpg.readthedocs.io/en/latest/?badge=latest)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/waterlandlab/CluBCpG.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/waterlandlab/CluBCpG/context:python)


## What is CluBCpG?
CluBCpG is a software package built to analyze whole genome bisulfite sequencing (WGBS) data. This toolkit will divide each chromosome into small user-defined intervals, extract all WGBS reads within those intervals, cluster them based on identity, and write a final report to the use containing all identified CpG methylation patterns.

## How do I use this?
Full documentation is available on [ReadTheDocs](https://clubcpg.readthedocs.io/en/latest/index.html)

### Requirements
* Python >= 3.5
* Installation of Samtools on your PATH

### Install
* Clone or download the master branch to your local machine
* __(Optional, but HIGHLY recommended)__ Create a new python virtual environment and activate that virualenv
* Navigate to folder containing setup.py
* Execute `pip install .` within that folder to install the package. Requirements will automatically be installed if not already present.

## Help! This isnt working.
Open an [Issue](https://github.com/waterlandlab/CluBCpG/issues/new/choose)

## Can you make it do this?
Open a [Feature Request](https://github.com/waterlandlab/CluBCpG/issues/new/choose)
