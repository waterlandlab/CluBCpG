================
Getting Started
================

.. contents:: Table of Contents

Requirements
=============

* Python 3.5+ (Fully tested on Python 3.5, 3.6, and 3.7)
* Installation of samtools on your PATH
* nix operating system (Ubuntu, CentOS, MacOS (presumably, but not tested on Mac)

.. NOTE::
    CluBCpG utilized Pysam, which is a wrapper around the htslib and samtools C-APIs, to read from BAM files. These tools
    do not support Windows, therefore, CluBCpG cannot be run on Windows.


Installation
=============

ClubCpG is not yet deposited onto PyPi, but can still be installed using pip. These steps will install CluBCpG.

.. NOTE::
    While optional, it is highly recommended to install CluBCpG in a virtual environment

1. Clone the public GitHub repo into a local directory

.. code-block:: bash

    git clone https://github.com/waterlandlab/CluBCpG.git

2. Navigate into the folder containing `setup.py` and run `pip install .`

.. code-block:: bash

    cd CluBCpg
    pip install .

Running Tests
==============

Included tests can be run to verify CluBCpG is correctly interacting with samtools to read form BAM files.

.. code-block:: bash

    python -m unittest -v clubcpg/tests/test_Module.py

Running this will download a small amount of test data in BAM format and will verify all dependencies are met
to correctly run CluBCpG
