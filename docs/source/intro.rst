================
Getting Started
================

Requirements
=============

* Python 3.5+ (Fully tested on Python 3.5, 3.6, and 3.7)
* Installation of samtools on your PATH
* Any \*nix operating system

    * (tested on Ubuntu, CentOS, RHEL)

.. NOTE::
    CluBCpG utilized Pysam, which is a wrapper around the htslib and samtools C-APIs, to read from BAM files. These tools
    do not support Windows, therefore, CluBCpG cannot be run on Windows.


Installation
=============

pip install
-------------

These steps will install CluBCpG.

.. NOTE::
    While optional, it is **highly recommended** to install CluBCpG in a virtual environment.

0. (Optional) Create and activate python virtualenv

.. code-block:: bash

    virtualenv ~/venv/clubcpg -p python3
    source activate ~/venv/clubcpg/bin/activate

1. Install CluBCpG with the following command

.. code-block:: bash

    pip install clubcpg


conda
------
Coming soon

Running Tests
==============

Included tests can be run to verify CluBCpG is correctly interacting with samtools to read form BAM files. The source code
can be downloaded from GitHub and the following can be run:

.. code-block:: bash

    python -m unittest -v clubcpg/tests/test_Module.py

Running this will download a small amount of test data in BAM format and will verify all dependencies are met
to correctly run CluBCpG
