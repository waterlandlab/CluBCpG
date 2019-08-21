from setuptools import setup

with open("README.md", "rt") as f:
    long_description = f.read()

setup(name="clubcpg",
      version="0.1.19",
      description="CluBCpG is a software package built to analyze whole genome bisulfite sequencing (WGBS) data",
      long_description=long_description,
      long_description_content_type="text/markdown",
      author="C. Anthony Scott, PhD",
      author_email="charles.scott@bcm.edu",
      url="https://github.com/waterlandlab/CluBCpG",
      license='MIT',
      packages=['clubcpg', 'clubcpg_prelim'],

      install_requires=[
          'pysam', 
          'numpy',
          'scikit-learn==0.21.2',
          'joblib',
          'scipy', 
          'pandas', 
          'fastcluster', 
          'pebble', 
          'tqdm'
      ],

      classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
        "Operating System :: POSIX",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
      ],

      scripts=[
          'bin/clubcpg-coverage',
            'bin/clubcpg-cluster',
            'bin/clubcpg-impute-train',
            'bin/clubcpg-impute-coverage',
            'bin/clubcpg-impute-cluster',
      ]

      )
