from setuptools import setup
from clubcpg.__init__ import __version__, __author__

setup(name="clubcpg",
      version=__version__,
      description="Package to identify epialleles using read clustering from WGBS data", # todo update this
      author=__author__,
      author_email="charles.scott@bcm.edu",
      license='MIT',
      packages=['clubcpg', 'clubcpg_prelim'],

      install_requires=[
          'pysam', 
          'numpy', 
          'matplotlib>3,<3.1', 
          'scikit-learn==0.21.2',
          'joblib',
          'seaborn', 
          'scipy', 
          'pandas', 
          'fastcluster', 
          'pebble', 
          'tqdm'],

      scripts=[
          'bin/clubcpg-coverage',
            'bin/clubcpg-cluster',
            'bin/clubcpg-impute-train',
            'bin/clubcpg-impute-coverage',
            'bin/clubcpg-impute-cluster',
      ]

      )
