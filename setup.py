from setuptools import setup

setup(name="clubcpg",
      version="0.1.13",
      description="Package to identify epialleles using read clustering from WGBS data", # todo update this
      author="C. Anthony Scott, PhD",
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
