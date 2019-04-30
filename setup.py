from setuptools import setup

setup(name="clubcpg",
      version="0.1.9",
      description="Package to identify epialleles using read clustering from WGBS data", # todo update this
      author="Anthony Scott, PhD",
      author_email="charles.scott@bcm.edu",
      license='MIT',
      packages=['clubcpg', 'clubcpg_prelim'],

      install_requires=[
          'pysam', 
          'numpy', 
          'matplotlib>3,<3.1', 
          'scikit-learn', 
          'joblib',
          'seaborn', 
          'scipy', 
          'pandas', 
          'fastcluster', 
          'pebble', 
          'tqdm'],

      scripts=['bin/CalculateBinCoverage',
               'bin/ClusterReads',
               'bin/TrainModels',
               'bin/ImputeBinCoverage',
               'bin/ImputeClusterReads'],

      # test_suite="clubcpg.tests.test_Module.py",
      )
