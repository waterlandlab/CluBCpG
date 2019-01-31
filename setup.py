from setuptools import setup

setup(name="MixtureAnalysis",
      version="0.1.6",
      description="Package to identify epialleles using read clustering from WGBS data",
      author="Anthony Scott, PhD",
      author_email="charles.scott@bcm.edu",
      license='MIT',
      packages=['MixtureAnalysis'],

      install_requires=['pysam', 'numpy', 'matplotlib', 'scikit-learn',
                        'seaborn', 'scipy', 'pandas', 'fastcluster', 'pebble'],

      scripts=['bin/CalculateCompleteBins',
               'bin/CombineClusterCompare',
               'bin/TrainModels',
               'bin/ImputeFromModels',
               'bin/ImputeCluster'],

      test_suite="MixtureAnalysis.tests.test_ParseBam",
      )
