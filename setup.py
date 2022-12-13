#!/usr/bin/env python

from distutils.core import setup, Extension
from distutils.command.build_ext import build_ext

setup(name='wormneuroatlas',
      version=0.1,
      description='Many neural atlases in the same place '+\
                  '(see README for citations)',
      author='Francesco Randi',
      author_email='francesco.randi@gmail.com',
      install_requires=[
          'numpy',
          'pycurl',
          'certifi',
      ],
      packages=['wormneuroatlas'],
      package_data={'wormneuroatlas': 
                     ['data/anatlas_neuron_positions.txt',
                      'data/aconnectome.json',
                      'data/aconnectome_ids.txt',
                      'data/aconnectome_white_1986_L4.csv',
                      'data/aconnectome_white_1986_A.csv',
                      'data/aconnectome_white_1986_whole.csv',
                      'data/aconnectome_witvliet_2020_7.csv',
                      'data/aconnectome_witvliet_2020_8.csv',
                      'data/aconnectome_ids_ganglia.json',
                      'data/sensoryintermotor_ids.json',
                      'data/esconnectome_monoamines_Bentley_2016.csv',
                      'data/esconnectome_neuropeptides_Bentley_2016.csv',
                      'data/GenesExpressing-unc-7-unc-9-inx-_-eat-5-thrs2.csv',
                      'data/GenesExpressing-neuropeptides.csv',
                      'data/GenesExpressing-neuropeptide-receptors.csv',
                      'data/GenesExpressing-daf-2-thrs2.csv',
                      'data/GenesExpressing-npr-4-thrs2.csv',
                      'data/GenesExpressing-npr-11-thrs2.csv',
                      'data/GenesExpressing-pdfr-1-thrs2.csv']}
     )
