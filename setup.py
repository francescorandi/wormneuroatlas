#!/usr/bin/env python

from distutils.core import setup, Extension
from distutils.command.build_ext import build_ext

setup(name='wormneuroatlas',
      version=0.1,
      description='Many neural atlases in the same place',
      author='Francesco Randi',
      author_email='francesco.randi@gmail.com',
      packages=['pumpprobe'],
      package_data={'pumpprobe': ['anatlas_neuron_positions.txt',
                                  'aconnectome.json',
                                  'aconnectome_ids.txt',
                                  'aconnectome_white_1986_L4.csv',
                                  'aconnectome_white_1986_A.csv',
                                  'aconnectome_white_1986_whole.csv',
                                  'aconnectome_witvliet_2020_7.csv',
                                  'aconnectome_witvliet_2020_8.csv',
                                  'aconnectome_ids_ganglia.json',
                                  'sensoryintermotor_ids.json',
                                  'esconnectome_monoamines_Bentley_2016.csv',
                                  'esconnectome_neuropeptides_Bentley_2016.csv',
                                  'GenesExpressing-unc-7-unc-9-inx-_-eat-5-thrs2.csv',
                                  'GenesExpressing-neuropeptides.csv',
                                  'GenesExpressing-neuropeptide-receptors.csv',
                                  'GenesExpressing-daf-2-thrs2.csv',
                                  'GenesExpressing-npr-4-thrs2.csv',
                                  'GenesExpressing-npr-11-thrs2.csv',
                                  'GenesExpressing-pdfr-1-thrs2.csv']}
     )
