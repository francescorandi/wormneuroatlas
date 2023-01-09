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
          'pdoc',
      ],
      packages=['wormneuroatlas'],
      package_data={'wormneuroatlas': 
                     ['data/SEE_README_FOR_CITATIONS',
                      'data/anatlas_neuron_positions.txt',
                      'data/funatlas.h5',
                      'data/cengen.h5',
                      'data/cengen_021821_conservative_threshold3.csv',
                      'data/cengen_021821_liberal_threshold1.csv',
                      'data/cengen_021821_medium_threshold2.csv',
                      'data/cengen_021821_stringent_threshold4.csv',
                      'data/c_elegans.PRJNA13758.WS286.geneIDs.txt',
                      'data/deorphanization_media_6.csv',
                      'data/froonikcx_peptide_gpcr.txt',
                      'data/cached_seq_ids.txt',
                      'data/aconnectome.json',
                      'data/neuron_ids.txt',
                      'data/aconnectome_white_1986_L4.csv',
                      'data/aconnectome_white_1986_A.csv',
                      'data/aconnectome_white_1986_whole.csv',
                      'data/aconnectome_witvliet_2020_7.csv',
                      'data/aconnectome_witvliet_2020_8.csv',
                      'data/aconnectome_ids_ganglia.json',
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
