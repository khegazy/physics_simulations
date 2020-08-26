import sys, os, glob
import logging
import argparse


def setup(parser):
  FLAGS = parser.parse_args()

  logging.basicConfig(
      level=logging.DEBUG,
      format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
      datefmt='%m-%d %H:%M',
      filename='./logs/diffraction.log',
      filemode='w')
  console = logging.StreamHandler()
  console.setLevel(logging.INFO)
  formatter = logging.Formatter('%(levelname)-5s %(message)s')
  console.setFormatter(formatter)
  logging.getLogger('').addHandler(console)

  return FLAGS

