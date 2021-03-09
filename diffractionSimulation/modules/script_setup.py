import sys, os, glob
import logging
import argparse
import numpy as np


def setup(parser, output_dir="./"):
  FLAGS = parser.parse_args()

  if FLAGS.LMK is not None:
    lmk = FLAGS.LMK.split(",")
    FLAGS.LMK = np.reshape(np.array(lmk).astype(int), (-1, 3))
    print("H0000", FLAGS.LMK)

  if not os.path.exists(os.path.join(output_dir, "logs")):
    os.makedirs(os.path.join(output_dir, "logs"))

  logging.basicConfig(
      level=logging.DEBUG,
      format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
      datefmt='%m-%d %H:%M',
      filename=os.path.join(output_dir, 'logs/diffraction.log'),
      filemode='w')
  console = logging.StreamHandler()
  console.setLevel(logging.INFO)
  formatter = logging.Formatter('%(levelname)-5s %(message)s')
  console.setFormatter(formatter)
  logging.getLogger('').addHandler(console)

  return FLAGS

