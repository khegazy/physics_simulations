import sys, os, glob
import logging
import argparse
import pickle as pl
import numpy as np
import scipy as sp
import h5py
import matplotlib.pyplot as plt
#import spherical_functions as sf
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
from scipy.special import factorial, jacobi
from multiprocessing import Pool
from functools import partial

sys.path.append("../../")
from parameters import get_parameters
from modules.script_setup import *
from modules.ADM import *

parser = argparse.ArgumentParser()

parser.add_argument(
    '--molName', type=str, default="NULL",
    help='Name of the molecule.'
)
parser.add_argument(
    '--run', type=str, default="NULL",
    help='Run name.'
)
parser.add_argument(
    '--calculation_type', type=str, default="analytic",
    help='Use parameters.py'
)
parser.add_argument(
    '--get_parameters', type=bool, default=True,
    help='Use parameters.py'
)
parser.add_argument(
    '--xyz_file', type=str, default="NULL",
    help='Address to XYZ file.'
)
parser.add_argument(
    '--detector_width', default=0.04, type=float,
    help='Detector width in meters.'
)
parser.add_argument(
    '--detector_height', default=0.04, type=float,
    help='Detector height in meters.'
)
parser.add_argument(
    '--imgBins', default=29, type=int,
    help='Number of height bins.'
)
parser.add_argument(
    '--beamline_length', default=4, type=float,
    help='Distance between the sample and the detector.'
)
parser.add_argument(
    '--smear_polar_bins', default=30, type=int,
    help='Number of polar bin sampling when smearing distribution.'
)
parser.add_argument(
    '--smear_azim_bins', default=30, type=int,
    help='Number of azimuthal bin sampling when smearing distribution.'
)
parser.add_argument(
    '--basis_folder', default="./bases/UED_best_fit_Nbases-4/", type=str,
    help='Folder containing the time bases'
)
parser.add_argument(
    '--basis_params', default=None, type=str,
    help='File name containing the fit results of the time bases'
)
parser.add_argument(
    '--basis_eval_params', default=None, nargs='+', type=int,
    help='Eval times [start, end, spacing]'
) 
parser.add_argument(
    '--get_bases_type', default=0, type=int,
    help='Specify which routine to use in get_bases.'
)
parser.add_argument(
    '--probe_FWHM', default=None, type=float,
    help='Electron pulse width (FWHM) in fs.'
)
parser.add_argument(
    '--scat_amps_dir', default="../scatteringAmplitudes/3.7MeV/", type=str,
    help='Path to folder containing the scattering amplitudes'
)
parser.add_argument(
    '--elEnergy', default=3.7e6, type=float,
    help='Electron energy in the e beam'
)

atom_names = {
  "H" : "hydrogen",
  "C" : "carbon",
  "O" : "oxygen",
  "N" : "nitrogen",
  "I" : "iodine",
  "S" : "sulfur"
}

"""
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
"""

def get_scattering_amplitudes(FLAGS, atom_types):

  scat_amps = {}
  dBrog = calculate_deBroglie_wavelength(FLAGS)
  print("DBW", dBrog)
  for atm in atom_types:
    if atm in scat_amps:
      continue

    angStr = []
    sctStr = []
    fName = os.path.join(FLAGS.scat_amps_dir, atom_names[atm] + "_dcs.dat")
    print(fName)
    with open(fName, 'r') as inpFile:
      ind=0
      for line in inpFile:
        if ind < 31:
          ind += 1
          continue

        angStr.append(line[2:11])
        sctStr.append(line[39:50])   
    
    angs = np.array(angStr).astype(np.float64)*np.pi/180
    q = 4*np.pi*np.sin(angs/2.)/dBrog
    scts = np.sqrt(np.array(sctStr).astype(np.float64))

    scat_amps[atm] = interp1d(q, scts, 'cubic')

  return scat_amps


def calculate_deBroglie_wavelength(FLAGS):
  C_AU = 1./0.0072973525664
  eV_to_au = 0.0367493
  angs_to_au = 1e-10/5.291772108e-11
  db_lambda = 2*np.pi*C_AU/\
    np.sqrt((FLAGS.elEnergy*eV_to_au + C_AU**2)**2\
    - (C_AU)**4) #au
  db_lambda /= angs_to_au
  
  return db_lambda
  #return 0.00296622 #angs


def get_molecule_distances(FLAGS, logging):
  fName = os.path.join(FLAGS.xyzDir, FLAGS.xyz_file)
  if not os.path.exists(fName):
    logging.fatal("Cannot find xyz file: " + fName)
    sys.exit(0)
  
  atom_types      = []
  atom_positions  = []
  with open(fName) as file:
    for i,ln in enumerate(file):
      if i == 0:
        Natoms = int(ln)
      elif i > 1:
        vals = ln.split()
        atom_types.append(vals[0])
        pos = [float(x) for x in vals[1:]]
        atom_positions.append(np.array([pos]))

  atom_positions = np.concatenate(atom_positions, axis=0)

  atom_distances_xyz    = []
  atom_distances_polar  = []
  atom_distances_types  = []
  for i,d1 in enumerate(atom_positions):
    for j_,d2 in enumerate(atom_positions[i+1:]):
      j = j_ + i + 1
      diff = d1-d2
      atom_distances_xyz.append(np.array([diff]))
      dist  = np.linalg.norm(diff)
      polar = np.arccos(diff[2]/dist)
      azim  = np.arctan2(diff[1], diff[0])
      atom_distances_polar.append(np.array([[dist, polar, azim]]))
      atom_distances_types.append([atom_types[i], atom_types[j]])
  atom_distances_xyz = np.concatenate(atom_distances_xyz, axis=0)
  atom_distances_polar = np.concatenate(atom_distances_polar, axis=0)
  
  return atom_types, atom_positions,\
      atom_distances_xyz, atom_distances_polar, atom_distances_types


def get_wignerD_comparison(phi_ea, theta_ea, chi_ea, l, m):
  lmm = np.ones((2*l+1, 3))*l
  lmm[:,1] -= np.arange(2*l+1)
  lmm[:,2] = m
  lmm = lmm.astype(int)

  D_ME = sf.Wigner_D_element(phi_ea, theta_ea, chi_ea, lmm)
  return D_ME, lmm


def small_factorial(start, end):
  start = np.max(start, 1)
  return np.prod(np.arange(start, end+1))

def nCk(n, k):
  return small_factorial(k,n)/small_factorial(1, n-k+1)

def nCk_np(n, k):
  return factorial(n)/(factorial(k)*factorial(n-k))

def Ylm_calc(m, l, azim, polar):
  print(m.shape, sp.special.sph_harm(m, l, azim, polar).shape)
  return (-1)**np.abs(m)*sp.special.sph_harm(m, l, azim, polar)

def get_wignerD_3d(phi_ea, theta_ea, chi_ea, l, m):
  # https://en.wikipedia.org/wiki/Wigner_D-matrix#Wigner_(small)_d-matrix

  phi_ea    = np.expand_dims(phi_ea, axis=-1)
  theta_ea  = np.expand_dims(theta_ea, axis=-1)
  chi_ea    = np.expand_dims(chi_ea, axis=-1)
  m_sum = np.ones((1, 2*l+1))
  m_sum[0,:] = np.ones(2*l+1)*l - np.arange(2*l+1)
  m_sum = m_sum.astype(int)

  # Calculate k, a, b, and lambda for jacobi polynomial
  k_array = np.array([l+m, l-m, l, l], dtype=int)
  k_matrix = np.tile(k_array, (2*l+1, 1)).transpose()

  k_matrix[2,:] += m_sum[0]
  k_matrix[3,:] -= m_sum[0]

  a_array = np.array([-1*m, m, m, -1*m])
  a_matrix = np.tile(a_array, (2*l+1, 1)).transpose()
  a_matrix += np.array([[1,-1,-1,1]]).transpose()*m_sum

  k_inds = np.argmin(k_matrix, axis=0)
  k, a = np.zeros_like(m_sum), np.zeros_like(m_sum)
  ii = np.arange(2*l+1)
  k[0,:] = k_matrix[k_inds[ii],ii]
  a[0,:] = a_matrix[k_inds[ii],ii]
  b = 2*l - 2*k - a
  lmbd = (m_sum - m)*np.maximum(0, (1-k_inds%3))

  # Calculate little d
  d_coeff = (-1)**lmbd*np.sqrt(nCk_np(2*l-k, k+a))/np.sqrt(nCk_np(k+b, b))\
      *np.sin(theta_ea/2)**a*np.cos(theta_ea/2)**b
  d = np.zeros_like(d_coeff)
  cos_theta = np.cos(theta_ea)
  for im in range(m_sum.shape[-1]):
    p = jacobi(k[0,im], a[0,im], b[0,im])
    d[:,im] = d_coeff[:,im]*np.polyval(p, cos_theta)[:,0]

  # Calculate D matrix element
  # In complex exp -i -> i to agree with sf package and agrees with active 
  #    rotation test
  D_ME = np.exp(np.complex(0,1)*m_sum*phi_ea)*d*np.exp(np.complex(0,1)*m*chi_ea)

  return D_ME.transpose(), l, m_sum[0,:], m


def rotate_Ylm_3d(phi_ea, theta_ea, chi_ea, l, m, polar, azim):
  D_ME, l, m_sum, m = get_wignerD_3d(phi_ea, theta_ea, chi_ea, l, m)

  Ylms = Ylm_calc(
    np.expand_dims(np.expand_dims(m_sum, -1), -1),
    np.ones((m_sum.shape[0], 1, 1))*l,
    np.expand_dims(azim,0),        
    np.expand_dims(polar,0))

  return np.sum(
      np.expand_dims(np.expand_dims(D_ME, -1), -1)\
      *np.expand_dims(Ylms, 1), axis=0)


def make_atomic_diffraction(
    detector_shape, detector_dist, s_map_polar,
    atom_types, scat_amps):

  diffraction = np.zeros((int(detector_shape[0]), int(detector_shape[1])))
  for atm in atom_types:
    diffraction += (1/detector_dist**2)*scat_amps[atm](s_map_polar[:,:,0])**2
  
  return diffraction


def make_diffraction(
    detector_shape, detector_dist, s_map_polar,
    atom_types, atom_distances_polar,
    smearing_ea, smearing_weights,
    scat_amps):

  smearing_weights = np.expand_dims(np.expand_dims(smearing_weights, -1), -1)
  diffraction = np.zeros((int(detector_shape[0]), int(detector_shape[1])))
  for ir, r_diff in enumerate(atom_distances_polar):

    # Rotate Y10 so it points in direction of r in lf
    #D_ME_, lmm = get_wignerD(r_diff[2], r_diff[1], 0, 1, 0)
    D_ME, l, m_sum, m = get_wignerD_3d(
        np.array([r_diff[2]]), np.array([r_diff[1]]), np.zeros(1), 1, 0)
    D_ME = D_ME[:,0]

    for im, m in enumerate(m_sum):
      if np.abs(D_ME[im]) <= 1e-10:
        #print("skipping")
        continue

      s_dot_r_lf_ = D_ME[im]*(s_map_polar[:,:,0]*r_diff[0]/(0.5*np.sqrt(3/np.pi)))\
          *rotate_Ylm_3d(smearing_ea[:,0], smearing_ea[:,1], smearing_ea[:,2],\
            l, m, s_map_polar[:,:,1], s_map_polar[:,:,2])
      s_dot_r_lf  = np.real(s_dot_r_lf_)
      smearing = np.sum(smearing_weights*np.cos(s_dot_r_lf), axis=0)
      diffraction += smearing*scat_amps[atom_types[ir][0]](s_map_polar[:,:,0])\
          *scat_amps[atom_types[ir][1]](s_map_polar[:,:,0])/(detector_dist**2)
  
  return diffraction


def numerical_alignment_diffraction(
    itm, times, bases,
    detector_shape, detector_dist, s_map_polar,
    atom_types, atom_distances_polar,
    smearing_ea, smearing_jacobian, scat_amps):

  sh_l = np.arange(bases.shape[0])*2
  sh_m = np.zeros_like(sh_l)
  smearing_weights = np.real(np.sum(np.expand_dims(bases[:,itm], axis=-1)\
      *sp.special.sph_harm(
        np.expand_dims(sh_m, axis=-1), np.expand_dims(sh_l, axis=-1),\
        np.expand_dims(smearing_ea[:,1], axis=0), np.expand_dims(smearing_ea[:,0], axis=0)),
      axis=0)).astype(np.float64)
  smearing_weights *= smearing_jacobian

  diffraction = make_diffraction(
      detector_shape, detector_dist, s_map_polar,
      atom_types, atom_distances_polar,
      smearing_ea, smearing_weights, scat_amps)

  # Saving Results
  fName = os.path.join("output",
      "mol_diffraction_Q-{0:.5g}_time-{1:.3g}.npy".format(
        s_map_polar[s_map_polar.shape[0]//2,-1,0], times[itm]))
  with open(fName, "wb") as file:
    np.save(file, diffraction)
  
  fName = os.path.join("output",
      "mol_diffraction_Q-{0:.5g}_time-{1:.3g}_shape[{2},{3}].dat".format(
        s_map_polar[s_map_polar.shape[0]//2,-1, 0], times[itm],
        diffraction.shape[0], diffraction.shape[1]))
  with open(fName, "wb") as file:
    diffraction.astype(np.double).tofile(file)



def main(FLAGS):


  #############################################
  #####  Lab Frame Coordinate System      #####
  #####                                   #####
  #####  Z: along the laser polarization  ##### 
  #####  Y: electron beam direction       #####
  #####  X: set by Z, Y                   #####
  #############################################


  atom_types, atom_positions,\
  atom_distances_xyz, atom_distances_polar,\
  atom_distances_types = get_molecule_distances(FLAGS, logging)
  scat_amps = get_scattering_amplitudes(FLAGS, atom_types)

  db_wvl  = calculate_deBroglie_wavelength(FLAGS) #angs
  k0      = 2*np.pi/db_wvl

  x_lab,z_lab = np.meshgrid(
      np.linspace(-1*FLAGS.detector_width/2, FLAGS.detector_width/2, FLAGS.imgBins),
      np.linspace(-1*FLAGS.detector_height/2, FLAGS.detector_height/2, FLAGS.imgBins))
  pixel_dist = np.sqrt(x_lab**2 + z_lab**2)
  detector_dist = np.sqrt(pixel_dist**2 + FLAGS.beamline_length**2)

  if FLAGS.calculation_type == "numerical":
  
    LMK, bases, norms, times  = get_bases(FLAGS, logging, plot=True)
    params, bases, times      = apply_fit_parameters(FLAGS, bases, times, logging)

    # Set up pixel and s maps
    x,z = x_lab,z_lab
    zero_mask = pixel_dist > 0

    signal_xyz_qf = np.concatenate([
        np.expand_dims(x, axis=-1),
        np.expand_dims(np.ones_like(x)*FLAGS.beamline_length, axis=-1),
        np.expand_dims(z, axis=-1)],
        axis=-1)

    incoming_e      = np.zeros((1,1,3))
    incoming_e[:,:,1] = k0
    s_xyz_lf    = k0*signal_xyz_qf/np.expand_dims(np.linalg.norm(signal_xyz_qf, axis=-1), axis=-1) - incoming_e
    s_map       = np.linalg.norm(s_xyz_lf, axis=-1)

    s_map = np.sqrt(x**2 + z**2)
    s_theta_lf  = np.zeros_like(z)
    s_theta_lf[zero_mask]  = np.arccos(s_xyz_lf[zero_mask,2]/s_map[zero_mask])
    s_phi_lf  = np.zeros_like(z)
    s_phi_lf[zero_mask]  = np.arctan2(s_xyz_lf[zero_mask,1], s_xyz_lf[zero_mask,0])

    s_map_polar = np.concatenate([
        np.expand_dims(s_map, axis=-1),
        np.expand_dims(s_theta_lf, axis=-1),
        np.expand_dims(s_phi_lf, axis=-1)],
        axis=-1)

    atm_diffraction = make_atomic_diffraction(
        (FLAGS.imgBins, FLAGS.imgBins), detector_dist,
        s_map_polar, atom_types, scat_amps)

    # Saving Results
    fName = os.path.join("output",
        "atm_diffraction_Q-{0:.5g}.npy".format(s_map[s_map.shape[0]//2,-1]))
    with open(fName, "wb") as file:
      np.save(file, atm_diffraction)
    
    print(atm_diffraction.shape)
    fName = os.path.join("output",
        "atm_diffraction_Q-{0:.5g}_shape[{1},{2}].dat".format(
          s_map[s_map.shape[0]//2,-1],
          atm_diffraction.shape[0], atm_diffraction.shape[1]))
    with open(fName, "wb") as file:
      atm_diffraction.astype(np.double).tofile(file)


    # Ensemble smearing
    smearing_polar  = np.arange(FLAGS.smear_polar_bins)*np.pi/FLAGS.smear_polar_bins
    smearing_azim   = np.arange(FLAGS.smear_azim_bins)*2*np.pi/FLAGS.smear_azim_bins
    smearing_spin   = np.arange(FLAGS.smear_azim_bins)*2*np.pi/FLAGS.smear_azim_bins
    azim_inds       = np.arange(len(smearing_azim)).astype(int)
    spin_inds       = np.arange(len(smearing_spin)).astype(int)
    smearing_ea     = np.zeros(
        (smearing_polar.shape[0], smearing_azim.shape[0], smearing_spin.shape[0], 3))

   
    for i in range(FLAGS.smear_polar_bins):
      smearing_ea[i,:,:,0]  = smearing_azim[i]
      for k in range(FLAGS.smear_azim_bins):
        smearing_ea[:,i,k,1] = smearing_polar[i]
        smearing_ea[i,azim_inds,:,2]  = smearing_spin[azim_inds]
    smearing_ea       = np.reshape(smearing_ea, (-1,3))
    smearing_jacobian = np.sin(smearing_ea[:,1])

    for tm in range(bases.shape[-1]):
      
      numerical_alignment_diffraction(
          tm, times, bases,
          (FLAGS.imgBins, FLAGS.imgBins),
          detector_dist, s_map_polar,
          atom_distances_types, atom_distances_polar,
          smearing_ea, smearing_jacobian, scat_amps)

  elif FLAGS.calculation_type == "analytic":
    if FLAGS.QperPix is None:
      logging.fatal("Must specify QperPix for the analytic option.")
      sys.exit(0)

    LMK, bases, norms, times  = get_bases(FLAGS, logging, plot=True)
    params, bases, times      = apply_fit_parameters(FLAGS, bases, times, logging)


    x,z = np.meshgrid(
        np.linspace(
          -1*FLAGS.QperPix*FLAGS.detector_width/2,
          FLAGS.QperPix*FLAGS.detector_width/2, FLAGS.imgBins),
        np.linspace(
          -1*FLAGS.QperPix*FLAGS.detector_height/2,
          FLAGS.QperPix*FLAGS.detector_height/2, FLAGS.imgBins))

    q = np.sqrt(x**2 + z**2)


  elif FLAGS.calculation_type == "azmAvg":
    q = np.arange(FLAGS.NradAzmBins)*FLAGS.QperPix

    atm_diffraction = np.zeros_like(q)
    mol_diffraction = np.zeros_like(q)
    for i,atm in enumerate(atom_types):
      atm_diffraction += scat_amps[atm](q)**2
     
    for i, atmTypes in enumerate(atom_distances_types):
      print(atmTypes[0], atmTypes[1], atom_distances_polar[i,0])
      mol_diffraction += np.sinc(q*atom_distances_polar[i,0])\
          *scat_amps[atmTypes[0]](q)*scat_amps[atmTypes[1]](q)


    fName_template_dat =\
        "{0}_sim_{1}Diffraction-lineout_align-random_Qmax-{2:.4g}_Bins[{3}].dat"
    #fName_template_npy =
    #    "{0}_sim_{1}Diffraction-lineout_align-random_Qmax-{2:.4g}.npy"
    fName = os.path.join(FLAGS.simOutputDir, fName_template_dat.format(
        FLAGS.molName, "atm", q[-1], FLAGS.NradAzmBins))
    with open(fName, "wb") as file:
      atm_diffraction.astype(np.double).tofile(file)
    """
    fName = os.path.join(FLAGS.simOutputDir, fName_template_npy.format(
        FLAGS.molName, "atm", q[-1]))
    with open(fName, "wb") as file:
      np.save(file, atm_diffraction)
    """

    fName = os.path.join(FLAGS.simOutputDir, fName_template_dat.format(
        FLAGS.molName, "mol", q[-1], FLAGS.NradAzmBins))
    with open(fName, "wb") as file:
      mol_diffraction.astype(np.double).tofile(file)
    """
    fName = os.path.join(FLAGS.simOutputDir, fName_template_npy.format(
        FLAGS.molName, "mol", q[-1]))
    with open(fName, "wb") as file:
      np.save(file, mol_diffraction)
    """

    fName = os.path.join(FLAGS.simOutputDir, fName_template_dat.format(
        FLAGS.molName, "sms", q[-1], FLAGS.NradAzmBins))
    with open(fName, "wb") as file:
      (q*mol_diffraction/atm_diffraction).astype(np.double).tofile(file)
    """
    fName = os.path.join(FLAGS.simOutputDir, fName_template_npy.format(
        FLAGS.molName, "sms", q[-1]))
    with open(fName, "wb") as file:
      np.save(file, q*mol_diffraction/atm_diffraction)
    """


    fName_template_h5 =\
        "{0}_sim_diffraction-azmAvg_align-random_Qmax-{1:.4g}.h5"
    fName = os.path.join(FLAGS.simOutputDir, fName_template_h5.format(
        FLAGS.molName, q[-1]))
    with h5py.File(fName, 'w') as hf:#open(fName, "wb") as file:
      hf.create_dataset("atm_diffraction",  data=atm_diffraction)
      hf.create_dataset("mol_diffraction",  data=mol_diffraction)
      hf.create_dataset("sms_diffraction",  data=q*mol_diffraction/atm_diffraction)


    # Plot Results
    fig, ax = plt.subplots(1, 2, figsize=(12,5))
    ax[0].plot(q, atm_diffraction)
    ax[0].plot(q, mol_diffraction)
    ax[1].plot(q, q*mol_diffraction/atm_diffraction)
    ax[0].set_xlim([0, q[-1]])
    ax[1].set_xlim([0, q[-1]])
    fig.savefig(os.path.join("./plots",
        "{}_sim-diffraction-azmAvg.png".format(FLAGS.molName)))


if __name__ == '__main__':
  FLAGS = setup(parser)
  if FLAGS.get_parameters:
    calc_type = FLAGS.calculation_type
    FLAGS = get_parameters(FLAGS.run)
    FLAGS.calculation_type = calc_type

  main(FLAGS)
  

