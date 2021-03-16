import sys, os, glob
import logging
import argparse
import pickle as pl
import numpy as np
import scipy as sp
import h5py
from copy import copy as copy
import matplotlib.pyplot as plt
#import spherical_functions as sf
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
from scipy.special import factorial, jacobi
from multiprocessing import Pool
from functools import partial

sys.path.append("./")
sys.path.append("../../")
from parameters import get_parameters
from modules.script_setup import *


parser = argparse.ArgumentParser()

parser.add_argument(
    '--molecule', type=str, default=None,
    help='Name of the molecule.'
)
parser.add_argument(
    '--run', type=str, default="simulation",
    help='Run name.'
)
parser.add_argument(
    '--eval_time', type=float, default=None,
    help='Time at which to simulate the diffraction pattern.'
)
parser.add_argument(
    '--time_ind', type=str, default=None,
    help='Time index at which to simulate the diffraction pattern.'
)
parser.add_argument(
    '--calculation_type', type=str, default=None,
    help='Use parameters.py'
)
parser.add_argument(
    '--xyz_file', type=str, default=None,
    help='Address to XYZ file.'
)
parser.add_argument(
    '--basis_folder', default=None, type=str,
    help='Folder containing the time bases'
)
parser.add_argument(
    '--output_fileName_suffix', default=None, type=str,
    help='Suffix for output file'
)
parser.add_argument(
    '--cluster', default=False, type=bool,
    help='Suffix for output file'
)
parser.add_argument(
    '--LMK', default=None, type=str,
    help='List of L0,M0,K0,...,Ln,Mn,Kn'
)

"""
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
"""

atom_names = {
  "H" : "hydrogen",
  "C" : "carbon",
  "O" : "oxygen",
  "N" : "nitrogen",
  "I" : "iodine",
  "F" : "flourine",
  "S" : "sulfur",
  "Br": "bromine"
}

"""
def setup(parser):
  params = parser.parse_args()

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

  return params
"""

def get_scattering_amplitudes(params, atom_types):

  scat_amps = {}
  dBrog = calculate_deBroglie_wavelength(params)
  print("DBW", dBrog)
  for imol in range(len(atom_types)):
    for atm in atom_types[imol]:
      if atm in scat_amps:
        continue

      angStr = []
      sctStr = []
      fName = os.path.join(params.scat_amps_dir, atom_names[atm] + "_dcs.dat")
      print(fName)
      with open(fName, 'r') as inpFile:
        ind=0
        for line in inpFile:
          if ind < 31:
            ind += 1
            continue

          angStr.append(line[2:11])
          sctStr.append(line[39:50])   
      
      angs = np.array(angStr).astype(np.float64)*np.pi/180.0
      q = 4*np.pi*np.sin(angs/2.)/dBrog
      scts = np.sqrt(np.array(sctStr).astype(np.float64))

      scat_amps[atm] = interp1d(q, scts, 'cubic')

  return scat_amps


def calculate_deBroglie_wavelength(params):
  C_AU = 1./0.0072973525664
  eV_to_au = 0.0367493
  angs_to_au = 1e-10/5.291772108e-11
  db_lambda = 2*np.pi*C_AU/\
    np.sqrt((params.elEnergy*eV_to_au + C_AU**2)**2\
    - (C_AU)**4) #au
  db_lambda /= angs_to_au
  
  return db_lambda
  #return 0.00296622 #angs


"""
def get_molecule_distances(params, logging):
  fName = os.path.join(params.xyzDir, params.xyz_file)
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
    for j,d2 in enumerate(atom_positions):
      if i == j:
        continue
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
"""

def get_molecule_distances(params, logging):
  
  atom_types      = []
  atom_positions  = []
  atom_distances_xyz    = []
  atom_distances_polar  = []
  atom_distances_types  = []

  if not isinstance(params.xyz_file, list):
    params.xyz_file = [params.xyz_file]

  for xyz_name in params.xyz_file:
    fName = os.path.join(params.xyzDir, xyz_name)
    print("Looking at file " + fName)
    if not os.path.exists(fName):
      logging.fatal("Cannot find xyz file: " + fName)
      sys.exit(0)
    
    atom_types.append([])
    atom_positions.append([])
    atom_distances_xyz.append([])
    atom_distances_polar.append([])
    atom_distances_types.append([])
   
    fill_atom_types = True
    with open(fName) as file:
      cur_atom = []
      for i,ln in enumerate(file):
        vals = ln.split()
        if len(vals) == 0:
          continue
        elif len(vals) == 1:
          Natoms = int(vals[0])
          if len(cur_atom) == Natoms:
            atom_positions[-1].append(np.array(cur_atom))
            cur_atom = []
            fill_atom_types = False
        else:
          if len(vals[0]) > 2:
            continue
          pos = [float(x) for x in vals[1:]]
          cur_atom.append(np.array(pos))
          if fill_atom_types:
            atom_types[-1].append(vals[0])
      atom_positions[-1].append(np.array(cur_atom))


    #atom_positions = np.concatenate(atom_positions, axis=0)

    if atom_positions[-1][0].shape[0] == 1:
      atom_distances_xyz[-1] = None
      atom_distances_polar[-1] = None
      atom_distances_types[-1] = None
    else:
      fill_types = True
      for iatm in range(len(atom_positions[-1])):
        cur_distances_xyz, cur_distances_polar = [], []
        for i,d1 in enumerate(atom_positions[-1][iatm]):
          for j,d2 in enumerate(atom_positions[-1][iatm]):
            if i == j:
              continue
            diff = d1-d2
            cur_distances_xyz.append(np.array(diff))
            dist  = np.linalg.norm(diff)
            polar = np.arccos(diff[2]/dist)
            azim  = np.arctan2(diff[1], diff[0])
            cur_distances_polar.append(np.array([dist, polar, azim]))
            if fill_types:
              atom_distances_types[-1].append([atom_types[-1][i], atom_types[-1][j]])
        atom_distances_xyz[-1].append(np.array(cur_distances_xyz))
        atom_distances_polar[-1].append(np.array(cur_distances_polar))
        fill_types = False

      atom_distances_xyz[-1] = np.array(atom_distances_xyz[-1])
      atom_distances_polar[-1] = np.array(atom_distances_polar[-1])
   
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
  return sp.special.sph_harm(m, l, azim, polar, dtype=np.complex64)
  #return np.abs(-1)**np.abs(m)*sp.special.sph_harm(m, l, azim, polar, dtype=np.complex64)

def get_wignerD_3d(phi_ea, theta_ea, chi_ea, l, m, k):
  # Return array of D matrices with L, M, K as l, [-l, l], m
  # Calculation based off of below derivation
  # https://en.wikipedia.org/wiki/Wigner_D-matrix#Wigner_(small)_d-matrix

  m = np.expand_dims(m, 0)
  phi_ea    = np.expand_dims(phi_ea, axis=-1)
  theta_ea  = np.expand_dims(theta_ea, axis=-1)
  chi_ea    = np.expand_dims(chi_ea, axis=-1)

  k_array = np.array([l+k, l-k, l, l], dtype=int)
  k_matrix = np.tile(np.expand_dims(k_array, -1), m.shape[1])

  k_matrix[2,:] += m[0]
  k_matrix[3,:] -= m[0]

  a_array = np.array([-1*k, k, k, -1*k])
  a_matrix = np.tile(np.expand_dims(a_array, -1), m.shape[1])
  a_matrix += np.array([[1,-1,-1,1]]).transpose()*m

  k_inds = np.argmin(k_matrix, axis=0)
  k_, a = np.zeros_like(m), np.zeros_like(m)
  ii = np.arange(k_inds.shape[0])
  k_[0,:] = k_matrix[k_inds[ii],ii]
  a[0,:] = a_matrix[k_inds[ii],ii]
  b = 2*l - 2*k_ - a
  lmbd = (m - k)*np.maximum(0, (1-(k_inds%3)))
  #lmbd = np.maximum(0, (m - k))
  #print(lmbd, (m-k), (1-(k_inds%3)), k_inds)
  #print("ll",lmbd)


  # Calculate little d
  d_coeff = (-1)**lmbd*np.sqrt(nCk_np(2*l-k_, k_+a))/np.sqrt(nCk_np(k_+b, b))\
      *np.sin(theta_ea/2)**a*np.cos(theta_ea/2)**b
  d = np.zeros_like(d_coeff)
  cos_theta = np.cos(theta_ea)
  jac = []
  for im in range(m.shape[-1]):
    p = jacobi(k_[0,im], a[0,im], b[0,im])
    jac.append(np.polyval(p, cos_theta)[:,0])
    d[:,im] = d_coeff[:,im]*np.polyval(p, cos_theta)[:,0]

  #print("d", d)
  #print("jac", jac)
  #print("l", lmbd)
  #print("mmmmmmmmmmm", m, np.complex(0,-1))
  #print("phi", phi_ea)
  #print("res",np.exp(np.complex(0,-1)*m*phi_ea))
  # Calculate D matrix element
  D_ME = np.exp(np.complex(0,-1)*m*phi_ea)*d*np.exp(np.complex(0,-1)*k*chi_ea)
  #D_ME = d*(np.cos(-1*m*phi_ea) + np.cos(-1*k*chi_ea)\
  #    + np.complex(0,1)*(np.sin(-1*m*phi_ea) + np.sin(-1*k*chi_ea)))
  #print("WTF", np.sin(-1*m*phi_ea))
  #print(-1*m*phi_ea)

  return D_ME.transpose().astype(np.complex64)

def rotate_Ylm_3d(phi_ea, theta_ea, chi_ea, l, m, polar, azim):
  m_sum = np.ones(2*l+1)*l - np.arange(2*l+1)
  m_sum = m_sum.astype(np.int16)
  #m_sum = np.ones((1, 2*l+1), dtype=np.int16)
  #m_sum[0,:] = np.ones(2*l+1)*l - np.arange(2*l+1, dtype=np.int16)
  D_ME = get_wignerD_3d(phi_ea, theta_ea, chi_ea, l, m_sum, m)

  """
  print("D SIZE", D_ME.shape)
  print(D_ME)
  """
  Ylms = Ylm_calc(
    np.expand_dims(np.expand_dims(m_sum, -1), -1),
    np.ones((m_sum.shape[0], 1, 1))*l,
    np.expand_dims(azim,0),        
    np.expand_dims(polar,0))

  return np.einsum('ia,ibc->abc', D_ME, Ylms)
  #return np.sum(
  #    np.expand_dims(np.expand_dims(D_ME, -1), -1)\
  #    *np.expand_dims(Ylms, 1), axis=0)


def make_atomic_diffraction(
    detector_shape, detector_dist, q_map_polar,
    atom_types, scat_amps):

  diffraction = np.zeros((int(detector_shape[0]), int(detector_shape[1])))
  for imol in range(len(atom_types)):
    for atm in atom_types[imol]:
      diffraction += (1/detector_dist**2)\
          *scat_amps[atm](q_map_polar[:,:,0])**2
  
  return diffraction


def make_diffraction(
    detector_shape, detector_dist, q_map_polar,
    atom_types, atom_distances_polar,
    smearing_ea, smearing_weights,
    scat_amps, logging=None):

  n = 4000
  N_points = smearing_weights.shape[0]
  """
  plt.plot(np.cos(smearing_ea[:,1]), smearing_weights)
  plt.savefig("testw.png")
  plt.close()
  """
  smearing_weights = np.expand_dims(np.expand_dims(smearing_weights, -1), -1)
  diffraction = np.zeros((int(detector_shape[0]), int(detector_shape[1])),
      dtype=np.complex)

  #####  Loop over smearing points  #####
  for imol in range(len(atom_distances_polar)):
    for igeo in range(len(atom_distances_polar[imol])):
      lc,hc = 0,np.min([n, N_points])
      while hc < N_points or lc == 0:
        if logging is not None:
          logging.info("Evaluating range {} - {} of {}".format(lc, hc, N_points))
        else:
          print("Evaluating range {} - {} of {}".format(lc, hc, N_points))

        ###  Loop over all pairwise distances  #####
        for ir, r_diff in enumerate(atom_distances_polar[imol][igeo]):
          """
          rr = 40
          cent = int(detector_shape[1]/2)
          print(cent, smearing_ea[lc:hc,1].shape)
          irr = (rr*np.cos(smearing_ea[lc:hc,1])).astype(int) + cent
          icc = (rr*np.cos(smearing_ea[lc:hc,0])).astype(int) + cent
          print("WEIGHT SH", smearing_weights.shape)
          diffraction[irr,icc] += smearing_weights[lc:hc,0,0]
          continue
          """

          # Rotate Y10 so it points in direction of r in lf
          #D_ME_, lmm = get_wignerD(r_diff[2], r_diff[1], 0, 1, 0)
          l = 1
          m_sum = np.ones(3) - np.flip(np.arange(3))
          m_sum = m_sum.astype(np.int16)
          D_ME = get_wignerD_3d(
              np.array([r_diff[2]]), np.array([r_diff[1]]), np.zeros(1),
              1, m_sum, 0)
          #    np.array([r_diff[2]]), np.array([r_diff[1]]), np.zeros(1), 1, 0)
          D_ME = D_ME[:,0]

          q_dot_r = np.zeros(
              (hc - lc, q_map_polar.shape[0], q_map_polar.shape[1]),
              dtype=np.complex64)
          for im, m in enumerate(m_sum):
            if logging is not None:
              logging.info("\tir / m: {} / {}".format(ir, m))
            else:
              print("\tir / m: {} / {}".format(ir, m))

            if np.abs(D_ME[im]) <= 1e-10:
              print("skipping")
              continue

            """
            rr,cc = 3,4
            Nb = 10
            tst = rotate_Ylm_3d(
              smearing_ea[:Nb,0], smearing_ea[:Nb,1], smearing_ea[:Nb,2],\
                  l, m, q_map_polar[:,:,1], q_map_polar[:,:,2])
            print("VAR", np.sum(np.abs(tst - np.mean(tst,axis=0,keepdims=True))**2))
            """

            q_dot_r +=  D_ME[im]*(q_map_polar[:,:,0]*r_diff[0]/(0.5*np.sqrt(3/np.pi)))\
                  *rotate_Ylm_3d(
                      smearing_ea[lc:hc,0], smearing_ea[lc:hc,1], smearing_ea[lc:hc,2],\
                      l, m, q_map_polar[:,:,1], q_map_polar[:,:,2])

          #print("QDR", np.amax(np.real(q_dot_r)), np.amax(np.imag(q_dot_r)))
          diffraction +=\
              scat_amps[atom_types[imol][ir][0]](q_map_polar[:,:,0])\
              *scat_amps[atom_types[imol][ir][1]](q_map_polar[:,:,0])/(detector_dist**2)\
              *np.sum(smearing_weights[lc:hc]*np.cos(np.real(q_dot_r)), axis=0)
          #diffraction +=\
          #    scat_amps[atom_types[ir][0]](q_map_polar[:,:,0])\
          #    *scat_amps[atom_types[ir][1]](q_map_polar[:,:,0])/(detector_dist**2)\
          #    *np.real(np.sum(smearing_weights[lc:hc]*np.cos(q_dot_r), axis=0))
          #    #*np.sum(smearing_weights[lc:hc]*np.cos(np.real(q_dot_r)), axis=0)
        lc = hc
        hc = np.min([hc+n, N_points])

        """
          diffraction +=\
              scat_amps[atom_types[ir][0]](q_map_polar[:,:,0])\
              *scat_amps[atom_types[ir][1]](q_map_polar[:,:,0])/(detector_dist**2)\
              *np.sum(smearing_weights*np.cos(
                np.real(D_ME[im]*(q_map_polar[:,:,0]*r_diff[0]/(0.5*np.sqrt(3/np.pi)))\
                *rotate_Ylm_3d(smearing_ea[:,0], smearing_ea[:,1], smearing_ea[:,2],\
                  l, m, q_map_polar[:,:,1], q_map_polar[:,:,2]))),
                axis=0)
        """
        """
          q_dot_r_lf = np.real(D_ME[im]*(q_map_polar[:,:,0]*r_diff[0]/(0.5*np.sqrt(3/np.pi)))\
              *rotate_Ylm_3d(smearing_ea[:,0], smearing_ea[:,1], smearing_ea[:,2],\
                l, m, q_map_polar[:,:,1], q_map_polar[:,:,2]))
          smearing = np.sum(smearing_weights*np.cos(q_dot_r_lf), axis=0)
          diffraction += smearing*scat_amps[atom_types[ir][0]](q_map_polar[:,:,0])\
              *scat_amps[atom_types[ir][1]](q_map_polar[:,:,0])/(detector_dist**2)
        """
  
  return np.real(diffraction), diffraction


def diffraction_calculation_numerical(
    molecule, itm, times, LMK, bases,
    detector_shape, detector_dist, q_map_polar,
    atom_types, atom_distances_polar,
    smearing_ea, smearing_jacobian, scat_amps, logging=None):

  smearing_weights = []
  for lmk in LMK:
    smearing_weights.append(
        ((2*lmk[0] + 1)/(8*np.pi**2))*get_wignerD_3d(
          smearing_ea[:,0],
          smearing_ea[:,1],
          smearing_ea[:,2],
          lmk[0], np.array([lmk[1]]), lmk[2]))
    #smearing_weights[-1] /= np.sqrt(4*np.pi/(2*lmk[0] + 1))

 
  """
  if molecule == "N2O":
    print("NORMING")
    smearing_weights = 2*np.pi**2*np.conj(smearing_weights)\
        /(4*np.pi/(2*lmk[0] + 1)
  elif molecule == "NO2":
    smearing_weights /= np.power(np.sqrt(4*np.pi/(2*lmk[0] + 1)), 1.5)
  """

  norm = np.sqrt(4*np.pi/(2*LMK[:,0] + 1))
  #smearing_weights = np.real(np.sum(
  smearing_weights = np.sum(
      np.expand_dims(norm*bases[:,itm], axis=-1)\
        *np.concatenate(smearing_weights, axis=0),
      axis=0)#.astype(np.float32)

  smearing_weights *= smearing_jacobian
  
  """
  print(smearing_ea.shape, smearing_weights.shape)
  for ang,w in zip(smearing_ea, smearing_weights):
    print(ang, w)
  """
 
  return make_diffraction(
      detector_shape, detector_dist, q_map_polar,
      atom_types, atom_distances_polar,
      smearing_ea, smearing_weights, scat_amps, logging)

def calc_dists(R):
  r     = np.expand_dims(R, 0) - np.expand_dims(R, 1)
  dR    = np.sqrt(np.sum(r**2, axis=-1))
  theta = np.arccos(r[:,:,2]/(dR + 1e-20))
  phi   = np.arctan2(r[:,:,1], r[:,:,0])
  mask = np.expand_dims((dR[:,:] > 0) , -1)

  return mask*np.concatenate([np.expand_dims(dR,-1),\
    np.expand_dims(theta, -1),\
    np.expand_dims(phi, -1)], axis=-1)


def diffraction_calculation_analytical(
    itm, times, LMK, bases,
    q_map, R, atom_types, scat_amps,
    detector_dist):

  L = np.reshape(LMK[:,0], (-1, 1, 1, 1))
  M = np.reshape(LMK[:,1], (-1, 1, 1, 1))
  K = np.reshape(LMK[:,2], (-1, 1, 1, 1))
  diffraction = np.zeros((detector_dist.shape[0], detector_dist.shape[1]),
      dtype=np.complex)

  for imol in range(len(atom_types)):
    for igeo in range(len(R[imol])):
      all_dists = calc_dists(R[imol][igeo])
      dist_inds1 = []
      dist_inds2 = []
      dist_scat_amps = []
      """
      for i, a1 in enumerate(atom_types):
        for j_, a2 in enumerate(atom_types[i+1:]):
          j = j_ + i+1
          dist_inds1.append(i)
          dist_inds2.append(j)
          dist_scat_amps.append(
              scat_amps[a1](q_map[:,:,0])*scat_amps[a2](q_map[:,:,0]))
      """
      for i, a1 in enumerate(atom_types[imol]):
        for j, a2 in enumerate(atom_types[imol]):
          if i == j:
            continue
          dist_inds1.append(i)
          dist_inds2.append(j)
          dist_scat_amps.append(
              scat_amps[a1](q_map[:,:,0])*scat_amps[a2](q_map[:,:,0]))

      dist_inds = (np.array(dist_inds1), np.array(dist_inds2))
      dist_scat_amps = np.array(dist_scat_amps)
      dists = all_dists[dist_inds]

      #print(lg.shape, q.shape, np.expand_dims(dists[:,0], axis=-1).shape)
      ADM_to_D = 8*np.pi**2/(2*LMK[:,0] + 1)
      """
      C = np.complex(0,1)**L*np.sqrt(4*np.pi*(2*L + 1))
      J = sp.special.spherical_jn(L,
          q_map[:,:,0]*np.expand_dims(np.expand_dims(dists[:,0], axis=-1), -1))
      Y = sp.special.sph_harm(-1*K, L,
          np.expand_dims(np.expand_dims(dists[:,2], axis=-1), axis=-1),
          np.expand_dims(np.expand_dims(dists[:,1], axis=-1), axis=-1))
      """
      
      C = (-1)**(2*L-K-M)*np.complex(0,1)**L*np.sqrt(4*np.pi*(2*L + 1))
      J = sp.special.spherical_jn(L,
          q_map[:,:,0]*np.expand_dims(np.expand_dims(dists[:,0], axis=-1), -1))
      Y = sp.special.sph_harm(-1*K, L,
          np.expand_dims(np.expand_dims(dists[:,2], axis=-1), axis=-1),
          np.expand_dims(np.expand_dims(dists[:,1], axis=-1), axis=-1))


      #diffraction = np.real(dist_scat_amps*C*J*Y)
      _diffraction = dist_scat_amps*C*J*Y
      diffraction += np.sum(np.sum(
          np.reshape(ADM_to_D*bases[:,itm], (-1, 1, 1, 1))*_diffraction\
          *sp.special.sph_harm(-1*M, L,
              np.expand_dims(q_map[:,:,2], axis=0),
              np.expand_dims(q_map[:,:,1], axis=0)),
              axis=1), axis=0)

  diffraction /= detector_dist
  #print("outp", calc_coeffs.shape)

  return np.real(diffraction), diffraction



def main(FLAGS, params):


  #############################################
  #####  Lab Frame Coordinate System      #####
  #####                                   #####
  #####  Z: along the laser polarization  ##### 
  #####  Y: electron beam direction       #####
  #####  X: set by Z, Y                   #####
  #############################################


  atom_types, atom_positions,\
  atom_distances_xyz, atom_distances_polar,\
  atom_distances_types = get_molecule_distances(params, logging)
  scat_amps = get_scattering_amplitudes(params, atom_types)

  db_wvl  = calculate_deBroglie_wavelength(params) #angs
  k0      = 2*np.pi/db_wvl

  """
  x_lab,z_lab = np.meshgrid(
      np.linspace(-1*params.detector_width/2, params.detector_width/2, params.imgBins),
      np.linspace(-1*params.detector_height/2, params.detector_height/2, params.imgBins))
  pixel_dist = np.sqrt(x_lab**2 + z_lab**2)
  detector_dist = np.sqrt(pixel_dist**2 + params.beamline_length**2)
  """
  q_x_lf,q_z_lf = np.meshgrid(
      -1*(np.arange(params.imgBins)-params.NradAzmBins+1)*params.QperPix,
      -1*(np.arange(params.imgBins)-params.NradAzmBins+1)*params.QperPix)
  q_x_lf = q_x_lf.astype(np.float32)
  q_z_lf = q_z_lf.astype(np.float32)
  q_map = np.sqrt(q_x_lf**2 + q_z_lf**2)

  # Detector xyz distance
  q_template = np.arange(params.NradAzmBins, dtype=np.float32)*params.QperPix
  deflection_angles = 2*np.arcsin(q_template*db_wvl/(4*np.pi))
  x_template = params.beamline_length*np.tan(deflection_angles)
  x_template = np.concatenate([np.flip(-1*x_template[1:]), x_template])
  x_lab, z_lab = np.meshgrid(x_template, x_template)
  pixel_dist = np.sqrt(x_lab**2 + z_lab**2).astype(np.float32)
  detector_dist = np.sqrt(pixel_dist**2 + params.beamline_length**2).astype(np.float32)

  """ 
  plt.pcolormesh(x_lab)
  plt.colorbar()
  plt.savefig("x_lab.png")
  plt.close()
  plt.pcolormesh(z_lab)
  plt.colorbar()
  plt.savefig("z_lab.png")
  plt.close()
  plt.pcolormesh(q_x_lf)
  plt.colorbar()
  plt.savefig("q_xlf.png")
  plt.close()
  plt.pcolormesh(q_z_lf)
  plt.colorbar()
  plt.savefig("q_zlf.png")
  plt.close()
  """

  if params.calculation_type == "debug":
    polar = np.linspace(0, np.pi, 5)[:-1]
    azim  = np.linspace(0, 2*np.pi, 5)[:-1]
    polar = np.array([np.pi/5])
    azim = np.array([np.pi/6])
    print(polar)
    print(azim)

    for p,a in zip(polar, azim):
      print(p,a)


    L = 2
    ylmtest = True
    for l in range(1,L):
      for m in range(-1*l, l+1):
        if True or np.abs(m) == l:
          #m *= 0
          if ylmtest:
            truth = np.conj(Ylm_calc(m, l, azim, polar))*np.sqrt(4*np.pi/(2*l+1))
            test = get_wignerD_3d(azim, polar, np.zeros_like(azim), l,
                np.array([m]), 0)
            print("\n######  L:{} M:{}  ######".format(l,m))
            #print(np.arange(-1*l,l+1))
            print(truth)
            print(test)
          else:
            truth = np.conj(Ylm_calc(m, l, azim, polar))
            test = get_wignerD_3d(azim, polar, np.zeros_like(azim), l,
                np.arange(-1*l,l+1), m)
            print("\n######  L:{} K:{}  ######".format(l,m))
            #print(np.arange(-1*l,l+1))
            print(truth)
            print(test)


  if params.calculation_type == "azmAvg":
    q = np.arange(params.NradAzmBins)*params.QperPix

    atm_diffractions = []
    mol_diffractions = []
    for imol in range(len(atom_distances_polar)):
      atm_diffractions.append(np.zeros_like(q))
      mol_diffractions.append(np.zeros_like(q))

      # Atomic diffraction
      for i,atm in enumerate(atom_types[imol]):
        print("ATM T: ",atm)
        atm_diffractions[imol] += np.abs(scat_amps[atm](q))**2
  
      # Molecular diffraction
      print("SIZE: ",atom_distances_polar[0].shape)
      if atom_distances_polar[imol] is not None:
        sincs = np.sum(
            np.sinc((q/np.pi)\
              *np.expand_dims(atom_distances_polar[imol][:,:,0], -1)),
            axis=0)
        for idst in range(len(atom_distances_types[imol])):
          mol_diffractions[imol] += sincs[idst]\
              *scat_amps[atom_distances_types[imol][idst][0]](q)\
              *scat_amps[atom_distances_types[imol][idst][1]](q)

        # Normalize by the number of geometries
        mol_diffractions[imol] /= atom_distances_polar[imol].shape[0]

    atm_diffraction = np.sum(np.array(atm_diffractions), axis=0)
    mol_diffraction = np.sum(np.array(mol_diffractions), axis=0)


    fName_template_dat =\
        "{0}_sim_{1}Diffraction-azmAvg_Qmax-{2:.4g}_Bins[{3}].dat"
    fName = os.path.join(params.simOutputDir, fName_template_dat.format(
        params.molecule, "atm", q[-1], params.NradAzmBins))
    with open(fName, "wb") as file:
      atm_diffraction.astype(np.double).tofile(file)

    fName = os.path.join(params.simOutputDir, fName_template_dat.format(
        params.molecule, "mol", q[-1], params.NradAzmBins))
    with open(fName, "wb") as file:
      mol_diffraction.astype(np.double).tofile(file)

    fName = os.path.join(params.simOutputDir, fName_template_dat.format(
        params.molecule, "sms", q[-1], params.NradAzmBins))
    with open(fName, "wb") as file:
      (q*mol_diffraction/atm_diffraction).astype(np.double).tofile(file)


    if len(mol_diffractions) > 1:
      for imol in range(len(mol_diffractions)):
        fName = os.path.join(params.simOutputDir, fName_template_dat.format(
            params.molecule+"_mol{}".format(imol), 
            "atm", q[-1], params.NradAzmBins))
        with open(fName, "wb") as file:
          atm_diffractions[imol].astype(np.double).tofile(file)

        fName = os.path.join(params.simOutputDir, fName_template_dat.format(
            params.molecule+"_mol{}".format(imol),
            "mol", q[-1], params.NradAzmBins))
        with open(fName, "wb") as file:
          mol_diffractions[imol].astype(np.double).tofile(file)

        fName = os.path.join(params.simOutputDir, fName_template_dat.format(
            params.molecule+"_mol{}".format(imol),
            "sms", q[-1], params.NradAzmBins))
        with open(fName, "wb") as file:
          (q*mol_diffraction[imol]/atm_diffraction[imol]).astype(np.double).tofile(file)



    fName_template_h5 =\
        "{0}_sim_diffraction-azmAvg_Qmax-{1:.4g}.h5"
    fName = os.path.join(params.simOutputDir, fName_template_h5.format(
        params.molecule, q_map[params.imgBins//2,-1]))
    with h5py.File(fName, 'w') as hf:
      hf.create_dataset("q",  data=q)
      hf.create_dataset("atm_diffraction",  data=atm_diffraction)
      hf.create_dataset("mol_diffraction",  data=mol_diffraction)
      hf.create_dataset("sms_diffraction",  data=q*mol_diffraction/atm_diffraction)


    """
    # Plot Results
    fig, ax = plt.subplots(1, 2, figsize=(12,5))
    ax[0].plot(q, atm_diffraction)
    ax[0].plot(q, mol_diffraction)
    ax[1].plot(q, q*mol_diffraction/atm_diffraction)
    ax[0].set_xlim([0, q[-1]])
    ax[1].set_xlim([0, q[-1]])
    fig.savefig(os.path.join("./plots",
        "{}_sim-diffraction-azmAvg.png".format(params.molecule)))
    """

  else:
  
    print("CLUSTER", FLAGS.cluster)
    LMK, bases, norms, times  = get_bases(params, logging,
        plot=(not FLAGS.cluster))
    #params, bases, times      = apply_fit_parameters(params, bases, times, logging)
    logging.info(LMK)

    # Set up pixel and s maps
    x,z = x_lab,z_lab
    zero_mask = pixel_dist > 0

    signal_xyz_qf = np.concatenate([
        np.expand_dims(x, axis=-1),
        np.expand_dims(np.ones_like(x)*params.beamline_length, axis=-1),
        np.expand_dims(z, axis=-1)],
        axis=-1)


    q_xyz_lf = np.concatenate([
      np.expand_dims(q_x_lf, -1),
      np.expand_dims(q_map*np.cos(
        (2*np.arcsin(q_map*db_wvl/(4*np.pi)) - np.pi)/2), -1),
      np.expand_dims(q_z_lf, -1)], -1)
    q_theta_lf  = np.zeros_like(z, dtype=np.float32)
    q_theta_lf[zero_mask]  = np.arccos(q_xyz_lf[zero_mask,2]/q_map[zero_mask])
    q_phi_lf  = np.zeros_like(z, dtype=np.float32)
    q_phi_lf[zero_mask]  = np.arctan2(q_xyz_lf[zero_mask,1], q_xyz_lf[zero_mask,0])

    """
    plt.pcolormesh(q_theta_lf)
    plt.colorbar()
    plt.savefig("q_theta.png")
    plt.close()
    plt.pcolormesh(q_phi_lf)
    plt.colorbar()
    plt.savefig("q_phi.png")
    plt.close()
    """

    q_map_polar = np.concatenate([
        np.expand_dims(q_map, axis=-1),
        np.expand_dims(q_theta_lf, axis=-1),
        np.expand_dims(q_phi_lf, axis=-1)],
        axis=-1)


    """ 
    incoming_e      = np.zeros((1,1,3))
    incoming_e[:,:,1] = k0
    s_xyz_lf    = k0*signal_xyz_qf/np.expand_dims(np.linalg.norm(signal_xyz_qf, axis=-1), axis=-1) - incoming_e
    q_map       = np.linalg.norm(s_xyz_lf, axis=-1)

    q_map = np.sqrt(x**2 + z**2)
    s_theta_lf  = np.zeros_like(z)
    s_theta_lf[zero_mask]  = np.arccos(s_xyz_lf[zero_mask,2]/q_map[zero_mask])
    s_phi_lf  = np.zeros_like(z)
    s_phi_lf[zero_mask]  = np.arctan2(s_xyz_lf[zero_mask,1], s_xyz_lf[zero_mask,0])

    q_map_polar = np.concatenate([
        np.expand_dims(q_map, axis=-1),
        np.expand_dims(s_theta_lf, axis=-1),
        np.expand_dims(s_phi_lf, axis=-1)],
        axis=-1)
    """

    atm_diffraction = make_atomic_diffraction(
        (params.imgBins, params.imgBins), detector_dist,
        q_map_polar, atom_types, scat_amps)

    # Saving Results
    """
    fName = os.path.join("output",
        "atm_diffraction_Q-{0:.5g}.npy".format(q_map[q_map.shape[0]//2,-1]))
    with open(fName, "wb") as file:
      np.save(file, atm_diffraction)
    
    print(atm_diffraction.shape)
    fName = os.path.join("output",
        "atm_diffraction_Q-{0:.5g}_shape[{1},{2}].dat".format(
          q_map[q_map.shape[0]//2,-1],
          atm_diffraction.shape[0], atm_diffraction.shape[1]))
    with open(fName, "wb") as file:
      atm_diffraction.astype(np.double).tofile(file)
    """


    if params.calculation_type == "numeric":
      ###  Ensemble smearing  ###
      smearing_polar  = np.linspace(0, np.pi, params.smear_polar_bins+1,
          dtype=np.float32)
      smearing_azim   = np.linspace(0, 2*np.pi, params.smear_azim_bins+1,
          dtype=np.float32)
      smearing_spin   = np.linspace(0, 2*np.pi, params.smear_spin_bins+1,
          dtype=np.float32)
      
      """
      smearing_polar = np.arange(params.smear_polar_bins+1)\
          *np.pi/(params.smear_polar_bins)
      smearing_azim = np.arange(params.smear_azim_bins+1)\
          *2*np.pi/(params.smear_azim_bins)
      smearing_spin = np.arange(params.smear_spin_bins+1)\
          *2*np.pi/(params.smear_spin_bins)
      #smearing_spin = np.concatenate([smearing_spin, smearing_spin[1:]+np.pi])
      """

      # Use midpoint for numerical integration
      smearing_polar  = (smearing_polar[1:] + smearing_polar[:-1])/2.
      smearing_azim   = (smearing_azim[1:] + smearing_azim[:-1])/2.
      smearing_spin   = (smearing_spin[1:] + smearing_spin[:-1])/2.

      """
      smearing_polar[params.smear_polar_bins//2:] =\
          smearing_polar[:params.smear_polar_bins//2]
      smearing_azim[params.smear_azim_bins//2:] =\
          smearing_azim[:params.smear_azim_bins//2]
      smearing_spin[params.smear_spin_bins//2:] =\
          smearing_spin[:params.smear_spin_bins//2]
      """
      print("VAR DELT", np.var(smearing_spin[1:] - smearing_spin[:-1]))

      #smearing_spin = np.array([4*np.pi/3])
      #print(smearing_spin)
      """
      smearing_polar  = np.arange(params.smear_polar_bins)\
          *np.pi/params.smear_polar_bins
      smearing_azim   = np.arange(params.smear_azim_bins)\
          *2*np.pi/params.smear_azim_bins
      smearing_spin   = np.arange(params.smear_azim_bins)\
          *2*np.pi/params.smear_azim_bins
      """

      smearing_ea     = np.zeros(
          (smearing_azim.shape[0], smearing_polar.shape[0],
            smearing_spin.shape[0], 3),
          dtype=np.float32)

     
      for i in range(len(smearing_azim)):
        smearing_ea[i,:,:,0]  = smearing_azim[i]
        for k in range(len(smearing_polar)):
          smearing_ea[:,k,:,1] = smearing_polar[k]
          smearing_ea[i,k,:,2]  = smearing_spin
      print("OSHAPE", smearing_ea.shape)
      smearing_ea       = np.reshape(smearing_ea, (-1,3))
      #smearing_ea = np.array([np.zeros(params.smear_azim_bins), smearing_azim, np.zeros(params.smear_azim_bins)]).transpose()
      print("SSSSSS", smearing_ea.shape)
      smearing_jacobian = np.sin(smearing_ea[:,1])
      smearing_jacobian *= np.pi/params.smear_polar_bins
      smearing_jacobian *= 2*np.pi/params.smear_azim_bins
      smearing_jacobian *= 2*np.pi/params.smear_spin_bins
      

      for tm in range(bases.shape[-1]):
        if FLAGS.time_ind is not None:
          if tm != FLAGS.time_ind:
            continue

        mol_diffraction, mol_diffraction_raw = diffraction_calculation_numerical(
            params.molecule, tm, times, LMK, bases.astype(np.float32),
            (params.imgBins, params.imgBins),
            detector_dist, q_map_polar,
            atom_distances_types, atom_distances_polar,
            smearing_ea, smearing_jacobian, scat_amps, logging)
        
        """
        mid = mol_diffraction.shape[0]//2
        mol_diffraction[mid,mid] = 0
        """
        rng = np.max(np.abs([np.amax(mol_diffraction), np.amin(mol_diffraction)]))
        #plt.pcolormesh(mol_diffraction, vmax=rng, vmin=-1*rng, cmap='seismic')
        #plt.colorbar()
        #plt.savefig("result.png")
        fName_template_h5 =\
            "{0}_sim_diffraction-numeric_Qmax-{1:.4g}{2}"
        fName = os.path.join(params.simOutputDir, fName_template_h5.format(
            params.molecule, q_map[params.imgBins//2,-1],
            "_time-{0:.6g}".format(times[tm])))
        if params.input_LMK is not None:
          suffix = "_LMK"
          for lmk in params.input_LMK:
            suffix += "-"+str(lmk[0])+"."+str(lmk[1])+"."+str(lmk[2])
          fName += suffix
        fName += ".h5"
        with h5py.File(fName, 'w') as hf:
          hf.create_dataset("q",  data=q_map)
          hf.create_dataset("atm_diffraction",  data=atm_diffraction)
          hf.create_dataset("mol_diffraction",  data=mol_diffraction)
          hf.create_dataset("mol_diffraction_raw",  data=mol_diffraction_raw)
          hf.create_dataset("mod_diffraction",  data=mol_diffraction/atm_diffraction)
          hf.create_dataset("detector_distance", data=detector_dist)


    elif params.calculation_type == "analytic":
      if params.QperPix is None:
        logging.fatal("Must specify QperPix for the analytic option.")
        sys.exit(0)
    
      #if params.molecule == "N2O":
      #  bases *= np.expand_dims(np.sqrt(4*np.pi/(2*LMK[:,0] + 1)), -1)**2
      for tm in range(bases.shape[-1]):
        if FLAGS.time_ind is not None:
          if tm != FLAGS.time_ind:
            continue

        mol_diffraction, mol_diffraction_raw = diffraction_calculation_analytical(
            tm, times, LMK, bases.astype(np.float32),
            q_map_polar, atom_positions, atom_types, scat_amps,
            detector_dist)

        rng = np.max(np.abs([np.amax(mol_diffraction), np.amin(mol_diffraction)]))
        #plt.pcolormesh(mol_diffraction, vmax=rng, vmin=-1*rng, cmap='seismic')
        #plt.colorbar()
        #plt.savefig("aresult.png")
        fName_template_h5 =\
            "{0}_sim_diffraction-analytic_Qmax-{1:.4g}{2}"
        fName = os.path.join(params.simOutputDir, fName_template_h5.format(
            params.molecule, q_map[params.imgBins//2,-1],
            "_time-{0:.6g}".format(times[tm])))
        if params.input_LMK is not None:
          suffix = "_LMK"
          for lmk in params.input_LMK:
            suffix += "-"+str(lmk[0])+"."+str(lmk[1])+"."+str(lmk[2])
          fName += suffix
        fName += ".h5"
        print(fName)
        with h5py.File(fName, 'w') as hf:
          hf.create_dataset("q",  data=q_map)
          hf.create_dataset("atm_diffraction",  data=atm_diffraction)
          hf.create_dataset("mol_diffraction",  data=mol_diffraction)
          hf.create_dataset("mol_diffraction_raw",  data=mol_diffraction_raw)
          hf.create_dataset("mod_diffraction",  data=mol_diffraction/atm_diffraction)
          hf.create_dataset("detector_distance", data=detector_dist)



if __name__ == '__main__':
  FLAGS = parser.parse_args()
  if FLAGS.molecule is not None:
    params = get_parameters(FLAGS.run, FLAGS.molecule)
  else:
    params = get_parameters(FLAGS.run)
  if FLAGS.molecule is not None:
    params.molecule = FLAGS.molecule
  FLAGS = setup(parser, output_dir=params.simOutputDir)

  # Flag handling
  if FLAGS.calculation_type is not None:
    params.calculation_type = FLAGS.calculation_type
  if FLAGS.xyz_file is not None:
    params.xyz_file = FLAGS.xyz_file
  if FLAGS.output_fileName_suffix is not None:
    params.output_fileName_suffix = FLAGS.output_fileName_suffix
  if FLAGS.basis_folder is not None:
    params.basis_folder = FLAGS.basis_folder
  if FLAGS.LMK is not None:
    params.input_LMK = FLAGS.LMK
  else:
    params.input_LMK = None
  if FLAGS.time_ind is not None:
    if FLAGS.time_ind == "nan":
      FLAGS.time_ind = None
    else:
      FLAGS.time_ind = int(FLAGS.time_ind)
  if FLAGS.eval_time is not None:
    params.sim_eval_times = np.array([float(FLAGS.eval_time)])
    FLAGS.time_ind = 0

  params.imgBins = 2*params.NradAzmBins - 1 
  
  if params.calculation_type != "azmAvg":
    from ADM import *

  main(FLAGS, params)
  

