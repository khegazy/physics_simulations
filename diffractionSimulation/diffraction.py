import sys, os, glob
import logging
import argparse
import pickle as pl
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
#import spherical_functions as sf
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
from scipy.special import factorial, jacobi
from multiprocessing import Pool
from functools import partial


#from absl import app
#from absl import flags
#from absl import logging

#import jax as jax
#import jax.numpy as jnp 

parser = argparse.ArgumentParser()

parser.add_argument(
    '--molecule', type=str, required=True,
    help='Name of the molecule.'
)
parser.add_argument(
    '--xyz_file', type=str, required=True,
    help='Address to XYZ file.'
)
parser.add_argument(
    '--detector_width', default=0.04, type=float,
    help='Detector width in meters.'
)
parser.add_argument(
    '--width_bins', default=29, type=int,
    help='Number of width bins.'
)
parser.add_argument(
    '--detector_height', default=0.04, type=float,
    help='Detector height in meters.'
)
parser.add_argument(
    '--height_bins', default=29, type=int,
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
    '--smear_time', default=None, type=float,
    help='Electron pulse width (STD) in fs.'
)
parser.add_argument(
    '--scat_amps_path', default="../scatteringAmplitudes/3.7MeV/", type=str,
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
  "I" : "iodine"
}

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


def get_scattering_amplitudes(FLAGS, atom_types):

  scat_amps = {}
  dBrog = calculate_deBroglie_wavelength(FLAGS)
  for atm in atom_types:
    if atm in scat_amps:
      continue

    angStr = []
    sctStr = []
    fName = os.path.join(FLAGS.scat_amps_path, atom_names[atm] + "_dcs.dat")
    with open(fName, 'r') as inpFile:
      ind=0
      for line in inpFile:
        if ind < 31:
          ind += 1
          continue

        angStr.append(line[2:11])
        sctStr.append(line[39:50])   
    
    angs = np.array(angStr).astype(np.float64)
    q = 4*np.pi*np.sin(angs/2.)/dBrog
    scts = np.sqrt(np.array(sctStr).astype(np.float64))

    scat_amps[atm] = interp1d(q, scts, 'cubic')

  return scat_amps


def get_bases(FLAGS, logging, normalize=False, plot=False):
  bases = {}
  norms = {}
  LMK = []
  if FLAGS.get_bases_type is 0:
    start_time, end_time = None, None
    files = "{}/*.dat".format(FLAGS.basis_folder)
    print("files", files)
    for fl in glob.glob(files):
      logging.info("Importing Bases: " + fl)
      L = int(fl[fl.find("_L-")+3:fl.find("_M-")])
      with open(fl, "rb") as file:
        bases[L] = np.fromfile(file, np.double)

        # Normalize bases
        bases_womean = {}
        if L != 0:
          bases_womean[L] = bases[L] - np.mean(bases[L])
        else:
          bases_womean[L] = bases[L]
        #print("norm",L,np.sqrt(np.sum(bases[L]**2)), np.amax(np.abs(bases[L])))
        norms[L] = np.sqrt(np.sum(bases_womean[L]**2))
        if normalize:
          bases[L] = bases_womean[L]
          bases[L] /= norms[L]
        
        # Get Time
        aind = fl.find("_time-")
        bind = fl.find("-", aind+7)
        cind = fl.find("_bins")
        sTime = float(fl[aind+6:bind])
        eTime = float(fl[bind+1:cind])
        if start_time is None:
          start_time  = sTime
          end_time    = eTime
        else:
          if (start_time != sTime) or (end_time != eTime):
            logging.critical("Bases do not have the same time range.")
            sys.exit(1)

    basisList = []
    normsList = []
    keys = np.sort(list(bases.keys()))
    for l in keys:
      basisList.append(np.expand_dims(bases[l], axis=0))
      normsList.append(norms[l])
      LMK.append(np.array([l,0,0]))
    tLen = basisList[0].shape[-1]
    times = start_time + np.linspace(0, 1, tLen)*(end_time - start_time)

  elif FLAGS.get_bases_type is 1:
    L_list = []
    files = "{}/*.npy".format(FLAGS.basis_folder)
    print("search", files)
    for fl in glob.glob(files):
      logging.info("Importing Bases: " + fl)
      D = fl[fl.find(" D")+2:-4]
      ln = len(D) + len(D)%2
      L = int(D[:ln//2])
      M = int(D[ln//2:])
      print(D,L,M)
      if L not in L_list:
        L_list.append(L)
      with open(fl, "rb") as file:
        bases[(L,M)] = np.load(file)
        print(bases[(L,M)])

        # Normalize bases
        bases_womean = {}
        if L != 0:
          bases_womean[(L,M)] = bases[(L,M)] - np.mean(bases[(L,M)])
        else:
          bases_womean[(L,M)] = bases[(L,M)]
        #print("norm",L,np.sqrt(np.sum(bases[L]**2)), np.amax(np.abs(bases[L])))
        norms[(L,M)] = np.sqrt(np.sum(bases_womean[(L,M)]**2))
        if normalize:
          bases[(L,M)] = bases_womean[(L,M)]
          bases[(L,M)] /= norms[(L,M)]
 
    sys.exit(0)
    basisList = []
    normsList = []
    print("LIST", L_list)
    L_list.sort()
    for l in L_list:
      for m in np.linspace(0, l, 2, dtype=int):
        basisList.append(np.expand_dims(bases[(l,m)], axis=0))
        normsList.append(norms[(l,m)])
        LMK.append(np.array([l,0,m]))

    fName = FLAGS.basis_folder.replace("/A/", "/times/").replace("temp-", "temp_")\
        + ".npy"
    logging.info("Time file: " + fName)
    with open(fName, "rb") as file:
      times = np.load(file)

  bases = np.concatenate(basisList, axis=0)
  if plot:
    inds = np.arange(len(times))# < 100
    for i, (L, M, K) in enumerate(LMK):
      plt.plot(times[inds], bases[i,inds])
      plt.savefig(os.path.join("./plots", "basis_input_{}-{}-{}.png".format(
          L, M, K)))
      plt.close()


  if FLAGS.smear_time is not None:
    delta_time = times[1] - times[0]
    bases = gaussian_filter1d(
        bases, FLAGS.smear_time/delta_time, axis=-1)


  if plot:
    inds = np.arange(len(times))# < 100
    for i, (L, M, K) in enumerate(LMK):
      plt.plot(times[inds], bases[i,inds])
      plt.savefig(os.path.join("./plots", "basis_smeared_{}-{}-{}.png".format(
          L, M, K)))
      plt.close()


  return np.array(LMK), bases, np.array(normsList), times


def apply_fit_parameters(FLAGS, bases, times, logging):
  if FLAGS.basis_params is None and FLAGS.basis_eval_params is None:
    return None, bases, None
  
  if FLAGS.basis_params is not None and FLAGS.basis_eval_params is not None:
    logging.fatal("Can only specify basis_params OR basis_eval_params")
    sys.exit(1)

  if FLAGS.basis_params is not None:
    with open(FLAGS.basis_params, 'rb') as file:
      params = pl.load(file)
  
    eval_times = params['startTime'] + np.arange(40)*0.1
    bases_interp = interp1d(times, bases, axis=-1)
    bases = bases_interp(eval_times)

    return params, bases, eval_times

  if FLAGS.basis_eval_params is not None:
    eval_times = np.linspace(
        FLAGS.basis_eval_params[0],
        FLAGS.basis_eval_params[1],
        FLAGS.basis_eval_params[2])
    
    bases_interp = interp1d(times, bases, axis=-1)
    bases = bases_interp(eval_times)
    
    return None, bases, eval_times

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
  if not os.path.exists(FLAGS.xyz_file):
    logging.fatal("Cannot find xyz file: " + FLAGS.xyz_file)
    sys.exit(1)
  
  atom_types      = []
  atom_positions  = []
  with open(FLAGS.xyz_file) as file:
    for i,ln in enumerate(file):
      if i == 0:
        Natoms = int(ln)
      elif i > 1:
        vals = ln.split()
        print("l",ln, vals)
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


def get_wignerD(phi_ea, theta_ea, chi_ea, l, m):
  lmm = np.ones((2*l+1, 3))*l
  lmm[:,1] -= np.arange(2*l+1)
  lmm[:,2] = m
  lmm = lmm.astype(int)

  D_ME = sf.Wigner_D_element(phi_ea, theta_ea, chi_ea, lmm)
  return D_ME, lmm


def rotate_Ylm(phi_ea, theta_ea, chi_ea, l, m, polar, azim):
  D_ME, lmm = get_wignerD(phi_ea, theta_ea, chi_ea, l, m)
  D_ME = np.expand_dims(np.expand_dims(D_ME, axis=-1), axis=-1)
  #Ylms = sp.special.sph_harm(lmm[:,1],1, s[:,1], s[:,2])
  Ylms = sp.special.sph_harm(
    np.expand_dims(np.expand_dims(lmm[:,1],-1), -1),
    np.expand_dims(np.expand_dims(lmm[:,0],-1), -1),
    np.expand_dims(azim,0),        
    np.expand_dims(polar,0))
  #print(lmm)
  #print(D_ME)
  #print(Ylms)
  #print("shapes ", D_ME.shape, Ylms.shape)
  return np.sum(D_ME*Ylms, axis=0)


def get_wignerD_map(phi_ea, theta_ea, chi_ea, l, m):
  lmm = np.ones((2*l+1, 3))*l
  lmm[:,1] -= np.arange(2*l+1)
  lmm[:,2] = m
  lmm = lmm.astype(int)

  print(phi_ea.shape, theta_ea.shape, chi_ea.shape, lmm.shape)
  D_ME = sf.Wigner_D_element(phi_ea, theta_ea, chi_ea, np.expand_dims(lmm, axis=0))
  return D_ME, lmm

def rotate_Ylm_map(phi_ea, theta_ea, chi_ea, l, m, polar, azim):
  D_ME, lmm = get_wignerD_map(phi_ea, theta_ea, chi_ea, l, m)
  D_ME = np.expand_dims(np.expand_dims(D_ME, axis=-1), axis=-1)
  #Ylms = sp.special.sph_harm(lmm[:,1],1, s[:,1], s[:,2])
  Ylms = sp.special.sph_harm(
    np.expand_dims(np.expand_dims(lmm[:,1],-1), -1),
    np.expand_dims(np.expand_dims(lmm[:,0],-1), -1),
    np.expand_dims(azim,0),        
    np.expand_dims(polar,0))
  #print(lmm)
  #print(D_ME)
  #print(Ylms)
  #print("shapes ", D_ME.shape, Ylms.shape)
  return np.sum(D_ME*Ylms, axis=0)

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
  print(a.shape, b.shape, k.shape, theta_ea.shape)
  d_coeff = (-1)**lmbd*np.sqrt(nCk_np(2*l-k, k+a))/np.sqrt(nCk_np(k+b, b))\
      *np.sin(theta_ea/2)**a*np.cos(theta_ea/2)**b
  d = np.zeros_like(d_coeff)
  cos_theta = np.cos(theta_ea)
  for im in range(m_sum.shape[-1]):
    p = jacobi(k[0,im], a[0,im], b[0,im])
    print("SHAPE", np.polyval(p, cos_theta).shape, d_coeff[:,im].shape)
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


def make_diffraction(
    FLAGS, s_map_polar,
    atom_types, atom_distances_polar,
    smearing_ea, smearing_weights,
    scat_amps):

  smearing_weights = np.expand_dims(np.expand_dims(smearing_weights, -1), -1)
  diffraction = np.zeros((FLAGS.height_bins, FLAGS.width_bins))
  for ir, r_diff in enumerate(atom_distances_polar):

    # Rotate Y10 so it points in direction of r in lf
    #D_ME_, lmm = get_wignerD(r_diff[2], r_diff[1], 0, 1, 0)
    D_ME, l, m_sum, m = get_wignerD_3d(
        np.array([r_diff[2]]), np.array([r_diff[1]]), np.zeros(1), 1, 0)#np.zeros_like(r_diff[2]), 2, 1)
    print(D_ME.shape, D_ME)
    D_ME = D_ME[:,0]

    print("D", D_ME)
    print(m_sum, D_ME.shape)
    #for il,(l,m,_) in enumerate(lmm):
    for im, m in enumerate(m_sum):
      #print("loop1", l,m,D_ME[il])
      if np.abs(D_ME[im]) <= 1e-10:
        #print("skipping")
        continue

      s_dot_r_lf_ = D_ME[im]*(s_map_polar[:,:,0]*r_diff[0]/(0.5*np.sqrt(3/np.pi)))\
          *rotate_Ylm_3d(smearing_ea[:,0], smearing_ea[:,1], smearing_ea[:,2],\
            l, m, s_map_polar[:,:,1], s_map_polar[:,:,2])
      s_dot_r_lf  = np.real(s_dot_r_lf_)
      print("SHAPE!!!!",(smearing_weights*np.cos(s_dot_r_lf)).shape)
      smearing = np.sum(smearing_weights*np.cos(s_dot_r_lf), axis=0)
      diffraction += smearing*scat_amps[atom_types[ir][0]](s_map_polar[:,:,0])\
          *scat_amps[atom_types[ir][1]](s_map_polar[:,:,0])
      #for i,(weight_ea, euler_angles) in enumerate(zip(smearing_weights, smearing_ea)):
      #  s_dot_r_lf_ = D_ME[il]*(s_map_polar[:,:,0]*r_diff[0]/(0.5*np.sqrt(3/np.pi)))\
      #      *rotate_Ylm(euler_angles[0], euler_angles[1], euler_angles[2],\
      #        l, m, s_map_polar[:,:,1], s_map_polar[:,:,2])
      #  s_dot_r_lf  = np.real(s_dot_r_lf_)
      #  diffraction += weight_ea*np.cos(s_dot_r_lf)\
      #      *scat_amps[atom_types[ir][0]](s_map_polar[:,:,0])*scat_amps[atom_types[ir][1]](s_map_polar[:,:,0])
  return diffraction



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

  LMK, bases, norms, times  = get_bases(FLAGS, logging, plot=True)
  params, bases, times      = apply_fit_parameters(FLAGS, bases, times, logging)

  db_wvl  = calculate_deBroglie_wavelength(FLAGS) #angs
  k0      = 2*np.pi/db_wvl

  # Set up pixel and s maps
  x,z = np.meshgrid(
      np.linspace(-1*FLAGS.detector_width/2, FLAGS.detector_width/2, FLAGS.width_bins),
      np.linspace(-1*FLAGS.detector_height/2, FLAGS.detector_height/2, FLAGS.height_bins))
  pixel_dist = np.sqrt(x**2 + z**2)
  zero_mask = pixel_dist > 0

  signal_xyz_qf = np.concatenate([
      np.expand_dims(x, axis=-1),
      np.expand_dims(np.ones_like(x)*FLAGS.beamline_length, axis=-1),
      np.expand_dims(z, axis=-1)],
      axis=-1)

  incoming_e      = np.zeros((1,1,3))
  incoming_e[:,:,1] = k0
  s_xyz_lf    = k0*signal_xyz_qf/np.expand_dims(np.linalg.norm(signal_xyz_qf, axis=-1), axis=-1) - incoming_e
  s_map    = np.linalg.norm(s_xyz_lf, axis=-1)
  #s_map       = 4*np.pi*np.sin(np.arctan(pixel_dist/FLAGS.beamline_length)/2)/db_wvl
  s_theta_lf  = np.zeros_like(z)
  #s_theta_lf[zero_mask]  = np.arctan2(x[zero_mask], z[zero_mask])
  s_theta_lf[zero_mask]  = np.arccos(s_xyz_lf[zero_mask,2]/s_map[zero_mask])
  s_phi_lf  = np.zeros_like(z)
  s_phi_lf[zero_mask]  = np.arctan2(s_xyz_lf[zero_mask,1], s_xyz_lf[zero_mask,0])

  s_map_polar = np.concatenate([
      np.expand_dims(s_map, axis=-1),
      np.expand_dims(s_theta_lf, axis=-1),
      np.expand_dims(s_phi_lf, axis=-1)],
      axis=-1)

  im = plt.imshow(s_map)
  plt.colorbar(im)
  plt.savefig("testS.png")
  plt.close()
  im = plt.imshow(s_theta_lf)
  plt.colorbar(im)
  plt.savefig("testTH.png")
  plt.close()
  im = plt.imshow(s_phi_lf)
  plt.colorbar(im)
  plt.savefig("testPH.png")
  plt.close()


  # Ensemble smearing

  #smearing_weights = np.array([1])
  #smearing_ea = np.array([[0,np.pi/2,0]])

  smearing_polar  = np.arange(FLAGS.smear_polar_bins)*np.pi/FLAGS.smear_polar_bins
  smearing_azim   = np.arange(FLAGS.smear_azim_bins)*2*np.pi/FLAGS.smear_azim_bins
  smearing_spin   = np.arange(FLAGS.smear_azim_bins)*2*np.pi/FLAGS.smear_azim_bins
  azim_inds       = np.arange(len(smearing_azim)).astype(int)
  spin_inds       = np.arange(len(smearing_spin)).astype(int)
  smearing_ea     = np.zeros(
      (smearing_polar.shape[0], smearing_azim.shape[0], smearing_spin.shape[0], 3))

 
  """
  for i in range(FLAGS.smear_polar_bins):
    smearing_ea[i,:,:,1] = smearing_polar[i]

  smearing_ea[:,azim_inds,:,0]  = smearing_azim[azim_inds]
  """
  for i in range(FLAGS.smear_polar_bins):
    smearing_ea[i,:,:,0]  = smearing_azim[i]
    for k in range(FLAGS.smear_azim_bins):
      smearing_ea[:,i,k,1] = smearing_polar[i]
      smearing_ea[i,azim_inds,:,2]  = smearing_spin[azim_inds]
  smearing_ea       = np.reshape(smearing_ea, (-1,3))
  smearing_jacobian = np.sin(smearing_ea[:,1])

  diffraction = make_diffraction(
      FLAGS, s_map_polar,
      atom_distances_types, atom_distances_polar,
      smearing_ea, smearing_jacobian, scat_amps)

  print(atom_distances_polar.shape)
  rand_diff = np.zeros_like(diffraction)
  for ir, r_diff in enumerate(atom_distances_polar):
      #scat_amps[atom_types[ir][0]](s_map_polar[:,:,0])\
      print(r_diff)
      rand_diff += scat_amps[atom_distances_types[ir][0]](s_map)*scat_amps[atom_distances_types[ir][1]](s_map)*np.sinc(s_map*r_diff[0]/np.pi)
      #rand_diff += scat_amps[atom_types[ir][0]](s_map)*scat_amps[atom_types[ir][1]](s_map)#*np.sinc(s_map*r_diff[0])

  iind = diffraction.shape[0]//2
  plt.plot(s_map[iind,:], diffraction[iind,:]/np.amax(diffraction[iind,:]), 'k-')
  plt.plot(s_map[iind,:], rand_diff[iind,:]/np.amax(rand_diff[iind,:]), 'b-')
  plt.savefig("testDiffLO.png")
  plt.close()

  print("PASSED FIRST")
  """
  
  smearing_polar_at   = np.arange(FLAGS.smear_polar_bins)*np.pi/FLAGS.smear_polar_bins
  smearing_azim_at    = np.arange(FLAGS.smear_azim_bins)*2*np.pi/FLAGS.smear_azim_bins
  azim_inds           = np.arange(len(smearing_azim_at)).astype(int)
  smearing_ea         = np.zeros((smearing_polar_at.shape[0], smearing_azim_at.shape[0], 3))

  for i in range(FLAGS.smear_polar_bins):
    smearing_ea[i,:,1] = smearing_polar[i]

  smearing_ea[:,azim_inds,0]  = smearing_azim[azim_inds]
  smearing_ea       = np.reshape(smearing_ea, (-1,3))
  smearing_weights  = np.sin(smearing_ea[:,1])
  """

  print("SHAPES", bases.shape)
  sh_l = np.arange(bases.shape[0])*2
  sh_m = np.zeros_like(sh_l)
  #bases[0,:] = 0 #FIX ME
  print("TIMES", times)
  for tm in range(bases.shape[-1]):
    """
    print("weights",bases[:,tm])
    smearing_weights = np.imag(np.sum(np.expand_dims(bases[:,tm], axis=-1)*sp.special.sph_harm(
          np.expand_dims(sh_m, axis=-1), np.expand_dims(sh_l, axis=-1),\
          np.expand_dims(smearing_ea[:,1], axis=0), np.expand_dims(smearing_ea[:,0], axis=0)),
        axis=0)).astype(np.float64)
    print("IMAGE", np.sum(np.abs(smearing_weights)))
    """
    smearing_weights = np.real(np.sum(np.expand_dims(bases[:,tm], axis=-1)*sp.special.sph_harm(
          np.expand_dims(sh_m, axis=-1), np.expand_dims(sh_l, axis=-1),\
          np.expand_dims(smearing_ea[:,1], axis=0), np.expand_dims(smearing_ea[:,0], axis=0)),
        axis=0)).astype(np.float64)
    smearing_weights *= smearing_jacobian

    diffraction = make_diffraction(
        FLAGS, s_map_polar,
        atom_distances_types, atom_distances_polar,
        smearing_ea, smearing_weights, scat_amps)

    #print(diffraction[int(FLAGS.width_bins/2),int(FLAGS.width_bins/2):])
    #print(np.sinc(s_map[int(FLAGS.width_bins/2),int(FLAGS.width_bins/2):]/np.pi)*diffraction[int(FLAGS.width_bins/2),int(FLAGS.width_bins/2)])
    #print(s_map[int(FLAGS.width_bins/2),int(FLAGS.width_bins/2):])
    #plt.imshow(np.real(s_dot_r_lf))
    #plt.savefig("testsdr.png")
    #plt.close()
    #s_dot_r_lf = s_map*s_dot_r_lf
    #s_dot_r_lf = np.real(s_dot_r_lf)
    #print(s_dot_r_lf)
    #diffraction = np.cos(s_dot_r_lf)

    # Saving Results
    fName = os.path.join("output",
        "mol_diffraction_time-{0:.3g}.npy".format(times[tm]))
    with open(fName, "wb") as file:
      np.save(file, diffraction)
    
    fName = os.path.join("output",
        "mol_diffraction_time-{0:.3g}_shape-{1}-{2}.dat".format(
          times[tm], diffraction.shape[0], diffraction.shape[1]))
    with open(fName, "wb") as file:
      diffraction.astype(np.double).tofile(file)


    plc = plt.pcolormesh(s_xyz_lf[:,:,0], s_xyz_lf[:,:,2], diffraction)
        #vmin=-250000, vmax=250000, cmap='seismic')
    plt.gca().set_aspect('equal')
    #im = plt.imshow(diffraction)
    plt.colorbar(plc)
    plt.savefig("./plots/molDiffraction_aniso_time-{}.png".format(tm))
    plt.close()
  sys.exit(0)






  """
  s_map       = 4*np.pi*np.sin(pixel_theta_lf/2)/db_wvl
  s_theta_lf  = np.arcsin(k0*pixel_theta_lf/(np.linalg.norm(pixel_theta_lf)*s_map))
  s_phi_lf    = pixel_phi_lf
  # Passive rotation to the mf
  lmm = np.ones((3,3))
  lmm[0,1] = -1
  lmm[1,1] = 0 
  a,b,c = -1*np.pi/2, -1*np.pi/2, 0
  euler_angles = (-1*a, -1*b, -1*c)
  Y11_mf = rotate_Ylm(*euler_angles, lmm, s_theta_lf, s_phi_lf)
  lmm[:,2] = 0
  Y10_mf = rotate_Ylm(*euler_angles, lmm, s_theta_lf, s_phi_lf)
  lmm[:,2] = -1
  Y1n1_mf = rotate_Ylm(*euler_angles, lmm, s_theta_lf, s_phi_lf)
  Y1m_norm = 0.5*np.sqrt(3/np.pi)
  s_theta_mf  = np.arccos(Y10_mf/Y1m_norm)
  s_phi_mf    = np.arcsin((-1*(Y11_mf + Y1n1_mf)/(np.sin(s_theta_mf)*Y1m_norm)))
  print(s_theta_lf)
  print(s_phi_lf)
  print("MF")
  print(s_theta_mf)
  print(s_phi_mf)

  s_polar_map = np.concatenate([
      np.reshape(s_map, -1),
      np.reshape(s_theta_map, -1),
      np.reshape(s_phi_map, -1)
  ])
  """

  """
  def sdr_lf(s, r, euler_angles):
    lmm = np.ones((3,3))
    lmm[0,1] = -1
    lmm[1,1] = 0 
    lmm[:,2] = 0
    return s[:,0]*r[0]*rotate_Ylm(*euler_angles, lmm, s[:,1], s[:2])

  print(sdr_lf(np.array([[1,0,0]]), [1], (0,0,0)))
  print(sdr_lf(np.array([[1,0,0], [1,np.pi, 0]]), [1], (0,0,0)))
  """
  """
  def sdr_lf(s_polar_lf, r_mag, mol_rotation_mf):
    sdr = r_mag*s_polar_lf[:,:,0]*(\
        2*np.sqrt(2)*np.sin(mol_rotation_mf[0])*np.sin(mol_rotation_mf[1])*np.sin(s_polar_lf[:,:,1])*np.cos(s_polar_lf[:,:,2])\
        + np.sqrt(2)*np.sin(mol_rotation_mf[0])*np.cos(mol_rotation_mf[1])*np.cos(s_polar_lf[:,:,1])\
        + np.sqrt(2)*np.cos(mol_rotation_mf[0])*np.sin(s_polar_lf[:,:,1])*np.sin(s_polar_lf[:,:,2]))
    
    return sdr


  signal = np.cos(sdr_lf(s_polar_lf, 1, np.array([0,0])))
  """

if __name__ == '__main__':
  FLAGS = setup(parser)
  main(FLAGS)
  
  
"""
def get_wignerD_3d_(phi_ea, theta_ea, chi_ea, l, m):

  m_sum[0,:,0] = np.ones(2*l+1)*l - np.arange(2*l+1)
  m_sum *= -1 # Must flip the sign for the spin-weighted SH
  n_sum = np.expand_dims(np.expand_dims(np.arange(l-m), 0), 0)
  theta_ea = np.expand_dims(np.expand_dims(theta_ea, -1), -1)
  phi_ea   = np.expand_dims(np.expand_dims(phi_ea, -1), -1)
  chi_ea   = np.expand_dims(chi_ea, -1)


  print(m_sum.shape, n_sum.shape, theta_ea.shape, phi_ea.shape, chi_ea.shape)
  print(theta_ea, phi_ea)
  # Factors of (-1)^m cancel
  D_coeff = np.sqrt(4*np.pi/(2*l+1))
  SWSH_coeff = (np.sin(theta_ea/2)**(2*l))*np.sqrt(
      factorial(l+m_sum)*factorial(l-m_sum)*factorial(2*l+1)\
      /(4*np.pi*factorial(1, l+m)*factorial(1, l-m)))\
      *np.exp(np.complex(0,1)*m_sum)
  #for n in range(0, l-m):
  SWSH_sum_array = (-1)**(l - n_sum - m)\
      *nCk_np(l-m, n_sum)*nCk_np(l+m, n_sum+m-m_sum)\
      /((np.tan(theta_ea/2) + 1e-10)**(2*n_sum+m-m_sum))
  print(SWSH_coeff.shape)
  print(SWSH_coeff)
  print(SWSH_sum_array.shape)
  print(SWSH_sum_array)
  print(1/((np.tan(theta_ea/2) + 1e-10)**(2*n_sum+m-m_sum)))
  print(nCk_np(l-m, n_sum),nCk_np(l+m, n_sum+m-m_sum))
  print((n_sum+m-m_sum).shape, n_sum+m-m_sum)
  SWSH = SWSH_coeff*SWSH_sum_array
  D_ME = np.sum(SWSH, -1)*np.exp(np.complex(0, -1*m*chi_ea))

  return D_ME, l, -1*m_sum, m
"""


