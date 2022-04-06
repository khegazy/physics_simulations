import sys, os, glob
import inspect
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

from modules.diffraction_simulation import diffraction_calculation_numeric, atomic_diffraction_calculation
from modules.diffraction_simulation import molecular_diffraction_calculation as diffraction_calculation_analytic

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
    '--output_folder', default=None, type=str,
    help='Folder for output file'
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
      fName = os.path.join(params["scat_amps_dir"], atom_names[atm] + "_dcs.dat")
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
    np.sqrt((params["elEnergy"]*eV_to_au + C_AU**2)**2\
    - (C_AU)**4) #au
  db_lambda /= angs_to_au
  
  return db_lambda
  #return 0.00296622 #angs


def get_molecule_distances(params, logging):
  
  atom_types      = []
  atom_positions  = []
  atom_distances_xyz    = []
  atom_distances_polar  = []
  atom_distances_types  = []

  if not isinstance(params["xyz_file"], list):
    params["xyz_file"] = [params["xyz_file"]]

  for xyz_name in params["xyz_file"]:
    fName = os.path.join(params["xyz_dir"], xyz_name)
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
      if len(cur_atom) != Natoms:
        raise RuntimeError("The number of atoms does not match the retrieved atoms!!!\n")
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




def main(params, on_cluster=False, time_ind=False):


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
      np.linspace(-1*params["detector_width"]/2, params["detector_width"]/2, params["imgBins"]),
      np.linspace(-1*params["detector_height"]/2, params["detector_height"]/2, params["imgBins"]))
  pixel_dist = np.sqrt(x_lab**2 + z_lab**2)
  detector_dist = np.sqrt(pixel_dist**2 + params["beamline_length"]**2)
  """
  q_x_lf,q_z_lf = np.meshgrid(
      -1*(np.arange(params["imgBins"])-params["NradAzmBins"]+1)*params["QperPix"],
      -1*(np.arange(params["imgBins"])-params["NradAzmBins"]+1)*params["QperPix"])
  q_x_lf = q_x_lf.astype(np.float32)
  q_z_lf = q_z_lf.astype(np.float32)
  q_map = np.sqrt(q_x_lf**2 + q_z_lf**2)

  # Detector xyz distance
  q_template = np.arange(params["NradAzmBins"], dtype=np.float32)*params["QperPix"]
  deflection_angles = 2*np.arcsin(q_template*db_wvl/(4*np.pi))
  x_template = params["beamline_length"]*np.tan(deflection_angles)
  x_template = np.concatenate([np.flip(-1*x_template[1:]), x_template])
  x_lab, z_lab = np.meshgrid(x_template, x_template)
  pixel_dist = np.sqrt(x_lab**2 + z_lab**2).astype(np.float32)
  detector_dist = np.sqrt(pixel_dist**2 + params["beamline_length"]**2).astype(np.float32)

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

  if params["calculation_type"] == "debug":
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


  if params["calculation_type"] == "azmAvg":
    q = np.arange(params["NradAzmBins"])*params["QperPix"]

    atm_diffractions = []
    mol_diffractions = []

    for imol in range(len(atom_distances_polar)):
      atm_diffractions.append(np.zeros_like(q))
      mol_diffractions.append(np.zeros_like(q))

      # Atomic diffraction
      for i,atm in enumerate(atom_types[imol]):
        atm_diffractions[imol] += np.abs(scat_amps[atm](q))**2
  
      # Molecular diffraction
      if atom_distances_polar[imol] is not None:
        if "smear_std_ratio" not in params:
          sincs = np.mean(
              np.sinc((q/np.pi)\
                *np.expand_dims(atom_distances_polar[imol][:,:,0], -1)),
              axis=0)
        else:
          if atom_distances_polar[imol].shape[0] != 1:
            raise RuntimeError("Cannot smear an ensemble of molecules from"\
                + "an xyz file, remove 'smear_std_ratio' from parameters")
          std_steps = np.linspace(-5, 5, 101)
          new_lens  = atom_distances_polar[imol][:,:,0]\
              + np.expand_dims(std_steps, -1)*params["smear_std_ratio"]\
              *np.sqrt(atom_distances_polar[imol][:,:,0])
          weights = np.exp(-0.5*std_steps**2)
          weights /= np.sum(weights)
          sincs = np.sum(np.reshape(weights, (-1, 1, 1))\
              *np.sinc((q/np.pi)*np.expand_dims(new_lens, -1)),
              axis=0)
        for itp in range(len(atom_distances_types[imol])):
          mol_diffractions[imol] += sincs[itp]\
              *scat_amps[atom_distances_types[imol][itp][0]](q)\
              *scat_amps[atom_distances_types[imol][itp][1]](q)

    atm_diffraction = np.sum(np.array(atm_diffractions), axis=0)
    mol_diffraction = np.sum(np.array(mol_diffractions), axis=0)


    fName_template_dat =\
        "{0}_sim_{1}Diffraction-azmAvg_Qmax-{2:.4g}_Bins[{3}].dat"
    fName = os.path.join(params["simOutputDir"], fName_template_dat.format(
        params["molecule"], "atm", q[-1], params["NradAzmBins"]))
    with open(fName, "wb") as file:
      atm_diffraction.astype(np.double).tofile(file)

    fName = os.path.join(params["simOutputDir"], fName_template_dat.format(
        params["molecule"], "mol", q[-1], params["NradAzmBins"]))
    with open(fName, "wb") as file:
      mol_diffraction.astype(np.double).tofile(file)

    fName = os.path.join(params["simOutputDir"], fName_template_dat.format(
        params["molecule"], "sms", q[-1], params["NradAzmBins"]))
    with open(fName, "wb") as file:
      (q*mol_diffraction/atm_diffraction).astype(np.double).tofile(file)


    if len(mol_diffractions) > 1:
      for imol in range(len(mol_diffractions)):
        fName = os.path.join(params["simOutputDir"], fName_template_dat.format(
            params["molecule"]+"_mol{}".format(imol), 
            "atm", q[-1], params["NradAzmBins"]))
        with open(fName, "wb") as file:
          atm_diffractions[imol].astype(np.double).tofile(file)

        fName = os.path.join(params["simOutputDir"], fName_template_dat.format(
            params["molecule"]+"_mol{}".format(imol),
            "mol", q[-1], params["NradAzmBins"]))
        with open(fName, "wb") as file:
          mol_diffractions[imol].astype(np.double).tofile(file)

        fName = os.path.join(params["simOutputDir"], fName_template_dat.format(
            params["molecule"]+"_mol{}".format(imol),
            "sms", q[-1], params["NradAzmBins"]))
        with open(fName, "wb") as file:
          (q*mol_diffraction[imol]/atm_diffraction[imol]).astype(np.double).tofile(file)



    fName_template_h5 =\
        "{0}_sim_diffraction-azmAvg_Qmax-{1:.4g}.h5"
    fName = os.path.join(params["simOutputDir"], fName_template_h5.format(
        params["molecule"], q[-1]))
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
        "{}_sim-diffraction-azmAvg.png".format(params["molecule"])))
    """

  else:
    if "sim_eval_times" in params:
      params["ADM_params"]["eval_times"] = params["sim_eval_times"]
    LMK, bases, norms, times  = get_ADMs(params["ADM_params"])#, logging, plot=(not on_cluster))
    #params, bases, times      = apply_fit_parameters(params, bases, times, logging)
    logging.info(LMK)

    # Set up pixel and s maps
    x,z = x_lab,z_lab
    zero_mask = pixel_dist > 0

    signal_xyz_qf = np.concatenate([
        np.expand_dims(x, axis=-1),
        np.expand_dims(np.ones_like(x)*params["beamline_length"], axis=-1),
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

    atm_diffraction = atomic_diffraction_calculation(
        q_map_polar, atom_types, scat_amps, detector_dist)

    print(atm_diffraction)
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


    if params["calculation_type"] == "numeric":
      ###  Ensemble smearing  ###
      smearing_polar  = np.linspace(0, np.pi, params["smear_polar_bins"]+1,
          dtype=np.float32)
      smearing_azim   = np.linspace(0, 2*np.pi, params["smear_azim_bins"]+1,
          dtype=np.float32)
      smearing_spin   = np.linspace(0, 2*np.pi, params["smear_spin_bins"]+1,
          dtype=np.float32)
      
      """
      smearing_polar = np.arange(params["smear_polar_bins"]+1)\
          *np.pi/(params["smear_polar_bins"])
      smearing_azim = np.arange(params["smear_azim_bins"]+1)\
          *2*np.pi/(params["smear_azim_bins"])
      smearing_spin = np.arange(params["smear_spin_bins"]+1)\
          *2*np.pi/(params["smear_spin_bins"])
      #smearing_spin = np.concatenate([smearing_spin, smearing_spin[1:]+np.pi])
      """

      # Use midpoint for numerical integration
      smearing_polar  = (smearing_polar[1:] + smearing_polar[:-1])/2.
      smearing_azim   = (smearing_azim[1:] + smearing_azim[:-1])/2.
      smearing_spin   = (smearing_spin[1:] + smearing_spin[:-1])/2.

      """
      smearing_polar[params["smear_polar_bins"]//2:] =\
          smearing_polar[:params["smear_polar_bins"]//2]
      smearing_azim[params["smear_azim_bins"]//2:] =\
          smearing_azim[:params["smear_azim_bins"]//2]
      smearing_spin[params["smear_spin_bins"]//2:] =\
          smearing_spin[:params["smear_spin_bins"]//2]
      """
      print("VAR DELT", np.var(smearing_spin[1:] - smearing_spin[:-1]))

      #smearing_spin = np.array([4*np.pi/3])
      #print(smearing_spin)
      """
      smearing_polar  = np.arange(params["smear_polar_bins"])\
          *np.pi/params["smear_polar_bins"]
      smearing_azim   = np.arange(params["smear_azim_bins"])\
          *2*np.pi/params["smear_azim_bins"]
      smearing_spin   = np.arange(params["smear_azim_bins"])\
          *2*np.pi/params["smear_azim_bins"]
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
      #smearing_ea = np.array([np.zeros(params["smear_azim_bins"]), smearing_azim, np.zeros(params["smear_azim_bins"])]).transpose()
      print("SSSSSS", smearing_ea.shape)
      smearing_jacobian = np.sin(smearing_ea[:,1])
      smearing_jacobian *= np.pi/params["smear_polar_bins"]
      smearing_jacobian *= 2*np.pi/params["smear_azim_bins"]
      smearing_jacobian *= 2*np.pi/params["smear_spin_bins"]
      

      for tm in range(bases.shape[-1]):
        if time_ind is not None:
          if tm != time_ind:
            continue

        mol_diffraction, mol_diffraction_raw = diffraction_calculation_numeric(
            params["molecule"], tm, times, LMK, bases.astype(np.float32),
            (params["imgBins"], params["imgBins"]),
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
        fName = os.path.join(params["simOutputDir"], fName_template_h5.format(
            params["molecule"], q_map[params["imgBins"]//2,-1],
            "_time-{0:.6g}".format(times[tm])))
        if params["input_LMK"] is not None:
          suffix = "_LMK"
          for lmk in params["input_LMK"]:
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


    elif params["calculation_type"] == "analytic":
      if params["QperPix"] is None:
        logging.fatal("Must specify QperPix for the analytic option.")
        sys.exit(0)
    
      #if params["molecule"] == "N2O":
      #  bases *= np.expand_dims(np.sqrt(4*np.pi/(2*LMK[:,0] + 1)), -1)**2
      dBrog = calculate_deBroglie_wavelength(params)
      for tm in range(bases.shape[-1]):
        if time_ind is not None:
          if tm != time_ind:
            continue

        print("INP SHAPPPPPPPPPPPPPPPPPPPp", bases.shape)
        print(atom_positions)
        mol_diffraction, mol_diffraction_raw = diffraction_calculation_analytic(
            LMK, bases.astype(np.float32)[:,tm],
            atom_positions, atom_types, scat_amps,
            q_map_polar, detector_dist=detector_dist, freq=dBrog)
        print("diff test", np.all(mol_diffraction == atm_diffraction), np.amax(np.abs(mol_diffraction - atm_diffraction)))

        rng = np.max(np.abs([np.amax(mol_diffraction), np.amin(mol_diffraction)]))
        #plt.pcolormesh(mol_diffraction, vmax=rng, vmin=-1*rng, cmap='seismic')
        #plt.colorbar()
        #plt.savefig("aresult.png")
        fName_template_h5 =\
            "{0}_sim_diffraction-analytic_Qmax-{1:.4g}{2}"
        fName = os.path.join(params["simOutputDir"], fName_template_h5.format(
            params["molecule"], q_map[params["imgBins"]//2,-1],
            "_time-{0:.6g}".format(times[tm])))
        if params["input_LMK"] is not None:
          suffix = "_LMK"
          for lmk in params["input_LMK"]:
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
  if hasattr(params, '__dict__'):
    params = params.__dict__

  if FLAGS.molecule is not None:
    params["molecule"] = FLAGS.molecule
  if FLAGS.output_folder is not None:
    params["simOutputDir"] = FLAGS.output_folder
  FLAGS = setup(parser, output_dir=params["simOutputDir"])

  print("IM HERE")
  # Flag handling
  if FLAGS.calculation_type is not None:
    params["calculation_type"] = FLAGS.calculation_type
  if FLAGS.xyz_file is not None:
    params["xyz_file"] = FLAGS.xyz_file
  if FLAGS.output_fileName_suffix is not None:
    params["output_fileName_suffix"] = FLAGS.output_fileName_suffix
  if FLAGS.basis_folder is not None:
    params["basis_folder"] = FLAGS.basis_folder
  if FLAGS.LMK is not None:
    params["input_LMK"] = FLAGS.LMK
  else:
    params["input_LMK"] = None
  if FLAGS.time_ind is not None:
    if FLAGS.time_ind == "nan":
      FLAGS.time_ind = None
    else:
      FLAGS.time_ind = int(FLAGS.time_ind)
  if FLAGS.eval_time is not None:
    params["sim_eval_times"] = np.array([float(FLAGS.eval_time)])
    FLAGS.time_ind = 0
  if "xyz_dir" not in params:
    params["xyz_dir"] = "XYZ"

  params["imgBins"] = 2*params["NradAzmBins"] - 1 
  
  if params["calculation_type"] != "azmAvg":
    from ADM import get_ADMs

  print("FFFFFFFFFIIIIIIIIIIIIIXXXXXXXXXXXXXX ME: Change from q_polar to freq calculation of theta/phi lf.")
  main(params, on_cluster=FLAGS.cluster, time_ind=FLAGS.time_ind)
  

