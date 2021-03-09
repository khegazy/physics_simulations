import sys, os, glob
import pickle as pl
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from collections import defaultdict
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d


def get_bases(params, logging=None, normalize=False, plot=False):
  bases = {}
  norms = {}
  LMK = []

  """
  if params.get_bases_type is 0:
    files = "{}/*.dat".format(folderName)

    # Get array of which bases to gather
    bases_Linds = fit_lgInds
    if FLAGS.setBases is not None:
      bases_Linds = np.unique(np.expand_dims(fit_lgInds, 0)\
          + np.expand_dims(FLAGS.setBases, -1).flatten())
      bases_Linds = bases_Linds[bases_Linds >=0]

    basisList = [None for x in range(len(bases_Linds))]
    normsList = [None for x in range(len(bases_Linds))]
    # Gather bases
    for fl in glob.glob(files):
      L = int(fl[fl.find("_L-")+3:fl.find("_M-")])
      ind = np.where(bases_Linds == L)[0]

      # Is L in the index list
      if len(ind) != 1:
        continue
      ind = ind[0]

      # Get the basis and normalize
      with open(fl, "rb") as file:
        if time_cuts is None:
          basisList[ind] = np.fromfile(file, np.double)
        else:
          basisList[ind] = np.fromfile(file, np.double)\
              [time_cuts[0]:time_cuts[1]]

        if subtract_mean:
          if L != 0:
            basisList[ind] -= np.mean(basisList[ind])
        #print("norm",L,np.sqrt(np.sum(bases[L]**2)), np.amax(np.abs(bases[L])))
        normsList[ind] = np.sqrt(np.sum(basisList[ind]**2))
        if normalize:
          basisList[ind] /= normsList[ind]

    allBases = np.array(basisList)
    allNorms = np.array(normsList)
    setBases = None
  """

  if params.molName == "N2O":
    start_time, end_time = None, None
    files = "{}/*.dat".format(params.basis_folder)
    print("files", files)
    for fl in glob.glob(files):
      if logging is not None:
        logging.info("Importing Bases: " + fl)
      L = int(fl[fl.find("_L-")+3:fl.find("_M-")])
      with open(fl, "rb") as file:
        bases[L] = np.fromfile(file, np.double)

        """
        if normalize:
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
        """

        # Analytic calculation expects Wigner D basis
        #if params.calculation_type == "analytic":
          #bases[L] *= np.expand_dims(np.sqrt(4*np.pi/(2*L + 1)), -1)**2
        bases[L] /= np.expand_dims(
            ((2*L + 1)/(8*np.pi**2))*np.sqrt(4*np.pi/(2*L+1)), -1)
        #if params.calculation_type == "analytic":
          #bases[L] *= np.expand_dims(np.sqrt(4*np.pi/(2*L + 1)), -1)**2

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
            if logging is not None:
              logging.critical("Bases do not have the same time range.")
            sys.exit(1)

    basisList = []
    #normsList = []
    keys = np.sort(list(bases.keys()))
    for l in keys:
      basisList.append(np.expand_dims(bases[l], axis=0))
      #normsList.append(norms[l])
      LMK.append(np.array([l,0,0]))
    tLen = basisList[0].shape[-1]
    times = start_time + np.linspace(0, 1, tLen)*(end_time - start_time)
    print("TIMES", times[0], times[-1])

  elif params.molName == "NO2":
    L_list = []
    M_dict = defaultdict(list)
    files = "{}/*.npy".format(params.basis_folder)
    print("search", files)
    for fl in glob.glob(files):
      if logging is not None:
        logging.info("Importing Bases: " + fl)
      D = fl[fl.find(" D")+2:-4]
      ln = len(D) + len(D)%2
      L = int(D[:ln//2])
      M = int(D[ln//2:])
      print(D,L,M)
      if L not in L_list:
        L_list.append(L)
      M_dict[L].append(M)
      with open(fl, "rb") as file:
        bases[(L,M)] = np.load(file)
      if M != 0:
        bases[(L,M)] *= 2

      """
      if normalize:
        # Normalize bases
        bases_womean = {}
        if L != 0:
          bases_womean[(L,M)] = bases[(L,M)] - np.mean(bases[(L,M)])
        else:
          bases_womean[L] = bases[L]
        #print("norm",L,np.sqrt(np.sum(bases[L]**2)), np.amax(np.abs(bases[L])))
        norms[L] = np.sqrt(np.sum(bases_womean[L]**2))
        if normalize:
          bases[L] = bases_womean[L]
          bases[L] /= norms[L]
      """

    basisList = []
    normsList = []
    L_list.sort()
    for l in L_list:
      for m in np.sort(M_dict[l]):
        basisList.append(np.expand_dims(bases[(l,m)], axis=0))
        LMK.append(np.array([l,0,m]))

    fName = params.basis_folder.replace("/A/", "/times/").replace("temp-", "temp_")\
        + ".npy"
    if logging is not None:
      logging.info("Time file: " + fName)
    with open(fName, "rb") as file:
      times = np.load(file)

  elif params.input_LMK is not None:
    basisList = []
    print("WAT", params.input_LMK)
    for i in range(len(params.input_LMK)):
      LMK.append(params.input_LMK[i,:])
      basisList.append(np.ones((1,1)))
    times = np.array([0])
  else:
    basisList = []
    for l in range(2,3):
      for k in range(-l,l+1):
        if k != 2:
          continue
        LMK.append(np.array([l,0,k]))
        basisList.append(np.complex(0,1)*np.ones((1,1)))

    times = np.array([0])

  bases = np.concatenate(basisList, axis=0)

  # Create evenly spaced values
  if bases.shape[-1] > 1:
    bases_interp = interp1d(times, bases, axis=-1, kind='cubic')
    times = np.linspace(times[0], times[-1], 1000)
    bases = bases_interp(times)

    if plot:
      inds = np.arange(len(times))# < 100
      for i, (L, M, K) in enumerate(LMK):
        plt.plot(times[inds], bases[i,inds])
        plt.savefig(os.path.join("./plots", "basis_input_{}-{}-{}.png".format(
            L, M, K)))
        plt.close()

    if params.probe_FWHM is not None:
      print("INTERPING")
      delta_time = times[1] - times[0]
      bases = gaussian_filter1d(
          bases, (0.001*params.probe_FWHM/2.355)/delta_time, axis=-1)

    # Select evaluation times 
    bases_interp = interp1d(times, bases, axis=-1, kind='cubic')
    times = params.sim_eval_times
    bases = bases_interp(times)

    if plot:
      inds = np.arange(len(times))# < 100
      for i, (L, M, K) in enumerate(LMK):
        plt.plot(times[inds], bases[i,inds])
        plt.savefig(os.path.join("./plots", "basis_smeared_{}-{}-{}.png".format(
            L, M, K)))
        plt.close()

  print("BASIS LMK", LMK)
  return np.array(LMK), bases, [], times



