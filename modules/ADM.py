import sys, os, glob
import pickle as pl
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d


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

    fName = FLAGS.basis_folder.replace("/A/", "/times/").replace("temp-", "temp_"    )\
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

