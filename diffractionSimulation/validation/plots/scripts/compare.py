import numpy as np
import h5py
import matplotlib.pyplot as plt

for i in range(1):
  #fName = "output/NO2_sim_diffraction-azmAvg_align-random_Qmax-13.27time-{}.h5".format(i)
  """
  fName = "output/test.h5"
  with h5py.File(fName, 'r') as hf:
    q = hf["q"][:]
    ad = hf["atm_diffraction"][:]
    mdn = hf["mol_diffractionN"][:]
    mda = hf["mol_diffractionA"][:]
  """
  folder = "/cds/group/ued/scratch/N2O/simulations/"
  mol = "test"
  tm="39.5"#"318.5"
  tm = "39.7626"#"39.1162"
  tm = "0"
  sq="20"
  fName = folder+mol+"_sim_diffraction-numeric_Qmax-"+sq+"_time-"+tm+".h5"
  with h5py.File(fName, 'r') as hf:
    q = hf["q"][:]
    ad = hf["atm_diffraction"][:]
    mdn = hf["mol_diffraction"][:]
  
  fName = folder+mol+"_sim_diffraction-analytic_Qmax-"+sq+"_time-"+tm+".h5"
  with h5py.File(fName, 'r') as hf:
    q = hf["q"][:]
    ad = hf["atm_diffraction"][:]
    mda = hf["mol_diffraction"][:]

  #print("mean", np.mean(mda/mdn))
  #mdn /= ad
  #mda /= ad
  #mdn /= np.sqrt(np.sum(mdn**2))
  #mda /= np.sqrt(np.sum(mda**2))
  #mdn *= q
  #mda *= q

  sh = ad.shape
  #print("NORMS", q[sh[0]//2, sh[1]//2], ad[sh[0]//2, sh[1]//2], md[sh[0]//2, sh[1]//2])

  fig, ax = plt.subplots()
  a = ax.pcolormesh(ad)
  fig.colorbar(a, ax=ax)
  ax.set_aspect('equal')
  fig.savefig("../atm_diff_{}.png".format(i))
  plt.close()
 
  fig, ax = plt.subplots()
  a = ax.pcolormesh(mdn)
  fig.colorbar(a, ax=ax)
  ax.set_aspect('equal')
  #print("NOT SAVING")
  fig.savefig("../mol_diff_numerical_{}.png".format(i))
  plt.close()
  
  fig, ax = plt.subplots()
  a = ax.pcolormesh(mda)
  fig.colorbar(a, ax=ax)
  ax.set_aspect('equal')
  fig.savefig("../mol_diff_analytic_{}.png".format(i))
  plt.close()
  
  sys.exit(0)
  fig, ax = plt.subplots()
  a = ax.pcolormesh(mdn - mda)
  fig.colorbar(a, ax=ax)
  ax.set_aspect('equal')
  fig.savefig("../mol_diff_{}.png".format(i))
  plt.close()
  
  fig, ax = plt.subplots()
  rng = np.max(np.abs([np.amax(mda), np.amin(mda)]))
  a = ax.pcolormesh(mda, vmin=-1*rng, vmax=rng, cmap='seismic')
  ax.xaxis.set_visible(False)
  ax.yaxis.set_visible(False)
  
  fig.colorbar(a, ax=ax)
  ax.set_aspect('equal')
  fig.savefig("../mol_diff_analytic_4fig.png".format(i))
  plt.close()

