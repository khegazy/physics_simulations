import os
import numpy as np
home_dir = "/cds/home/k/khegazy/"

def get_parameters(run, molecule="N2O"):


  class parameters: 

    def __init__(self, run, molecule):
      self.molName     = molecule
      self.run         = run
      self.data_format = 'ts'
      self.elEnergy    = 3.7e6
      self.beamline_length = 4
      self.probe_FWHM  = 180


      self.QperPix         = 0.028156*(137.0/112.0)*1.03
      self.NmaxRadBins     = 750
      self.NradAzmBins     = 375 

      self.xyz_file            = "N2O.xyz"
      
      self.xyzDir              = "./XYZfiles/"
      self.scat_amps_dir       = home_dir+"simulation/scatteringAmplitudes/3.7MeV/"
      self.axisDistDir         = "/cds/group/ued/scratch/N2O/axis_distributions/" + self.molName + "/"
      self.simOutputDir        = "./output"
  
      if run == "sim_validation":
        self.molName  = "CF2IBr"
        self.xyz_file = "CF2IBr.xyz"
        self.NradAzmBins = 75
        self.QperPix = 20.0/(self.NradAzmBins-1)
        N = 1
        self.sim_eval_times = np.arange(N)
        #self.sim_times = np.arange(5)*0.25 + 318
        self.calculation_type = "analytic"
        self.xyz_file = os.path.join(self.xyzDir, self.xyz_file)
        self.basis_folder = os.path.join(self.axisDistDir, "A/temp-100K")
        self.smear_polar_bins = 5 
        self.smear_azim_bins  = 5
        self.smear_spin_bins  = 5

  return parameters(run, molecule)
