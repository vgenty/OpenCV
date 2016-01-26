#%matplotlib inline

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image  as mpimg
import matplotlib.cm     as cm
import matplotlib.collections as mc
import numpy as np
import cv2

import sys

from larlite import larlite as fmwk
from ROOT import fcllite

cfg_mgr = fcllite.ConfigManager()
cfg_mgr.AddCfgFile("ClusterBasicsPSet.fcl")

config = cfg_mgr.Config()
pset = config.get_pset("ClusterBasics");

my_proc = fmwk.ana_processor()

for x in xrange(len(sys.argv)-1):    
    my_proc.add_input_file(sys.argv[x+1])

my_proc.set_io_mode(fmwk.storage_manager.kBOTH)
my_proc.set_output_file("deleteme_basics_mu.root")
my_proc.set_data_to_write(fmwk.data.kHit,'gaushit')
my_proc.set_data_to_write(fmwk.data.kCluster,'hullcluster')
my_proc.set_data_to_write(fmwk.data.kAssociation,'hullcluster')

my_unit         = fmwk.HitImageManager()
my_imagemaker   = fmwk.AlgoImageMaker()
my_clusterbasics  = fmwk.AlgoClusterBasics(pset)

my_unit.SetImageMaker(my_imagemaker)
my_unit.SetClusterMaker(my_clusterbasics)
my_proc.add_process(my_unit)
my_proc.run()


