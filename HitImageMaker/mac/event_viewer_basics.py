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
cfg_mgr.AddCfgFile("../mac/ClusterBasicsPSet.fcl")

config = cfg_mgr.Config()
pset = config.get_pset("ClusterBasics");

my_proc = fmwk.ana_processor()
my_proc.add_input_file(sys.argv[1])

my_proc.set_io_mode(fmwk.storage_manager.kBOTH)
my_proc.set_output_file("deleteme.root")
my_proc.set_data_to_write(fmwk.data.kHit,'gaushit')
my_proc.set_data_to_write(fmwk.data.kCluster,'hullcluster')
my_proc.set_data_to_write(fmwk.data.kAssociation,'hullcluster')

my_unit         = fmwk.HitImageManager()
my_imagemaker   = fmwk.AlgoImageMaker()
my_clusterbasics  = fmwk.AlgoClusterBasics(pset)

my_unit.SetImageMaker(my_imagemaker)
my_unit.SetClusterMaker(my_clusterbasics)
my_proc.add_process(my_unit)





def getview(plane,ax) :

    hulls = my_clusterbasics.hulls(plane)
    for x in xrange(hulls.size()):
        hull = hulls[x]

        x_array = []
        y_array = []
        for p in xrange(hull.size()):
            x_array.append(hull[p].first)
            y_array.append(hull[p].second)

        x_array.append(x_array[0])
        y_array.append(y_array[0])
        
        ax.plot(x_array,y_array,color='black',lw=2)

    pts = my_clusterbasics.other_hits(plane)
    xx= []
    yy= []
    for x in xrange(pts.size()):
        p = pts[x]

        xx.append(p.first)
        yy.append(p.second)

    ax.plot(xx,yy,'o',markersize=1)

    pclusters = my_clusterbasics.p_clusters(plane)

    xx= []
    yy= []
    for x in xrange(pclusters.size()):
        pcluster = pclusters[x]

        hits = pcluster.hits()
        xx= []
        yy= []

        for i in xrange(hits.size()):
            hit = hits[i]
            xx.append(hit.first)
            yy.append(hit.second)

        ax.plot(xx,yy,'o',markersize=3,c=np.random.rand(3,1))

    return ax



while ( my_proc.process_event() ):
    fig = plt.subplots(figsize=(10,6))
    ax = plt.subplot(3,1,1)
    ax = getview(0,ax)
    ax = plt.subplot(3,1,2)
    ax = getview(1,ax)
    ax = plt.subplot(3,1,3)
    ax = getview(2,ax)

    plt.show()

