#%matplotlib inline

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image  as mpimg
import matplotlib.cm     as cm
import matplotlib.collections as mc
import numpy as np
import cv2

from larlite import larlite as fmwk
from ROOT import fcllite

cfg_mgr = fcllite.ConfigManager()
cfg_mgr.AddCfgFile("../mac/ClusterHullPSet.fcl")

config = cfg_mgr.Config()
pset = config.get_pset("ClusterHull");

my_proc = fmwk.ana_processor()
my_proc.add_input_file("/Users/vgenty/Desktop/single_e.root")

my_proc.set_io_mode(fmwk.storage_manager.kBOTH)
my_proc.set_output_file("deleteme.root")
my_proc.set_data_to_write(fmwk.data.kHit,'gaushit')
my_proc.set_data_to_write(fmwk.data.kCluster,'hullcluster')
my_proc.set_data_to_write(fmwk.data.kAssociation,'hullcluster')

my_unit         = fmwk.HitImageManager()
my_imagemaker   = fmwk.AlgoImageMaker()
my_clusterhull  = fmwk.AlgoClusterHull(pset)

my_unit.SetImageMaker(my_imagemaker)
my_unit.SetClusterMaker(my_clusterhull)
my_proc.add_process(my_unit)





def getview(plane,ax) :
    
    houghs = my_clusterhull.houghs(plane)
    hh=[]
    for h in xrange(houghs.size()):
        l = houghs[h]
        hh.append([(l[0],l[1]),(l[2],l[3])])

    cc = mc.LineCollection(hh, linewidths=1)

    hulls = my_clusterhull.hulls(plane)
    for x in xrange(hulls.size()):
        hull = hulls[x]

        x_array = []
        y_array = []
        for p in xrange(hull.size()):
            x_array.append(hull[p].first)
            y_array.append(hull[p].second)

        ax.plot(x_array,y_array,color='black',lw=2)

    pts = my_clusterhull.other_hits(plane)
    xx= []
    yy= []
    for x in xrange(pts.size()):
        p = pts[x]

        xx.append(p.first)
        yy.append(p.second)

    ax.plot(xx,yy,'o',markersize=1)

    pclusters = my_clusterhull.p_clusters(plane)

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

    ax.add_collection(cc)

    return x



while ( my_proc.process_event() ):
    fig = plt.subplots(figsize=(10,6))
    ax = plt.subplot(3,1,1)
    ax = getview(0,ax)
    ax = plt.subplot(3,1,2)
    ax = getview(1,ax)
    ax = plt.subplot(3,1,3)
    ax = getview(2,ax)

    plt.show()

