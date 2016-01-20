import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image  as mpimg
import matplotlib.cm     as cm
import matplotlib.collections as mc
import numpy as np
import cv2
import sys

from larlite import larlite as fmwk
storage = fmwk.storage_manager()

filename = str(sys.argv[1])

storage.add_in_filename(filename)

storage.set_io_mode(storage.kREAD)
storage.open()

while ( storage.next_event() ):


    clusters = storage.get_data(fmwk.data.kCluster,"hullcluster")
    hits     = storage.get_data(fmwk.data.kHit,"gaushit")
    ass_data = storage.get_data(fmwk.data.kAssociation,"hullcluster")


    cluster_to_hit_ass = ass_data.association(clusters.id(),hits.id());
    

    fig, ax = plt.subplots(figsize=(10,6))

    for c in xrange(cluster_to_hit_ass.size()) :
        xx= []; yy= [];
        for h in xrange(cluster_to_hit_ass[c].size()):
            hit_idx = cluster_to_hit_ass[c][h]
            xx.append(hits[hit_idx].WireID().Wire); yy.append(hits[hit_idx].PeakTime())
        
        ax.plot(xx,yy,'.')
        
    plt.show()
