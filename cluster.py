#####################################################
# Imports
#####################################################

import glob
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from mdtraj import xtc
import rmsd

#####################################################
# Globals
#####################################################
# this parameter controls how many trajectories will
# be loaded. The computation is O(N^2), so loading
# trajectories will slow it down substantially

MAX_TRAJECTORIES = 10

#####################################################
# Script
#####################################################

print 'Loading all of the files...'
files = glob.glob('xtc/RUN*/CLONE*/*.xtc')
xyz = []
for i,f in enumerate(files):
    xyz.append(xtc.read(f).xyz)
    print 'Loading trajectory %s' % i

    if i >= MAX_TRAJECTORIES:
        break

print 'Frames loaded successfully'

# concatenate all of the data together into a
# single array array
xyz = np.concatenate(xyz)
print 'Number of frames of data loaded %s' % len(xyz)

# compute the full pairwise distance matrix
# between all of the observed conformations
# using the rmsd_qcp algorithm
print 'Computing the pairwise distance matrix'
distances = np.empty(len(xyz)*(len(xyz)-1)/2)
k = 0
for i in xrange(len(xyz)):
    for j in xrange(i+1, len(xyz)):
        distances[k] = rmsd.rmsd_qcp(xyz[j], xyz[j])
        k += 1
    print 'RMSD computation %s/%s' % (k, len(xyz)*(len(xyz)-1)/2)


print 'mean pairwise distance  ', np.mean(distances)
print 'stddev pairwise distance', np.std(distances)

print 'Running hierarchical clustering (UPGMA)...'
# run hierarchical clustering on the distance matrix
Z = linkage(distances, method='average')

# get flat clusters from the linkage matrix corresponding
# to states
print 'Flattening the clusters...'
assignments = fcluster(Z, t=5000, criterion='maxclust')
print 'Number of clusters found', len(np.unique(assignments))
