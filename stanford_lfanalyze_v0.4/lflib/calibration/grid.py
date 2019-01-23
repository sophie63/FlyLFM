import numpy as np
from collections import deque

# Build a 2D grid from an unlabeled set of lenslet centers.  This
# method is ad-hoc, but seems to work fairly reliably.
#
# Here we iterate over the putative matches, looking for pixels that
# at within 2 pixels of the EXPECTED_LENSLET_SIZE.  This seems to be a
# very reliable way of finding neighboring lenslets.  Note that we
# throw away any pixels here that don't have exactly four matches.
#
def build_grid(putative_centers, expected_lenslet_size, adjacency_tolerance, b_estimate_rotation = True):

    # Optionally, first estimate the rotation of the lenslet grid using random sampling (+-45 degrees)
    lenslet_rotation = 0
    ROTATION_SAMPLE_SIZE = 400
    if b_estimate_rotation:
        import random
        if len(putative_centers) < ROTATION_SAMPLE_SIZE:
            print "ERROR: there are not enough lenslet centers to perform calibration.  Check you calibration image and try again."
            raise SystemExit
        rand_centers = random.sample(putative_centers, ROTATION_SAMPLE_SIZE)
        rotation_samples = []
        for rc in rand_centers:
            diff = np.cast['float32'](putative_centers-np.tile(rc,
                                                               (putative_centers.shape[0],1)))
            dist = np.abs(np.sqrt(np.sum(diff**2,1))-expected_lenslet_size)
            adj = np.nonzero(dist < adjacency_tolerance)[0]
            for a in adj:
                dr = putative_centers[a][0] - rc[0]
                dc = putative_centers[a][1] - rc[1]
                theta = np.arctan2(dc,dr)
                theta = ((theta + (5*np.pi/4)) % (np.pi/2)) - (np.pi/4)
                rotation_samples.append(theta)
        lenslet_rotation = np.median(rotation_samples)
        print( '\t--> Estimating lenslet rotation (n=%d): %f degrees' % 
               (len(rotation_samples),
                np.rad2deg(lenslet_rotation)) )
    else:
        print( '\t--> Assuming lenslet basis aligned with sensor.')      

    # First we find the lenslet closest to the center.
    num_putative_lenslets = putative_centers.shape[0]
    center_of_mass = np.mean(putative_centers,0)
    diff = putative_centers - np.tile(center_of_mass, (num_putative_lenslets,1))
    dist = np.abs(np.sqrt(np.sum(diff**2,1)))
    center_lenslet = np.argmin(dist)
    print ( '\t--> Found center lenslet %d near pixel [ %d %d ]' %
            (center_lenslet,
             putative_centers[center_lenslet,1],
             putative_centers[center_lenslet,0]) )

    # Pull off the first entry to get things started and mark it as touched
    touched = np.zeros(putative_centers.shape[0])
    pending = deque( [(putative_centers[center_lenslet,:], 0, 0)] )
    touched[center_lenslet] = 1

    # Walk the graph, numbering lenslets as they are reached from
    # the central node.  Our goals here is to give each lenslet a
    # unique 2D coordinate on a 2D grid, where the center lenslet
    # is (0,0).  This is also fairly efficient, since we rapidly
    # narrow down the number of lenslets we have to search through
    # as the algorithm proceeds.
    iter_count = 0
    lenslets = []
    while (len(pending) > 0):
        current_match = pending.popleft()
        current_center = current_match[0]
        r = current_match[1]
        c = current_match[2]

        # OPTIMIZATION: We periodically prune the putative_centers
        # list and reset the touched list, since this
        # significantly increases the speed of matching towards
        # the end.
        if (iter_count % 300 == 0):
            untouched_entries = np.nonzero(touched == 0)[0]
            putative_centers = putative_centers[untouched_entries,:]
            touched = np.zeros(putative_centers.shape[0])

        if (iter_count % 1000 == 0):
            print ('\t    Iteration %d at [%d, %d]  ( %d pending, %d left in match pool )' %
                   (iter_count, c, r, len(pending), len(putative_centers)))
        iter_count = iter_count + 1

        # Appending this lenslet to the master list along with its location in the grid
        lenslets.append((current_center[0], current_center[1], r, c))

        # Now find the neighbors of this lenslet, and add them to the pending queue.
        diff = np.cast['float32'](putative_centers-np.tile(current_center,
                                                           (putative_centers.shape[0],1)))
        dist = np.abs(np.sqrt(np.sum(diff**2,1))-expected_lenslet_size)
        matches = np.nonzero(dist < adjacency_tolerance)[0]
        for m in matches:
            if not touched[m]:
                touched[m] = 1
                dr = putative_centers[m][0] - current_center[0]
                dc = putative_centers[m][1] - current_center[1]
                dr_rot = dr*np.cos(-lenslet_rotation) - dc*np.sin(-lenslet_rotation)
                dc_rot = dr*np.sin(-lenslet_rotation) + dc*np.cos(-lenslet_rotation)
                #print dr, dc, dr_rot, dc_rot
                if np.abs(dc_rot) < adjacency_tolerance:
                    if dr_rot > 0:
                        # Up (row+1, col)
                        pending.append((putative_centers[m],r+1,c))
                    else:
                        # Down (row-1, col)
                        pending.append((putative_centers[m],r-1,c))
                elif np.abs(dr_rot) < adjacency_tolerance:
                    if dc_rot > 0:
                        # Right (row, col+1)
                        pending.append((putative_centers[m],r,c+1))
                    else:
                        # Left (row, col-1)
                        pending.append((putative_centers[m],r,c-1))
                else:
                    # Silently discard any lenslets that aren't
                    # either immediately above or below this
                    # lenslet (to within a 2 pixel margin).
                    pass

    # Pick one of the lenslets and compute roughly what lenslet it is
    # in the original image.  We use this information to recenter the
    # lenslet transform.  This is not perfect, btw, but does roughly
    # center the rectified image where the old image was.
    l = lenslets[0]
    import math
    minr = -math.floor(l[0] / expected_lenslet_size)
    minc = -math.floor(l[1] / expected_lenslet_size)

    # Recenter the lenslet coordinates.  Note that we add an
    # additional 0.5 lenslet offset so that the lenslet centers
    # themselves appear at [0.5,0.5], [0.5, 1.5], ... etc.
    #
    print '\t--> Recentering lenslets by [%d, %d].' % (minc, minr)
    adjusted_lenslets = np.zeros((len(lenslets),4))
    m = 0
    for l in lenslets:
        adjusted_lenslets[m,:] = np.array([l[0], l[1], float(l[2]-minr-0.5), float(l[3]-minc-0.5)])
        m += 1

    return adjusted_lenslets
