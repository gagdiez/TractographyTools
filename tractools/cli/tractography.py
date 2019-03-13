''' Given a set of seeds, computes streamlines for each seeds and returns
    a visits' map over a mask '''
import itertools
from functools import partial
import os
import multiprocessing
import logging

import numpy

from dipy.direction import ProbabilisticDirectionGetter as probabilistic
from dipy.direction import DeterministicMaximumDirectionGetter as deterministic

from .. import utils

def tracking(shm_file, mask_file, outdir, force_overwrite,
             particles, step_size, max_lenght, max_angle, algorithm,
             wpid_seeds_info):
    ''' Tracking function that will run in parallel

        Params:
            shm_file: SHM file computed from the dwi file
            mask_file: mask were to perform tractography
            outdir: Directory were to save streamlines
            force_overwrite: if True, existing files will be overwriten
            step_size: size in mm of each step in the tracking
            max_lenght: maximum lenght of each streamline
            max_angle: maximum angle at each step of tracking
            algoright: either 'probabilistic' or 'deterministic'
            wpid_seeds_info: tuple which contains:
                - wpid: The id of this worker
                - seeds: One list for each seed with points to track from
                - info: CIFTI information for each seed:
                    -mtype: A valid CIFTI MODELTYPE
                    -name: A valid CIFTI BRAINSTRUCTURE
                    -coord: Voxel or vertex to which the seed makes reference
                    -size: size of the CIFTI SURFACE (if applies)
        Returns:
            list of streamlines '''
    import citrix
    import streamlines as sl

    from dipy.data import default_sphere
    from dipy.tracking.local import LocalTracking
    from dipy.tracking.local import BinaryTissueClassifier

    wpid, (seeds, cifti_info) = wpid_seeds_info

    logging.debug("Worker {} started".format(wpid))

    # Check if file exists
    outfile = os.path.join(outdir, "stream_{}.trk".format(wpid))

    if os.path.isfile(outfile) and not force_overwrite:
        print("File already exists, use the -f flag to overwrite it")
        return

    shm = citrix.load(shm_file)
    shm_data = shm.get_data()

    mask_nib = citrix.load(mask_file)
    mask = mask_nib.get_data()

    if algorithm == 'deterministic':
        directions = deterministic.from_shcoeff(shm_data, max_angle,
                                                default_sphere)
    else:
        directions = probabilistic.from_shcoeff(shm_data, max_angle,
                                                default_sphere)

    classifier = BinaryTissueClassifier(mask)

    percent = max(1, len(seeds)/5)
    streamlines = []
    used_seeds = []

    for i, s in enumerate(seeds):
        if i % percent == 0:
            logging.debug("{}, {}/{} seeds".format(wpid, i, len(seeds)))

        # Repeat the seeds as long as needed
        if len(s) == 3:
            # It's one point
            s = [s]
        repeated_seeds = itertools.cycle(s)

        res = LocalTracking(directions, classifier, repeated_seeds, shm.affine,
                            step_size=step_size, maxlen=max_lenght)
        it = res._generate_streamlines()  # This is way faster, just remember
                                          #  after to move them into mm space
                                          #  again.
        for streamline in itertools.islice(it, particles*len(s)):
            if streamline is not None and len(streamline) > 1:
                streamlines.append(streamline)

            if cifti_info[i][0] == 'CIFTI_MODEL_TYPE_SURFACE':
                used_seeds.append(cifti_info[i][2])
            else:
                used_seeds.append([int(cf) for cf in cifti_info[i][2]])

    streamlines = sl.Streamlines(streamlines, affine=shm.affine)

    numpy.savetxt(os.path.join(outdir, "info_{}.txt".format(wpid)),
                  used_seeds)
    
    sl.io.save(streamlines, outfile, shm.shape[:3], shm.header.get_zooms()[:3])

    logging.debug("Worker {} finished".format(wpid))
    return


def tractography(shm_file, mask_file, seeds_file, outdir,
                 algorithm='probabilistic', particles=5000, step_size=1,
                 max_lenght=200, max_angle=30,
                 nbr_process=0, spp=None, verbose=0, force_overwrite=False):
    if verbose:
        logging.basicConfig(level=logging.DEBUG)

    #Start multiprocessing environment
    if not nbr_process:
        nbr_process = multiprocessing.cpu_count()

    cifti_info, seeds_pnts = utils.seeds.load(seeds_file)
    s = 500 if spp is None else spp
    seed_chunks = [seeds_pnts[i:i+s] for i in range(0, len(seeds_pnts), s)]
    info_chunks = [cifti_info[i:i+s] for i in range(0, len(cifti_info), s)]

    logging.debug("Chunks: {}, Seeds: {}, Particles: {}".format(
        len(seed_chunks), len(seed_chunks[0]), particles))

    pool = multiprocessing.Pool(nbr_process)
    tracking_ = partial(tracking,
                        shm_file, mask_file, outdir, force_overwrite,
                        particles, step_size, max_lenght, max_angle, algorithm)

    # Numerated chunks and info: [(0, (sds0, cks0)), (1, (sds1, cks1))..]
    numerated_seeds_chunks = list(enumerate(zip(seed_chunks, info_chunks)))

    logging.debug("Starting multiprocessing environment")
    pool.map(tracking_, numerated_seeds_chunks)
    pool.close()
    pool.join()

    logging.debug("vmgenerator finished")
    finished_file = os.path.join(outdir, 'vmgenerator.end')
    open(finished_file, 'w').close()
