#!/usr/bin/env python
import logging

import citrix
import streamlines as st
import numpy as np


def compute_visit_map(reference_shape, streamlines):

    visit_map = np.zeros(reference_shape)

    for streamline in streamlines:
        voxels = np.floor(streamline).astype(int)
        visit_map[tuple(np.transpose(voxels))] += 1

    return visit_map

def visit_map(reference_volume_file, streamlines_file, outfile,
              log_transform=False, normalize=False, binary=False):
    logging.basicConfig(level=logging.DEBUG)
    reference = citrix.load(reference_volume_file)

    streamlines = st.io.load(streamlines_file)
    streamlines.transform(np.linalg.inv(reference.affine))

    visit_map = compute_visit_map(reference.shape, streamlines)

    if log_transform:
        logging.debug("Applying log transform")

        visit_map += 1
        visit_map = np.log(visit_map)/np.log(visit_map.max())

    elif normalize:
        logging.debug("Normalizing")
        visit_map /= visit_map.max()

    elif binary:
        visit_map = visit_map > 0

    citrix.save(outfile, visit_map, reference.header, reference.affine, 1)
