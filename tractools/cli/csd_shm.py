''' Given a dwi, fits an spherical harmonic model '''
from dipy.io import read_bvals_bvecs
from dipy.core.gradients import gradient_table
from dipy.reconst.csdeconv import auto_response
from dipy.reconst.csdeconv import ConstrainedSphericalDeconvModel

import nibabel
import numpy
import logging
import citrix

def csd_shm(dwi_file, bvals_file, bvecs_file, outfile,
            roi_center=None, roi_radius=10, fa_threshold=0.75, verbose=0):
    if verbose:
        logging.basicConfig(level=logging.DEBUG)

    if roi_center is not None and len(roi_center) != 3:
        raise ValueError("roi_center should be a 3-D position")

    bvals, bvecs = read_bvals_bvecs(bvals_file, bvecs_file)
    gtab = gradient_table(bvals, bvecs, b0_threshold=bvals.min())

    diffusion_img = nibabel.load(dwi_file)
    diffusion_data = diffusion_img.get_data()

    logging.debug("Fitting CSD model")
    response, ratio = auto_response(gtab, diffusion_data,
                                    roi_center, roi_radius,
                                    fa_threshold)
    logging.debug("Response: {response}, ratio: {ratio}")

    csd_model = ConstrainedSphericalDeconvModel(gtab, response, sh_order=6)
    csd_fit = csd_model.fit(diffusion_data)

    logging.debug("Saving SHM Coefficients")
    citrix.save(outfile, csd_fit.shm_coeff,
                affine=diffusion_img.affine, header=diffusion_img.header,
                version=1)
