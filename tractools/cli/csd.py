''' Given a dwi, fits an spherical harmonic model '''
import logging

from dipy.io import read_bvals_bvecs
from dipy.data import get_sphere
from dipy.direction import peaks_from_model
from dipy.core.gradients import gradient_table
from dipy.reconst.csdeconv import auto_response
from dipy.reconst.csdeconv import ConstrainedSphericalDeconvModel

import nibabel
import citrix

def csd(dwi_file, bvals_file, bvecs_file,
        outfile_shm=None, roi_center=None, roi_radius=10, fa_threshold=0.75,
        outfile_peaks_dir=None, outfile_peaks_val=None, mibrain_file=None,
        peak_mask=None, npeaks=5, min_separation=30, normalize=0,
        verbose=0):

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
    logging.debug("Response: {}, ratio: {}".format(response, ratio))

    csd_model = ConstrainedSphericalDeconvModel(gtab, response, sh_order=6)

    if outfile_shm is not None:
        csd_fit = csd_model.fit(diffusion_data)

        logging.debug("Saving SHM Coefficients")
        citrix.save(outfile_shm, csd_fit.shm_coeff,
                    affine=diffusion_img.affine, header=diffusion_img.header,
                    version=1)

    if outfile_peaks_dir is None and outfile_peaks_val is None and mibrain_file is None:
        return

    logging.debug("Computing peaks")
    sphere = get_sphere()

    if peak_mask is not None:
        peak_mask = nibabel.load(peak_mask).get_data().astype(bool)

    peaks = peaks_from_model(csd_model, diffusion_data, sphere, 0.5,
                             min_separation, peak_mask,
                             normalize_peaks=normalize, npeaks=npeaks,
                             sh_order=6)

    logging.debug("Saving peaks")
    if outfile_peaks_dir is not None:
        citrix.save(outfile_peaks_dir, peaks.peak_dirs,
                    affine=diffusion_img.affine, header=diffusion_img.header,
                    version=1)

    if outfile_peaks_val is not None:
        citrix.save(outfile_peaks_val, peaks.peak_values,
                    affine=diffusion_img.affine, header=diffusion_img.header,
                    version=1)

    if mibrain_file is not None:
        peaks = peaks.peak_dirs * peaks.peak_values[...,None]
        peaks = peaks.reshape(*peaks.shape[:-2], -1)
        citrix.save(mibrain_file, peaks,
                    affine=diffusion_img.affine, header=diffusion_img.header,
                    version=1)
