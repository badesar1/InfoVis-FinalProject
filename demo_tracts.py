from dipy.core.gradients import gradient_table
from dipy.data import get_fnames
from dipy.io.gradients import read_bvals_bvecs
from dipy.io.image import load_nifti, load_nifti_data
from dipy.reconst.csdeconv import auto_response_ssst
from dipy.reconst.shm import CsaOdfModel
from dipy.data import default_sphere
from dipy.direction import peaks_from_model
from dipy.viz import window, actor
from dipy.tracking.stopping_criterion import ThresholdStoppingCriterion
from dipy.tracking import utils
from dipy.tracking.local_tracking import LocalTracking
from dipy.tracking.streamline import Streamlines
from dipy.viz import colormap

import numpy as np
import mpdenoise as mp
import matplotlib.pyplot as plt

# load example data
hardi_fname, hardi_bval_fname, hardi_bvec_fname = get_fnames('stanford_hardi')
label_fname = get_fnames('stanford_labels')

data, affine, img = load_nifti(hardi_fname, return_img=True)
labels = load_nifti_data(label_fname)
bvals, bvecs = read_bvals_bvecs(hardi_bval_fname, hardi_bvec_fname)
gtab = gradient_table(bvals, bvecs)

print(img.shape)
# run the denoising
imgdata = img.get_fdata()
imgslice = imgdata[:,:,30:35,:]
Signal, Sigma, Npars, _, _ = mp.denoise(imgslice, kernel=[5,5,5],patchtype='nonlocal',patchsize=100,shrinkage='frob',algorithm='cordero-grande')
f, ax = plt.subplots(1, 3, sharey=True)
ax[0].imshow(imgslice[:,:,3,20])
ax[0].set_title('original')
ax[1].imshow(Signal[:,:,3,20])
ax[1].set_title('denoised')
ax[2].imshow(Sigma[:,:,3])
ax[2].set_title('sigma')
plt.show()

white_matter = (labels == 1) | (labels == 2)
response, ratio = auto_response_ssst(gtab, Signal, roi_radii=10, fa_thr=0.7)
csa_model = CsaOdfModel(gtab, sh_order=6)
csa_peaks = peaks_from_model(csa_model, Signal, default_sphere,
                             relative_peak_threshold=.8,
                             min_separation_angle=45,
                             mask=white_matter[:,:,30:35])

interactive = True
scene = window.Scene()
scene.add(actor.peak_slicer(csa_peaks.peak_dirs,
                            csa_peaks.peak_values,
                            colors=None))
window.record(scene, out_path='csa_direction_field.png', size=(900, 900))
if interactive:
    window.show(scene, size=(800, 800))

stopping_criterion = ThresholdStoppingCriterion(csa_peaks.gfa, .25)
sli = csa_peaks.gfa.shape[2] // 2
plt.figure('GFA')
plt.subplot(1, 2, 1).set_axis_off()
plt.imshow(csa_peaks.gfa[:, :, sli].T, cmap='gray', origin='lower')

plt.subplot(1, 2, 2).set_axis_off()
plt.imshow((csa_peaks.gfa[:, :, sli] > 0.25).T, cmap='gray', origin='lower')
plt.show()

seed_mask = (labels[:,:,30:35])
seeds = utils.seeds_from_mask(seed_mask, affine, density=[2, 2, 2])
streamlines_generator = LocalTracking(csa_peaks, stopping_criterion, seeds,
                                      affine=affine, step_size=.5)
streamlines = Streamlines(streamlines_generator)
# Prepare the display objects.
color = colormap.line_colors(streamlines)

streamlines_actor = actor.line(streamlines,
                                colormap.line_colors(streamlines))

# Create the 3D display.
scene = window.Scene()
scene.add(streamlines_actor)

# Save still images for this static example. Or for interactivity use
window.record(scene, out_path='tractogram_EuDX.png', size=(800, 800))
if interactive:
    window.show(scene)