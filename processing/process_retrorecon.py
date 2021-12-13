import os
import numpy as np
import nibabel as nib
import mpdenoise as mp
import scipy.io as sio

root = '/mnt/labspace/Projects/HR/highres_1020/nifti_retrorecon'
#twixroot = '/mnt/labspace/Projects/HR/highres_1020/twix/1mmiso/1mmiso_recon'
output = '/mnt/labspace/Projects/HR/highres_1020/recon'
nii1_mag = nib.load(os.path.join(root,'dicom_retrorecon_DIFF_MESO_b2_AP_1iso_te93_mag_20211020203424_24.nii'))
img1_mag = np.array(nii1_mag.get_data())
nii1_phi = nib.load(os.path.join(root,'dicom_retrorecon_DIFF_MESO_b2_AP_1iso_te93_phase_20211020203424_25.nii'))
img1_phi = np.array(nii1_phi.get_data())
nii2_mag = nib.load(os.path.join(root,'dicom_retrorecon_DIFF_MESO_b1_AP_1iso_te93_mag_20211020203424_12.nii'))
img2_mag = np.array(nii2_mag.get_data())
nii2_phi = nib.load(os.path.join(root,'dicom_retrorecon_DIFF_MESO_b1_AP_1iso_te93_phase_20211020203424_13.nii'))
img2_phi = np.array(nii2_phi.get_data())

print(np.min(img1_phi))
print(np.max(img1_phi))

mag = np.concatenate((img1_mag, img2_mag), axis=3)
phi = np.concatenate((img1_phi, img2_phi), axis=3)/4096*np.pi
img = mag*np.exp(1j*phi)

kernel_phase = [15,15,1]
kernel_mp = [5,5,5]

nii1_bval = np.loadtxt(os.path.join(root, 'dicom_retrorecon_DIFF_MESO_b2_AP_1iso_te93_mag_20211020203424_24.bval'))
nii1_bvec = np.loadtxt(os.path.join(root, 'dicom_retrorecon_DIFF_MESO_b2_AP_1iso_te93_mag_20211020203424_24.bvec'))
nii2_bval = np.loadtxt(os.path.join(root, 'dicom_retrorecon_DIFF_MESO_b1_AP_1iso_te93_mag_20211020203424_12.bval'))
nii2_bvec = np.loadtxt(os.path.join(root, 'dicom_retrorecon_DIFF_MESO_b1_AP_1iso_te93_mag_20211020203424_12.bvec'))
bvec = np.concatenate((nii1_bvec,nii2_bvec),axis=1)
bval = np.concatenate((nii1_bval,nii2_bval),axis=0)
np.savetxt(os.path.join(output,'2mmiso.bvec'),bvec,delimiter=' ')
np.savetxt(os.path.join(output,'2mmiso.bval'),bval,delimiter=' ')


# noMPmagnitude:
# nii_noMP = nib.Nifti1Image(mag, nii1_mag.affine, nii1_mag.header)
# nib.save(nii_noMP, os.path.join(output,'2mmiso_noMPmagnitude.nii'))
#imdict = sio.loadmat(os.path.join(twixroot,'noMPmagnitude.mat'))
#img_proc = imdict["noMPmagnitude"]
#nii_noMP = nib.Nifti1Image(img_proc, nii1_mag.affine, nii1_mag.header)
#nib.save(nii_noMP, os.path.join(output,'1mmiso_noMPmagnitude_RMRrecon.nii'))

# print(np.min(mag))
# print(np.max(mag))
# print(mag.shape)

# MPmagnitude
# Signal, Sigma, Npars, _, _ = mp.denoise(mag, kernel=kernel_mp,patchtype='nonlocal',patchsize=100,shrinkage='threshold',algorithm='cordero-grande')
# nii_MP = nib.Nifti1Image(Signal, nii1_mag.affine, nii1_mag.header)
# nib.save(nii_MP, os.path.join(output,'2mmiso_MPmagnitude_noshrink.nii'))
# nii_MP = nib.Nifti1Image(Sigma, nii1_mag.affine, nii1_mag.header)
# nib.save(nii_MP, os.path.join(output,'2mmiso_MPmagnitude_sigma_noshrink.nii'))
# nii_MP = nib.Nifti1Image(Npars, nii1_mag.affine, nii1_mag.header)
# nib.save(nii_MP, os.path.join(output,'2mmiso_MPmagnitude_npars_noshrink.nii'))
#imdict = sio.loadmat(os.path.join(twixroot,'MPmagnitude.mat'))
#img_proc = imdict["MPmagnitude"]
#nii_MP = nib.Nifti1Image(img_proc, nii1_mag.affine, nii1_mag.header)
#nib.save(nii_MP, os.path.join(output,'1mmiso_MPmagnitude_RMRrecon.nii'))

# # MPcomplex
Signal, Sigma, Npars, _, _ = mp.denoise(img, kernel=kernel_phase,patchtype='box',shrinkage='threshold',algorithm='cordero-grande',crop=75)
phi = np.angle(Signal)
print(Signal.dtype)

nii_MP = nib.Nifti1Image(phi, nii1_mag.affine, nii1_mag.header)
nib.save(nii_MP, os.path.join(output,'1mmiso_MPcomplex_phaseEst.nii'))
img_np = img*np.exp(-1j*phi)
Signal, Sigma, Npars, _, _ = mp.denoise(np.real(img_np), kernel=kernel_mp,patchtype='nonlocal',patchsize=100,shrinkage='frob',algorithm='cordero-grande')
Signal = np.absolute(Signal).astype('uint16')
nii_MP = nib.Nifti1Image(Signal, nii1_mag.affine, nii1_mag.header)
nib.save(nii_MP, os.path.join(output,'1mmiso_MPcomplex_crop.nii'))
nii_MP = nib.Nifti1Image(Sigma, nii1_mag.affine, nii1_mag.header)
nib.save(nii_MP, os.path.join(output,'1mmiso_MPcomplex_sigma.nii'))
nii_MP = nib.Nifti1Image(Npars, nii1_mag.affine, nii1_mag.header)
nib.save(nii_MP, os.path.join(output,'1mmiso_MPcomplex_npars.nii'))
#imdict = sio.loadmat(os.path.join(twixroot,'MPcomplex.mat'))
#img_proc = imdict["MPcomplex"]
#img_proc = np.absolute(img_proc*1000).astype('uint16')
#nii_MPc = nib.Nifti1Image(img_proc, nii1_mag.affine, nii1_mag.header)
#nib.save(nii_MPc, os.path.join(output,'1mmiso_MPcomplex_RMRrecon.nii'))
