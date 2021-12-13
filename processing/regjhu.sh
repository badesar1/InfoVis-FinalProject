#!/bin/sh

#  regJHU.sh

export FSLOUTPUTTYPE=NIFTI

diffusionDir=/mnt/labspace/Projects/HR/highres_1020/recon
subjs=( designer_MPcomplex designer_MPcomplex_2mmiso )

template=${FSLDIR}/data/atlases/JHU/JHU-ICBM-FA-1mm.nii.gz
label=${FSLDIR}/data/atlases/JHU/JHU-ICBM-labels-1mm.nii.gz

for i in ${subjs[@]}; do
    i="${i//[$'\t\r\n ']}"
    echo $i

    outdir=${diffusionDir}/$i
    paramdir=${diffusionDir}/$i

    # create output directory
    # mkdir $outdir
    # in case there are NaN values in fa.nii, remove them
    mrcalc -force $paramdir/fa.nii -finite $paramdir/fa.nii 0 -if $paramdir/fa_noNaN.nii
    # compute the affine registration to JHU template
    flirt -in $paramdir/fa_noNaN.nii -ref $template -omat $outdir/fa2jhu_affine.txt -dof 12 -cost corratio
    # compute the nonlinear warp to JHU template, initialize with affine
    fnirt --in=$paramdir/fa_noNaN.nii --ref=${template} --aff=$outdir/fa2jhu_affine.txt --cout=$outdir/fa2jhu_nonlin --iout=$outdir/fa2jhu_nonlin_fa --config=FA_2_FMRIB58_1mm
    # invert the nonlinear warp
    invwarp -w $outdir/fa2jhu_nonlin -o $outdir/jhu2fa_nonlin -r $paramdir/fa_noNaN.nii
    # apply the JHU-to-subject warp on ROIs
    applywarp --in=${label} --ref=$paramdir/fa_noNaN.nii --warp=$outdir/jhu2fa_nonlin --out=$outdir/jhulabels --interp=nn
done
