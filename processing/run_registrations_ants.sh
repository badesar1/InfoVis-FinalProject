#export ANTSPATH=/data/apps/AdvancedNormTools/install/bin

subjects=(SMR009_19JAN21
SMR010_20JAN21
SMR011_21JAN21
SMR012_24JAN21
SMR013_25JAN21
SMR014_26JAN21
SMR015_27JAN21
SMR016_01FEB21
SMR017_03FEB21
SMR018_03FEB21
SMR019_09FEB21
SMR020_10FEB21
SMR021_10FEB21
SMR022_14FEB21
SMR023_14FEB21
SMR024_17FEB21
SMR025_18FEB21
SMR026_18FEB21
SMR027_19FEB21
SMR028_22FEB21)

scans=(PRISMA_1 PRISMA_2 SKYRA_1 SKYRA_2 SKYRA_3)
params=(md mk fa mw ad rd aw rw)
#methods=(MPcomplex MPmagnitude noMPmagnitude)
methods=(MPcomplex MPmagnitude noMPmagnitude)
root=/mnt/labspace/Projects/SM_reproducibility/SMfixedTE

## subject specific population templates
# rm -r $root/subject_population_templates_nov21
stroot=$root/subject_population_templates_nov21_2
mkdir $stroot
for subj in ${subjects[@]}; do
    #mkdir $root/subject_population_templates/${subj}_template
    mkdir $stroot/${subj}_b0
    mkdir $stroot/${subj}_fa
    mkdir $stroot/${subj}_masks
    for scan in ${scans[@]}; do
        mrcalc -force $root/$subj/YARRApreproc_sept21/$scan/MPcomplex_nov21/wlls_smooth/b0_dti.nii -finite $root/$subj/YARRApreproc_sept21/$scan/MPcomplex_nov21/wlls_smooth/b0_dti.nii 0 -if -abs - | mrtransform -identity - $stroot/${subj}_b0/${scan}_b0.nii
        mrcalc -force $root/$subj/YARRApreproc_sept21/$scan/MPcomplex_nov21/wlls_smooth/fa_dti.nii -finite $root/$subj/YARRApreproc_sept21/$scan/MPcomplex_nov21/wlls_smooth/fa_dti.nii 0 -if -abs - | mrtransform -identity - $stroot/${subj}_fa/${scan}_fa.nii
        #cp $root/$subj/YARRApreproc_sept21/$scan/MPcomplex/processing/post_eddy_mask_mask.nii.gz $root/subject_population_templates_nov21/${subj}_masks/${scan}_mask.nii.gz 
        mrcalc -force $root/$subj/YARRApreproc_sept21/$scan/MPcomplex_nov21/wlls_smooth/fa_dti.nii -finite 1 0 -if - | mrtransform -identity - $stroot/${subj}_masks/${scan}_mask.nii
    done
done
for subj in ${subjects[@]}; do
    #population_template -voxel_size 2 -force $root/subject_population_templates_nov21/${subj}_b0 $root/subject_population_templates_nov21/${subj}_b0/template_b0.nii $root/subject_population_templates_nov21/${subj}_fa $root/subject_population_templates_nov21/${subj}_fa/template_fa.nii -type rigid -linear_transformations_dir $root/subject_population_templates_nov21/${subj}_transforms -mask_dir $root/subject_population_templates_nov21/${subj}_masks -template_mask $root/subject_population_templates_nov21/${subj}_masks/template_mask.nii -transformed_dir $root/subject_population_templates_nov21/${subj}_b0/transformed,$root/subject_population_templates_nov21/${subj}_fa/transformed
    population_template -voxel_size 2 -force $stroot/${subj}_fa $stroot/${subj}_fa/template_fa.nii -type rigid -linear_transformations_dir $stroot/${subj}_transforms -mask_dir $stroot/${subj}_masks -template_mask $stroot/${subj}_masks/template_mask.nii -transformed_dir $stroot/${subj}_fa/transformed
done

# ## compute subject wise CVs
for subj in ${subjects[@]}; do
    mkdir $stroot/${subj}_CV_wlls_smooth
    for p in ${params[@]}; do
        for scan in ${scans[@]}; do
            mrtransform -force -identity $root/$subj/YARRApreproc_sept21/$scan/MPcomplex_nov21/wlls_smooth/${p}.nii $root/$subj/YARRApreproc_sept21/$scan/MPcomplex_nov21/wlls_smooth/${p}_i.nii
            mrtransform -force -identity $root/$subj/YARRApreproc_sept21/$scan/MPmagnitude/wlls_smooth/${p}.nii $root/$subj/YARRApreproc_sept21/$scan/MPmagnitude/wlls_smooth/${p}_i.nii
            mrtransform -force -identity $root/$subj/YARRApreproc_sept21/$scan/noMPmagnitude/wlls_smooth/${p}.nii $root/$subj/YARRApreproc_sept21/$scan/noMPmagnitude/wlls_smooth/${p}_i.nii
            mrtransform -force -linear $stroot/${subj}_transforms/${scan}_fa.txt -template $stroot/${subj}_fa/template_fa.nii $root/$subj/YARRApreproc_sept21/$scan/MPcomplex_nov21/wlls_smooth/${p}_i.nii $stroot/${subj}_CV_wlls_smooth/${p}_${scan}_MPcomplex.nii
            mrtransform -force -linear $stroot/${subj}_transforms/${scan}_fa.txt -template $stroot/${subj}_fa/template_fa.nii $root/$subj/YARRApreproc_sept21/$scan/MPmagnitude/wlls_smooth/${p}_i.nii $stroot/${subj}_CV_wlls_smooth/${p}_${scan}_MPmagnitude.nii
            mrtransform -force -linear $stroot/${subj}_transforms/${scan}_fa.txt -template $stroot/${subj}_fa/template_fa.nii $root/$subj/YARRApreproc_sept21/$scan/noMPmagnitude/wlls_smooth/${p}_i.nii $stroot/${subj}_CV_wlls_smooth/${p}_${scan}_noMPmagnitude.nii
        done
        cd $root/subject_population_templates_nov21/${subj}_CV_wlls_smooth
        for meth in ${methods[@]}; do
            P1=${p}_PRISMA_1_${meth}.nii
            P2=${p}_PRISMA_2_${meth}.nii
            S1=${p}_SKYRA_1_${meth}.nii
            S2=${p}_SKYRA_2_${meth}.nii
            S3=${p}_SKYRA_3_${meth}.nii
            mrcalc -force $P1 $P2 -subtract -abs $P1 $P2 -add 2 -div -div $stroot/${subj}_CV_wlls_smooth/${p}_P1P2_${meth}.nii
            mrcalc -force $S1 $S2 -subtract -abs $S1 $S2 -add 2 -div -div $stroot/${subj}_CV_wlls_smooth/${p}_S1S2_${meth}.nii
            mrcalc -force $S1 $S3 -subtract -abs $S1 $S3 -add 2 -div -div $stroot/${subj}_CV_wlls_wmooth/${p}_S1S3_${meth}.nii
            mrcalc -force $P1 $S3 -subtract -abs $P1 $S3 -add 2 -div -div $stroot/${subj}_CV_wlls_smooth/${p}_P1S3_${meth}.nii
        done
    done
done

# # ## overall subject template
## rm -r $root/overall_population_template_nov21
otroot=$root/overall_population_template_nov21_2
mkdir $otroot
mkdir $otroot/template_b0_HR
mkdir $otroot/template_fa_HR
mkdir $otroot/template_masks_HR
mkdir $otroot/warps_HR
for subj in ${subjects[@]}; do
    #cp $root/subject_population_templates_nov21/${subj}_b0/template_b0.nii $root/overall_population_template_nov21/template_b0_HR/${subj}_template_b0.nii
    cp $stroot/${subj}_fa/template_fa.nii $otroot/template_fa_HR/${subj}_template_fa.nii
    cp $stroot/${subj}_masks/template_mask.nii $otroot/template_masks_HR/${subj}_template_mask.nii
done
#population_template -force -voxel_size 0.5 $root/overall_population_template_nov21/template_b0_HR $root/overall_population_template_nov21/template_b0_HR/population_template_b0_HR.nii $root/overall_population_template_nov21/template_fa_HR $root/overall_population_template_nov21/template_fa_HR/population_template_fa_HR.nii -warp_dir $root/overall_population_template_nov21/warps_HR -mask_dir $root/overall_population_template_nov21/template_masks_HR -template_mask $root/overall_population_template_nov21/template_masks_HR/population_template_mask_HR.nii -transformed_dir $root/overall_population_template_nov21/template_b0_HR/transformed,$root/overall_population_template_nov21/template_fa_HR/transformed
population_template -force -scratch $otroot/scratch -voxel_size 0.5 $otroot/template_fa_HR_orig $otroot/template_fa_HR/population_template_fa_HR.nii -warp_dir $otroot/warps_HR -mask_dir $otroot/template_masks_HR -template_mask $otroot/template_masks_HR/population_template_mask_HR.nii -transformed_dir $otroot/template_fa_HR/transformed

# # # ## concatenate subject to population warp
mkdir $root/warped_maps_nov21_wlls_smooth
for subj in ${subjects[@]}; do
    warpconvert -template $otroot/template_fa_HR/population_template_fa_HR.nii $otroot/warps_HR/${subj}_template_fa.mif warpfull2deformation $otroot/warps_HR/${subj}_template_deformation.mif
    for scan in ${scans[@]}; do
        transformcompose $stroot/${subj}_transforms/${scan}_fa.txt $otroot/warps_HR/${subj}_template_deformation.mif $stroot/${subj}_transforms/${scan}_warp_HR.mif
        for p in ${params[@]}; do
            mrtransform -force -warp $stroot/${subj}_transforms/${scan}_warp_HR.mif $root/$subj/YARRApreproc_sept21/$scan/MPcomplex_nov21/wlls_smooth/${p}_i.nii $root/warped_maps_nov21_wlls_smooth/${subj}_${scan}_MPcomplex_${p}_HR.nii &
            mrtransform -force -warp $stroot/${subj}_transforms/${scan}_warp_HR.mif $root/$subj/YARRApreproc_sept21/$scan/MPmagnitude/wlls_smooth/${p}_i.nii $root/warped_maps_nov21_wlls_smooth/${subj}_${scan}_MPmagnitude_${p}_HR.nii &
            mrtransform -force -warp $stroot/${subj}_transforms/${scan}_warp_HR.mif $root/$subj/YARRApreproc_sept21/$scan/noMPmagnitude/wlls_smooth/${p}_i.nii $root/warped_maps_nov21_wlls_smooth/${subj}_${scan}_noMPmagnitude_${p}_HR.nii 
        done
    done
done

# # ## bring subject wise CVs to population space
mkdir $root/overall_population_template_nov21/template_CV_wlls_smooth
for p in ${params[@]}; do
    for meth in ${methods[@]}; do
        for subj in ${subjects[@]}; do
            mrtransform -force -template $otroot/template_fa_HR/population_template_fa_HR.nii -warp $otroot/warps_HR/${subj}_template_deformation.mif $stroot/${subj}_CV_wlls_smooth/${p}_P1P2_${meth}.nii $otroot/template_CV_wlls_smooth/${subj}_${p}_${meth}_P1P2.nii &
            mrtransform -force -template $otroot/template_fa_HR/population_template_fa_HR.nii -warp $otroot/warps_HR/${subj}_template_deformation.mif $stroot/${subj}_CV_wlls_smooth/${p}_S1S2_${meth}.nii $otroot/template_CV_wlls_smooth/${subj}_${p}_${meth}_S1S2.nii &
            mrtransform -force -template $otroot/template_fa_HR/population_template_fa_HR.nii -warp $otroot/warps_HR/${subj}_template_deformation.mif $stroot/${subj}_CV_wlls_smooth/${p}_S1S3_${meth}.nii $otroot/template_CV_wlls_smooth/${subj}_${p}_${meth}_S1S3.nii &
            mrtransform -force -template $otroot/template_fa_HR/population_template_fa_HR.nii -warp $otroot/warps_HR/${subj}_template_deformation.mif $stroot/${subj}_CV_wlls_smooth/${p}_P1S3_${meth}.nii $otroot/template_CV_wlls_smooth/${subj}_${p}_${meth}_P1S3.nii
        done
        mrmath -force $otroot/template_CV_wlls_smooth/*_${p}_${meth}_P1P2.nii mean $otroot/template_CV_wlls_smooth/mean_${p}_${meth}_P1P2.nii &
        mrmath -force $otroot/template_CV_wlls_smooth/*_${p}_${meth}_S1S2.nii mean $otroot/template_CV_wlls_smooth/mean_${p}_${meth}_S1S2.nii &
        mrmath -force $otroot/template_CV_wlls_smooth/*_${p}_${meth}_S1S3.nii mean $otroot/template_CV_wlls_smooth/mean_${p}_${meth}_S1S3.nii &
        mrmath -force $otroot/template_CV_wlls_smooth/*_${p}_${meth}_P1S3.nii mean $otroot/template_CV_wlls_smooth/mean_${p}_${meth}_P1S3.nii
    done
done



# # ## population to JHU warp
# flirt -in $root/overall_population_template_nov21/template_fa_HR/population_template_fa_HR.nii -ref /data/apps/fsl/data/atlases/JHU/JHU-ICBM-FA-1mm.nii.gz -dof 12 -omat $root/overall_population_template_nov21/template_fa_HR/population_template_fa2jhu_flirt.mat
# fnirt --in=$root/overall_population_template_nov21/template_fa_HR/population_template_fa_HR.nii --ref=/data/apps/fsl/data/atlases/JHU/JHU-ICBM-FA-1mm.nii.gz --aff=$root/overall_population_template_nov21/template_fa_HR/population_template_fa2jhu_flirt.mat --config=FA_2_FMRIB58_1mm --cout=$root/overall_population_template_nov21/template_fa_HR/population_template_fa2jhu_fnirt
# invwarp -w $root/overall_population_template_nov21/template_fa_HR/population_template_fa2jhu_fnirt -o $root/overall_population_template_nov21/template_fa_HR/population_template_jhu2fa_fnirt -r $root/overall_population_template_nov21/template_fa_HR/population_template_fa_HR.nii

# warpinit -force /data/apps/fsl/data/atlases/JHU/JHU-ICBM-FA-1mm.nii.gz $root/overall_population_template_nov21/template_fa_HR/population_template_jhu2fa_identity_warp.nii
# applywarp --in=$root/overall_population_template_nov21/template_fa_HR/population_template_jhu2fa_identity_warp.nii --ref=$root/overall_population_template_nov21/template_fa_HR/population_template_fa_HR.nii --warp=$root/overall_population_template_nov21/template_fa_HR/population_template_jhu2fa_fnirt --out=$root/overall_population_template_nov21/template_fa_HR/population_template_jhu2fa_fnirt2mrtrix_warp
#mrtransform -force -interp nearest -warp $root/overall_population_template_nov21/template_fa_HR/population_template_fa2jhu_fnirt2mrtrix_warp.mif.nii.gz $root/overall_population_template_nov21/template_fa_HR/population_template_fa_HR.nii $root/overall_population_template_nov21/template_fa_HR/population_template_2JHU.nii
#mrtransform -force -interp nearest -warp $root/overall_population_template_nov21/template_fa_HR/population_template_jhu2fa_fnirt2mrtrix_warp.nii.gz /mnt/labspace/Projects/ADDF/PROCESSED/WM_ROI_analysis/JHU-ICBM-labels/JHU_labels_indthr/JHU-ICBM-labels-indthr.nii $root/overall_population_template_nov21/template_fa_HR/population_template_JHU_labels.nii

# mkdir $root/subject_rois_nov21
# for subj in ${subjects[@]}; do
#     for scan in ${scans[@]}; do
#         warpinvert $root/subject_population_templates_nov21/${subj}_transforms/${scan}_warp_HR.mif $root/subject_population_templates_nov21/${subj}_transforms/${scan}_invwarp_HR.mif
#         transformcompose -force $root/overall_population_template_nov21/template_fa_HR/population_template_jhu2fa_fnirt2mrtrix_warp.nii.gz $root/subject_population_templates_nov21/${subj}_transforms/${scan}_invwarp_HR.mif $root/subject_population_templates_nov21/${subj}_transforms/${scan}_warp_jhu2subj.mif
#         mrtransform -force -template $root/subject_population_templates_nov21/${subj}_fa/${scan}_fa.nii -interp nearest -warp $root/subject_population_templates_nov21/${subj}_transforms/${scan}_warp_jhu2subj.mif /mnt/labspace/Projects/ADDF/PROCESSED/WM_ROI_analysis/JHU-ICBM-labels/JHU_labels_indthr/JHU-ICBM-labels-indthr.nii $root/subject_rois_nov21/${subj}_${scan}_rois_HR.nii
#         mrcalc $root/subject_population_templates_nov21/${subj}_masks/${scan}_mask.nii $root/subject_rois_nov21/${subj}_${scan}_rois_HR.nii -mult $root/subject_rois_nov21/${subj}_${scan}_rois_thr.nii
#     done
# done




