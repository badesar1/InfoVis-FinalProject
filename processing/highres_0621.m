clear

basefolder='/cbi05data/data1/Hamster/Ben/highres_0621/HR_MESO_22JUNE21_twix/1p25mmiso/mat_files_all';
addpath(genpath('/data/Ben/scripts'))
addpath(genpath('/cbi05data/data1/Hamster/Ben/coilstuff'))
cd(basefolder);

s = 1;
for SL=1:100
    load(['imgreconX_RMRFULL',num2str(s),'_D1_TRG.mat'])
    I(:,:,s,:)=imgrecon(:,:,1,:,1);
    Idn(:,:,s,:) = imgrecon(:,:,1,:,2);
s=s+1;
end
%I=single(I);
%save('/cbi05data/data1/Hamster/Ben/highres_0621/I.mat','I');
niftiwrite(flip(flip(abs(I),1),2),fullfile('/cbi05data/data1/Hamster/Ben/highres_0621/HR_MESO_22JUNE21_twix/1p25mmiso/1p25mmiso_recon_trg.nii'))
niftiwrite(flip(flip(abs(Idn),1),2),fullfile('/cbi05data/data1/Hamster/Ben/highres_0621/HR_MESO_22JUNE21_twix/1p25mmiso/1p25mmiso_recon_RMT_trg.nii'))
system('fslcpgeom /cbi05data/data1/Hamster/Ben/highres_0621/HR_MESO_22JUNE21_nifti/1p2mmiso_processing/tmp_dwi.nii /cbi05data/data1/Hamster/Ben/highres_0621/HR_MESO_22JUNE21_twix/1p25mmiso/1p25mmiso_recon_trg.nii');
system('fslcpgeom /cbi05data/data1/Hamster/Ben/highres_0621/HR_MESO_22JUNE21_nifti/1p2mmiso_processing/tmp_dwi.nii /cbi05data/data1/Hamster/Ben/highres_0621/HR_MESO_22JUNE21_twix/1p25mmiso/1p25mmiso_recon_RMT_trg.nii');

%%
addpath('/mnt/labspace/Projects/SM_reproducibility/scripts/reconstruction');
kernelPhase = [15,15,1];
[dn_phase, Ns_phase, Ps_phase]=MP_Loop4_crop(single(I), kernelPhase);
[Cx,Cy,Cz,Cd]=size(I);
imgrecon_P = 0.*I;
for nz=1:Cz
    for nd=1:Cd
        imgrecon_P(:,:,nz,nd) = single(I(:,:,nz,nd).*exp(-1i*angle(dn_phase(:,:,nz,nd))));
    end
end
%save('/cbi05data/data1/Hamster/Ben/highres_0621/dn_phase.mat','dn_phase');


%%
kernelMP = [5,5,5];
psize = 125;
bval = importdata('/cbi05data/data1/Hamster/Ben/highres_0621/HR_MESO_22JUNE21_nifti/1p2mmiso_parameters/dwi_designer.bval')'./1000;
bvec = importdata('/cbi05data/data1/Hamster/Ben/highres_0621/HR_MESO_22JUNE21_nifti/1p2mmiso_parameters/dwi_designer.bvec')';
bvec1 = bvec(bval==0|bval==1,:);
bval1 = bval(bval==0|bval==1);

rmpath('/mnt/labspace/Projects/SM_reproducibility/scripts/reconstruction');
addpath(genpath('/cbi05data/data1/Hamster/DESIGNER/utils'));

[Idn, S, P] = MPnonlocal(real(imgrecon_P), kernelMP,psize,'h');
niftiwrite(flip(flip(abs(Idn),1),2),fullfile('/cbi05data/data1/Hamster/Ben/highres_0621/HR_MESO_22JUNE21_twix/1p25mmiso/1mmiso_recon_MPcomplex_trg.nii'))
niftiwrite(flip(flip(abs(S),1),2),fullfile('/cbi05data/data1/Hamster/Ben/highres_0621/HR_MESO_22JUNE21_twix/1p25mmiso/1mmiso_sigma_MPcomplex_trg.nii'))
niftiwrite(flip(flip(abs(P),1),2),fullfile('/cbi05data/data1/Hamster/Ben/highres_0621/HR_MESO_22JUNE21_twix/1p25mmiso/1mmiso_npars_MPcomplex_trg.nii'))
system('fslcpgeom /cbi05data/data1/Hamster/Ben/highres_0621/HR_MESO_22JUNE21_nifti/1p2mmiso_processing/tmp_dwi.nii /cbi05data/data1/Hamster/Ben/highres_0621/HR_MESO_22JUNE21_twix/1p25mmiso/1mmiso_recon_MPcomplex_trg.nii');
system('fslcpgeom /cbi05data/data1/Hamster/Ben/highres_0621/HR_MESO_22JUNE21_nifti/1p2mmiso_processing/tmp_dwi.nii /cbi05data/data1/Hamster/Ben/highres_0621/HR_MESO_22JUNE21_twix/1p25mmiso/1mmiso_sigma_MPcomplex_trg.nii');
system('fslcpgeom /cbi05data/data1/Hamster/Ben/highres_0621/HR_MESO_22JUNE21_nifti/1p2mmiso_processing/tmp_dwi.nii /cbi05data/data1/Hamster/Ben/highres_0621/HR_MESO_22JUNE21_twix/1p25mmiso/1mmiso_npars_MPcomplex_trg.nii');


Idn1 = Idn(:,:,:,bval==0|bval==1);
% [b0, dt, md, rd, ad, fa, fe] = dti_fit(Idn1,[bvec1,bval1]);
% save('/cbi05data/data1/Hamster/Ben/highres_0621/Idn.mat','Idn');
% save('/cbi05data/data1/Hamster/Ben/highres_0621/fe.mat','fe');

%%

% bval1 = importdata('/cbi05data/data1/Hamster/Ben/highres_0621/HR_MESO_22JUNE21_nifti/HR_MESO_22JUNE21_dicom_DIFF_1mmiso_b1000_20dir_te85_20210622154709_21.bval');
% bvec1 = importdata('/cbi05data/data1/Hamster/Ben/highres_0621/HR_MESO_22JUNE21_nifti/HR_MESO_22JUNE21_dicom_DIFF_1mmiso_b1000_20dir_te85_20210622154709_21.bvec');
% 
% bval2 = importdata('/cbi05data/data1/Hamster/Ben/highres_0621/HR_MESO_22JUNE21_nifti/HR_MESO_22JUNE21_dicom_DIFF_1mmiso_b2000_64dir_te85_20210622154709_15.bval');
% bvec2 = importdata('/cbi05data/data1/Hamster/Ben/highres_0621/HR_MESO_22JUNE21_nifti/HR_MESO_22JUNE21_dicom_DIFF_1mmiso_b2000_64dir_te85_20210622154709_15.bvec');
% 
% bval = cat(2,bval2,bval1);
% dlmwrite('/cbi05data/data1/Hamster/Ben/highres_0621/HR_MESO_22JUNE21_twix/1mmiso/b.bval',bval,' ');
% bvec = cat(2,bvec2,bvec1);
% dlmwrite('/cbi05data/data1/Hamster/Ben/highres_0621/HR_MESO_22JUNE21_twix/1mmiso/b.bvec',bvec,' ');
