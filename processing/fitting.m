addpath(genpath('/cbi05data/data1/Hamster/DESIGNER/utils'))
cd /mnt/labspace/Projects/HR/HIGHRES_080321/

bvec1 = importdata('NIFTI_online_recon/DICOM_online_recon_DIFF_MESO2_1p2iso_te73_no_SMS_20210803135345_50.bvec');
bvec2 = importdata('NIFTI_online_recon/DICOM_online_recon_DIFF_MESO1_1p2iso_te73_no_SMS_20210803135345_44.bvec');
bval1 = importdata('NIFTI_online_recon/DICOM_online_recon_DIFF_MESO2_1p2iso_te73_no_SMS_20210803135345_50.bval');
bval2 = importdata('NIFTI_online_recon/DICOM_online_recon_DIFF_MESO1_1p2iso_te73_no_SMS_20210803135345_44.bval');
bval = cat(2,bval1,bval2)'/1000;
bvec = cat(2,bvec1,bvec2)';
grad = cat(2,bvec,bval);

data_mpcomplex = load('twix/1p2_no_SMS/1p2mmiso_te73_no_SMS/MPcomplex.mat');
dwi_mpc = data_mpcomplex.MPcomplex;
% data_rmt = load('twix/1p2_no_SMS/1p2mmiso_te73_no_SMS/RMT.mat');
% dwi_rmt = data_rmt.RMT;
data_mpmag = load('twix/1p2_no_SMS/1p2mmiso_te73_no_SMS/MPmagnitude.mat');
dwi_mag = data_mpmag.MPmagnitude;
data_nompmag = load('twix/1p2_no_SMS/1p2mmiso_te73_no_SMS/noMPmagnitude.mat');
dwi_no = data_nompmag.noMPmagnitude;

D = {dwi_mpc,dwi_mag,dwi_no};
FE=[];
MD = [];
MK = [];
for i = 1:3
%[b0, dt, md, rd, ad, fa, fe] = dti_fit(D{i}, grad);
[b0, dt] = dki_fit(D{i}, grad);
[fa,md,rd,ad,fe,mk,rk,ak,l1,l2,l3] = dki_parameters(dt);
FE = cat(2,FE,rot90(fe));
MD = cat(2,MD,rot90(md));
MK = cat(2,MK,rot90(mk));
end

figure;
imagesc(abs(squeeze(FE(:,:,55,:))));
figure; 
imagesc(abs(MD(:,:,55)))

niftiwrite(fe,'fe_1p2mmiso_mpcomplex.nii');

%%%%
%%
clear
addpath(genpath('/cbi05data/data1/Hamster/DESIGNER/utils'))
addpath('/cbi05data/data1/Hamster/Ben');
addpath('/data/apps/mrtrix3/matlab');

root = '/mnt/labspace/Projects/MESO_v2.0/fitting_test_reproducability_NOV21';
subjects = {'SMR011_21JAN21','SMR015_27JAN21'};
scans = {'PRISMA_1','PRISMA_2','SKYRA_1','SKYRA_2','SKYRA_3'};

for i = 1:numel(subjects)
for j = 1:numel(scans)

dpath = fullfile(root,subjects{i},scans{j});
data = niftiread(fullfile(dpath,'data.nii'));
%sigma = niftiread(fullfile(dpath,'sigma.nii'));
%mif = read_mrtrix(fullfile(dpath,'noise.mif'));
%sigma = flip(mif.data,1);

%data_rc = sqrt(data.^2 - sigma.^2);
outdir = fullfile('/mnt/labspace/Ben/fitting_nov_repro_akc_robust2',subjects{i},scans{j});
mkdir(outdir);

detectoutliers = 0;
cumulants = 0;
dti = 1;
dki = 1;
wmti = 0;
fitconstraints = '0,0,0';
akc = 1;
dkiroot = '/cbi05data/data1/Hamster/DESIGNER/utils';
fitwdki = 0;
smooth = 0;
robust = 1;

tensorfitting_sm(data,dpath,outdir,detectoutliers,cumulants,dti,dki,wmti,fitconstraints,akc,dkiroot,fitwdki,smooth,robust);
copyfile(fullfile(outdir,'dkt.nii'),fullfile(dpath,'DKIfits','dkt_ben_akc_robust2.nii'));
end
end
%%
%clear
addpath(genpath('/cbi05data/data1/Hamster/DESIGNER/utils'))
addpath('/cbi05data/data1/Hamster/Ben');

root = '/mnt/shepht01lab/Shepherd-Lab-Users/Ben/110621_clinical_case/recon';
dpath = fullfile(root,'designer_mpcomplex_proc2');
data = niftiread(fullfile(dpath,'dwi_designer.nii'));

outdir = fullfile(dpath,'smooth');
mkdir(outdir);

detectoutliers = 0;
cumulants = 0;
dti = 1;
dki = 1;
wmti = 0;
fitconstraints = '0,0,0';
akc = 1;
dkiroot = '/cbi05data/data1/Hamster/DESIGNER/utils';
fitwdki = 1;
smooth = 1;
tensorfitting(dpath,outdir,detectoutliers,cumulants,dti,dki,wmti,fitconstraints,akc,dkiroot,fitwdki,smooth);

%%
subj='SMR015_27JAN21';
scans = {'PRISMA_1','PRISMA_2','SKYRA_3','SKYRA_1','SKYRA_2'};
root = '/mnt/labspace/Projects/SM_reproducibility/SMfixedTE';
out = '/mnt/labspace/Projects/MESO_v2.0/fitting_test_reproducability_NOV21';
for i = 1:numel(scans)
    %mkdir(fullfile(out,subj,scans{i},'DKIfits'));
    b = importdata(fullfile(out,subj,scans{i},'data.bval'));
    delete(fullfile(out,subj,scans{i},'data.bval'));
    b = b./1000;
    dlmwrite(fullfile(out,subj,scans{i},'data.bval'), b, 'delimiter',' ');
%     copyfile(fullfile(root,subj,'YARRApreproc_sept21',scans{i},'MPcomplex','processing','dwiec_lte.nii.gz'),...
%         fullfile(out,subj,scans{i},'data.nii.gz'));
%     copyfile(fullfile(root,subj,'YARRApreproc_sept21',scans{i},'MPcomplex','processing','dwi_lte.bval'),...
%         fullfile(out,subj,scans{i},'data.bval'));
%     copyfile(fullfile(root,subj,'YARRApreproc_sept21',scans{i},'MPcomplex','processing','dwiec_lte.eddy_rotated_bvecs'),...
%         fullfile(out,subj,scans{i},'data.bvec'));
%     copyfile(fullfile(root,subj,'YARRApreproc_sept21',scans{i},'MPcomplex','processing','post_eddy_mask_mask.nii.gz'),...
%         fullfile(out,subj,scans{i},'mask.nii.gz'));
end




