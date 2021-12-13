clear
close all
clc
root = '/Volumes/Research/fieree01lab/labspace/Santiago/Preprocessed_volunteers';
mkdir(root,'results')
scans = {'SKYRA_1', 'SKYRA_2', 'SKYRA_3', 'PRISMA_1', 'PRISMA_2'};
params = {'fa'};%,'awf'};

subjects = {'Reproducibility_2020_11_29_mor_shanea','Reproducibility_2020_12_01_eve_ricardo','Reproducibility_2020_12_04_eve_ralph'};


imno = [];
immp = [];
imrmt = [];
for s = 1%:numel(subjects)
    wm = niftiread(fullfile(root,subjects{s},'seg','jhulabels_ithr.nii.gz'));
    for p = 1:numel(params)
        for sc = 1:numel(scans)
            parami_mp(:,:,:,sc) = niftiread(fullfile(root,subjects{s},'reg2',params{p},[scans{sc},'_1_',params{p},'.nii']));
            parami_nomp(:,:,:,sc) = niftiread(fullfile(root,subjects{s},'reg2',params{p},[scans{sc},'_2_',params{p},'.nii']));
            parami_rmt(:,:,:,sc) = niftiread(fullfile(root,subjects{s},'reg2',params{p},[scans{sc},'_3_',params{p},'.nii']));

        end

        P12 = [2,4]; %% change this to define which scanners we are comapring - 1/2 for skyra1 v skyra2, 4,5 for prisma1 v prisma2, see "scans" variable
        cov_prisma_mp = abs(parami_mp(:,:,:,P12(1)) - parami_mp(:,:,:,P12(2))) .* 6 ./ abs(parami_mp(:,:,:,P12(1)) + parami_nomp(:,:,:,P12(1)) + parami_rmt(:,:,:,P12(1)) + parami_mp(:,:,:,P12(2)) + parami_nomp(:,:,:,P12(2)) + parami_rmt(:,:,:,P12(2)));
        cov_prisma_nomp = abs(parami_nomp(:,:,:,P12(1)) - parami_nomp(:,:,:,P12(2))) .* 6 ./abs(parami_mp(:,:,:,P12(1)) + parami_nomp(:,:,:,P12(1)) + parami_rmt(:,:,:,P12(1)) + parami_mp(:,:,:,P12(2)) + parami_nomp(:,:,:,P12(2)) + parami_rmt(:,:,:,P12(2)));
        cov_prisma_rmt = abs(parami_rmt(:,:,:,P12(1)) - parami_rmt(:,:,:,P12(2))) .* 6 ./ abs(parami_mp(:,:,:,P12(1)) + parami_nomp(:,:,:,P12(1)) + parami_rmt(:,:,:,P12(1)) + parami_mp(:,:,:,P12(2)) + parami_nomp(:,:,:,P12(2)) + parami_rmt(:,:,:,P12(2)));
        %keyboard
        imno=cat(4,imno,cov_prisma_nomp);
        immp=cat(4,immp,cov_prisma_mp);
        imrmt=cat(4,imrmt,cov_prisma_rmt);
        

    end
end

cov_prisma_nomp = mean(imno,4);
cov_prisma_mp = mean(immp,4);
cov_prisma_rmt = mean(imrmt,4);
% cov_prisma_nomp = wm;
% 
% [min(cov_prisma_nomp(wm==19)), max(cov_prisma_nomp(wm==19)), median(cov_prisma_nomp(wm==19)), mean(cov_prisma_nomp(wm==19))]
% [min(cov_prisma_mp(wm==19)), max(cov_prisma_mp(wm==19)), median(cov_prisma_mp(wm==19)), mean(cov_prisma_mp(wm==19))]
% [min(cov_prisma_rmt(wm==19)), max(cov_prisma_rmt(wm==19)), median(cov_prisma_rmt(wm==19)), mean(cov_prisma_rmt(wm==19))]
roi=5;
[min(cov_prisma_nomp(wm==roi)), max(cov_prisma_nomp(wm==roi)), median(cov_prisma_nomp(wm==roi)), mean(cov_prisma_nomp(wm==roi))]
[min(cov_prisma_mp(wm==roi)), max(cov_prisma_mp(wm==roi)), median(cov_prisma_mp(wm==roi)), mean(cov_prisma_mp(wm==roi))]
[min(cov_prisma_rmt(wm==roi)), max(cov_prisma_rmt(wm==roi)), median(cov_prisma_rmt(wm==roi)), mean(cov_prisma_rmt(wm==roi))]


% 
% mean(cov_prisma_nomp(wm==21))
% mean(cov_prisma_mp(wm==21))
% mean(cov_prisma_rmt(wm==21))
% 
% 
% mean(cov_prisma_nomp(wm==5))
% mean(cov_prisma_mp(wm==5))
% mean(cov_prisma_rmt(wm==5))
% 
% 
% std(cov_prisma_rmt(wm==20))

mask = niftiread(fullfile(root,subjects{s},'seg','fast_seg.nii.gz'));

cov_prisma_mp(mask==0) = 0;
cov_prisma_nomp(mask==0) = 0;
cov_prisma_rmt(mask==0) = 0;

figure;
imagesc(rot90(wm(30:100,20:110,50)))

%%
close all
figure('color','k','position',[0 0 750 250]); 
fimg = cat(2,rot90(cov_prisma_nomp(30:100,20:110,50)),rot90(cov_prisma_mp(30:100,20:110,50)),rot90(cov_prisma_rmt(30:100,20:110,50)));

% figure;
% fimg1 = rot90(cov_prisma_mp(30:100,20:110,50)) - rot90(cov_prisma_nomp(30:100,20:110,50));
% fimg2 = rot90(cov_prisma_rmt(30:100,20:110,50)) - rot90(cov_prisma_nomp(30:100,20:110,50));
% 
% sum(fimg1(:)>0)
% sum(fimg1(:)<0)

% f = cat(2, fimg1, fimg2);
% imagesc(f,[-.2, .2])



%fimg(fimg>0.03) = 0;
%fimg(fimg<0) = 0;

imagesc(fimg,[0 0.3]); 
cmap = jet;
cmap(1,:) = [0,0,0];
colormap(cmap);

axis off
axis image
set(gca,'fontsize',16)
c = colorbar;
c.Color = 'w';
axis off

    