clear
close all
clc

subjects = {'Reproducibility_2020_11_29_mor_shanea','Reproducibility_2020_12_01_eve_ricardo','Reproducibility_2020_12_04_eve_ralph'};
root = '/Volumes/Research/fieree01lab/labspace/Santiago/Preprocessed_volunteers';
scans = {'SKYRA_1', 'SKYRA_2', 'SKYRA_3', 'PRISMA_1', 'PRISMA_2'};
%params = {'ad','rd','ad','fa','mk','rk','ak','awf'};
params = {'rotinv/S02','rotinv/S22', 'rotinv/S03','rotinv/S23','md','ad','rd','fa','mk','ak','rk','awf'};

row=1;
data = zeros(numel(subjects), numel(scans), numel(params), 2, 3, 3);
%data = cell(0);
figure('color','w');
s=3;
for p = 1:numel(params)
    subplot(3,4,p)
    for sc = 1:numel(scans)
        scdir = fullfile(root,subjects{s},scans{sc},'MB_recon','LTE');
        wm = niftiread(fullfile(root, subjects{s},'reg2',[scans{sc},'_b02MPR'],'jhulabels.nii'));
        gm = niftiread(fullfile(root, subjects{s},'reg2',[scans{sc},'_b02MPR'],'aparc+aseg.nii.gz'));
        allgm = gm == 1030;
        allwm = wm == 20 | wm == 19;

        nodn = niftiread(fullfile(scdir,'designer_noMP_noRic',[params{p},'.nii']));
        mp = niftiread(fullfile(scdir,'designer_MP',[params{p},'.nii']));
        rmt = niftiread(fullfile(scdir,'designer_complexMP',[params{p},'.nii']));
            

        L = lines(3);
        if (sc == 1 || sc == 2 || sc==3) 
        if sc == 1
            l = '-';
        elseif sc ==2
            l = '--';
        elseif sc ==3
            l = ':';
        end

        [N,edges] = histcounts(nodn(allwm),30, 'Normalization','probability'); hold on
        edges = edges(2:end) - (edges(2)-edges(1))/2;
        plot(edges, smooth(N),'color',L(1,:),'linewidth',1,'LineStyle',l);
        [N,edges] = histcounts(mp(allwm),30, 'Normalization','probability');
        edges = edges(2:end) - (edges(2)-edges(1))/2;
        plot(edges, smooth(N),'color',L(2,:),'linewidth',1,'LineStyle',l);
        [N,edges] = histcounts(rmt(allwm),30, 'Normalization','probability');
        edges = edges(2:end) - (edges(2)-edges(1))/2;
        plot(edges, smooth(N),'color',L(3,:),'linewidth',1,'LineStyle',l);
%             
%             histogram(nodn(allwm),30,'Normalization','pdf','DisplayStyle','stairs','linewidth',1,'EdgeColor',l,'LineStyle','-'); hold on
%             histogram(mp(allwm),30,'Normalization','pdf','DisplayStyle','stairs','linewidth',1,'EdgeColor',l,'LineStyle','--');
%             histogram(rmt(allwm),30,'Normalization','pdf','DisplayStyle','stairs','linewidth',1,'EdgeColor',l,'LineStyle',':');
%            
%             
            h1 = xline(nanmedian(nodn(allwm)),'color',L(1,:),'LineStyle',l);
            h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
            h2 = xline(nanmedian(mp(allwm)),'color',L(2,:),'LineStyle',l);
            h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
            h3 = xline(nanmedian(rmt(allwm)),'color',L(3,:),'LineStyle',l);
            h3.Annotation.LegendInformation.IconDisplayStyle = 'off';
            %title([upper(scans{sc})])
            if p==1
                xlabel('R_0 b1000')
            elseif p ==2
                xlabel('R_2 b1000')
            elseif p==3
                xlabel('R_0 b2000')
            elseif p==4
                xlabel('R_2 b2000')
            else
            xlabel(upper(params{p}))
            end
            
           
            set(gca,'fontsize',16)
            %keyboard
            end
          
 


            
       

    end
end
h = gca;
%legend('no denoising - PRISMA_1','MPPCA - PRISMA_1','RMT - SKYRA_1','no denoising - PRISMA_2','MPPCA - PRISMA_2','RMT - PRISMA_2');

legend('no denoising - SKYRA_1','MPPCA - SKYRA_1','RMT - SKYRA_1','no denoising - SKYRA_2','MPPCA - SKYRA_2','RMT - SKYRA_2','no denoising - PRISMA_3','MPPCA - PRISMA_3','RMT - PRISMA_3')
legend boxoff
%% medians
%close all
sclabel = {'P_1 v P_2','S_1 v S_2','P_1 v S_3','P_2 v S_3','S_1 v S_3','S_2 v S_3'};
med = data(:,:,:,:,2,:);

for p = 1:numel(params)
%     cov_prisma(p,:) = getcov(data,4,5,p);
%     cov_skyra(p,:) = getcov(data,1,2,p);
%     cov_p1s3_te(p,:) = getcov(data,3,4,p);
%     cov_p2s3_te(p,:) = getcov(data,3,5,p);
%     cov_s1s3_sc(p,:) = getcov(data,1,3,p);
%     cov_s2s3_sc(p,:) = getcov(data,2,3,p);
        
    
%     cov_prisma = mean(abs(med(:,4,:) - med(:,5,:)) ./ mean(med(:,[4,5],:),[2,3]));
%     cov_skyra = mean(abs(med(:,1,:) - med(:,2,:)) ./ mean(med(:,[1,2],:),[2,3]));
%     cov_p1s3_te = mean(abs(med(:,3,:) - med(:,4,:)) ./ mean(med(:,[3,4],:),[2,3]));
%     cov_p2s3_te = mean(abs(med(:,3,:) - med(:,5,:)) ./ mean(med(:,[3,5],:),[2,3]));
%     cov_s1s3_sc = mean(abs(med(:,1,:) - med(:,3,:)) ./ mean(med(:,[1,3],:),[2,3]));
%     cov_s2s3_sc = mean(abs(med(:,2,:) - med(:,3,:)) ./ mean(med(:,[2,3],:),[2,3]));

    cov_prisma(p,:,:) = mean(abs(med(:,4,p,:,:,:) - med(:,5,p,:,:,:)) ./ mean(med(:,[4,5],p,:,:,:),[2,6])) ./ sqrt(136);
    cov_skyra(p,:,:) = mean(abs(med(:,1,p,:,:,:) - med(:,2,p,:,:,:)) ./ mean(med(:,[1,2],p,:,:,:),[2,6])) ./ sqrt(136);
    cov_p1s3_te(p,:,:) = mean(abs(med(:,3,p,:,:,:) - med(:,4,p,:,:,:)) ./ mean(med(:,[3,4],p,:,:,:),[2,6]))./ sqrt(136);
    cov_p2s3_te(p,:,:) = mean(abs(med(:,3,p,:,:,:) - med(:,5,p,:,:,:)) ./ mean(med(:,[3,5],p,:,:,:),[2,6]))./ sqrt(136);
    cov_s1s3_sc(p,:,:) = mean(abs(med(:,1,p,:,:,:) - med(:,3,p,:,:,:)) ./ mean(med(:,[1,3],p,:,:,:),[2,6]))./ sqrt(136);
    cov_s2s3_sc(p,:,:) = mean(abs(med(:,2,p,:,:,:) - med(:,3,p,:,:,:)) ./ mean(med(:,[2,3],p,:,:,:),[2,6]))./ sqrt(136);

end

%%
figure('color','w');
for p = 1:numel(params)
    cov_wm = squeeze(cat(1,cov_prisma(p,1,:),cov_skyra(p,1,:),cov_p1s3_te(p,1,:),cov_p2s3_te(p,1,:),cov_s1s3_sc(p,1,:),cov_s2s3_sc(p,1,:)));
    %cov_wm = squeeze(cat(1,cov_prisma,cov_skyra,cov_p1s3_te,cov_p2s3_te,cov_s1s3_sc,cov_s2s3_sc));

    %keyboard
    subplot(1,5,p)
    bar(cov_wm,'grouped')
    title(upper(params{p}))
    set(gca,'fontsize',16)
    ylabel('CoV - WM')
    if p == 1
       legend('no denoising','magnitude MP','RMT')
    end
    set(gca,'XTickLabel', sclabel,'LineWidth',1,'FontSize',16);
    xtickangle(45)
%    suptitle('White Matter');
end

function cov = getcov(data, sind1, sind2, p)
        scan_wm_nodn_s1 = cat(1,data{1,sind1,p,1,1}, data{1,sind2,p,1,1});
        scan_wm_mp_s1 = cat(1,data{1,sind1,p,1,2}, data{1,sind2,p,1,2});
        scan_rmt_s1 = cat(1,data{1,sind1,p,1,3}, data{1,sind2,p,1,3});
        prisma_12_wm_nodn_s2 = cat(1,data{2,sind1,p,1,1}, data{2,sind2,p,1,1});
        prisma_12_wm_mp_s2 = cat(1,data{2,sind1,p,1,2}, data{2,sind2,p,1,2});
        prisma_12_wm_rmt_s2 = cat(1,data{2,sind1,p,1,3}, data{2,sind2,p,1,3});
        prisma_12_wm_nodn_s3 = cat(1,data{3,sind1,p,1,1}, data{3,sind2,p,1,1});
        prisma_12_wm_mp_s3 = cat(1,data{3,sind1,p,1,2}, data{3,sind2,p,1,2});
        prisma_12_wm_rmt_s3 = cat(1,data{3,sind1,p,1,3}, data{3,sind2,p,1,3});
        
        figure;
        L = lines(3);
        histogram(data{1,sind1,p,1,1},30,'DisplayStyle','stairs','linewidth',1,'EdgeColor',L(1,:),'LineStyle','-'); hold on
        histogram(data{1,sind2,p,1,1},30,'DisplayStyle','stairs','linewidth',1,'EdgeColor',L(1,:),'LineStyle','--');
        histogram(data{1,sind1,p,1,2},30,'DisplayStyle','stairs','linewidth',1,'EdgeColor',L(2,:),'LineStyle','-'); hold on
        histogram(data{1,sind2,p,1,2},30,'DisplayStyle','stairs','linewidth',1,'EdgeColor',L(2,:),'LineStyle','--');
        histogram(data{1,sind1,p,1,3},30,'DisplayStyle','stairs','linewidth',1,'EdgeColor',L(3,:),'LineStyle','-'); hold on
        histogram(data{1,sind2,p,1,3},30,'DisplayStyle','stairs','linewidth',1,'EdgeColor',L(3,:),'LineStyle','--');
        legend('no denoising - scanner1','no denoising - scanner2','MP - scanner1','MP - scanner2','RMT - scanner1','RMT scanner2')
        keyboard
        
        cov_nodn_s1 = nanstd(prisma_12_wm_nodn_s1)./nanmean(cat(1,prisma_12_wm_nodn_s1, prisma_12_wm_mp_s1, prisma_12_wm_rmt_s1));
        cov_mp_s1 = nanstd(prisma_12_wm_mp_s1)./nanmean(cat(1,prisma_12_wm_nodn_s1, prisma_12_wm_mp_s1, prisma_12_wm_rmt_s1));
        cov_rmt_s1 = nanstd(prisma_12_wm_rmt_s1)./nanmean(cat(1,prisma_12_wm_nodn_s1, prisma_12_wm_mp_s1, prisma_12_wm_rmt_s1));
        
        cov_nodn_s2 = nanstd(prisma_12_wm_nodn_s2)./nanmean(cat(1,prisma_12_wm_nodn_s2, prisma_12_wm_mp_s2, prisma_12_wm_rmt_s2));
        cov_mp_s2 = nanstd(prisma_12_wm_mp_s2)./nanmean(cat(1,prisma_12_wm_nodn_s2, prisma_12_wm_mp_s2, prisma_12_wm_rmt_s2));
        cov_rmt_s2 = nanstd(prisma_12_wm_rmt_s2)./nanmean(cat(1,prisma_12_wm_nodn_s2, prisma_12_wm_mp_s2, prisma_12_wm_rmt_s2));
        
        cov_nodn_s3 = nanstd(prisma_12_wm_nodn_s3)./nanmean(cat(1,prisma_12_wm_nodn_s3, prisma_12_wm_mp_s3, prisma_12_wm_rmt_s3));
        cov_mp_s3 = nanstd(prisma_12_wm_mp_s3)./nanmean(cat(1,prisma_12_wm_nodn_s3, prisma_12_wm_mp_s3, prisma_12_wm_rmt_s3));
        cov_rmt_s3 = nanstd(prisma_12_wm_rmt_s3)./nanmean(cat(1,prisma_12_wm_nodn_s3, prisma_12_wm_mp_s3, prisma_12_wm_rmt_s3));
        
        cov_nodn = nanmean(cat(1,cov_nodn_s1,cov_nodn_s2,cov_nodn_s3));
        cov_mp = nanmean(cat(1,cov_mp_s1,cov_mp_s2,cov_mp_s3));
        cov_rmt = nanmean(cat(1,cov_rmt_s1,cov_rmt_s2,cov_rmt_s3));
        cov = cat(1, cov_nodn, cov_mp, cov_rmt);
end

% figure('color','w','position',[0 0 2000 1000]);
% for p = 1:numel(params)
%     cov_wm = squeeze(cat(2,cov_prisma(p,2,:),cov_skyra(p,2,:),cov_p1s3_te(p,2,:),cov_p2s3_te(p,2,:),cov_s1s3_sc(p,2,:),cov_s2s3_sc(p,2,:)));
%     %subplot(2,4,p)
%     bar(cov_wm,'grouped')
%     title(upper(params{p}))
%     set(gca,'fontsize',16)
%     ylabel('CoV - GM')
%     if p == 1
%         legend('no denoising','magnitude MP','RMT')
%     end
%     set(gca,'XTickLabel', sclabel,'LineWidth',1,'FontSize',16);
%     xtickangle(45)
% %    suptitle('Gray Matter');
% end
% 
%     
            
            
            





