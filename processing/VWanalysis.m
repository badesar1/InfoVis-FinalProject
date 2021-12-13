clear
close all
clc
root = '/Volumes/Research/fieree01lab/labspace/Santiago/Preprocessed_volunteers/Reproducibility_2020_12_04_eve_ralph';
mkdir(root,'results')
scans = {'SKYRA_1', 'SKYRA_2', 'SKYRA_3', 'PRISMA_1', 'PRISMA_2'};
params = {'md'};%,'awf'};

wm = niftiread(fullfile(root,'seg','jhulabels_ithr.nii.gz'));
gm = niftiread(fullfile(root,'seg','aparc+aseg.nii'));

for p = 1:numel(params)
    wmdata = [];
    gmdata = [];
    groups_wm = cell(0);
    scanner_wm = cell(0);
    groups_gm = cell(0);
    scanner_gm = cell(0);
    for sc = 1:numel(scans)
        parami_mp = niftiread(fullfile(root,'reg2',params{p},[scans{sc},'_1_',params{p},'.nii']));
        parami_nomp = niftiread(fullfile(root,'reg2',params{p},[scans{sc},'_2_',params{p},'.nii']));
        parami_rmt = niftiread(fullfile(root,'reg2',params{p},[scans{sc},'_3_',params{p},'.nii']));
        
        parami_mp(parami_mp>3) = NaN;
        parami_nomp(parami_nomp>3) = NaN;
        parami_rmt(parami_rmt>3) = NaN;
        
        wm_mp(:,sc) = parami_mp(wm==5); 
        wm_nomp(:,sc) = parami_nomp(wm==5);
        wm_rmt(:,sc) = parami_rmt(wm==5);
        
        gm_mp(:,sc) = parami_mp(gm==1030);
        gm_nomp_ = parami_nomp(gm==1030);
        gm_nomp_(gm_nomp_>10) = NaN;
        gm_nomp(:,sc) = gm_nomp_;
        gm_rmt(:,sc) = parami_rmt(gm==1030);
        
        wmdata = cat(1,wmdata,wm_mp(:,sc),wm_nomp(:,sc),wm_rmt(:,sc));
        gmdata = cat(1,gmdata,gm_mp(:,sc),gm_nomp(:,sc),gm_rmt(:,sc));
        wmdatasize = sum(sum(sum(wm==5)));
        gmdatasize = sum(sum(sum(gm==1030)));
        
        groups1 = cell(size(wm_mp,1),1);
        groups1(:) = {'MP'};
        groups2 = cell(size(wm_mp,1),1);
        groups2(:) = {'noMP'};
        groups3 = cell(size(wm_mp,1),1);
        groups3(:) = {'RMT'};
        groups_wm = cat(1,groups_wm,groups1,groups2,groups3);
        scanner1 = cell(wmdatasize*3,1);
        scanner1(:) = {scans{sc}};
        scanner_wm = cat(1,scanner_wm, scanner1);
        
        groups1 = cell(size(gm_mp,1),1);
        groups1(:) = {'MP'};
        groups2 = cell(size(gm_mp,1),1);
        groups2(:) = {'noMP'};
        groups3 = cell(size(gm_mp,1),1);
        groups3(:) = {'RMT'};
        groups_gm = cat(1,groups_gm,groups1,groups2,groups3);
        scanner1 = cell(gmdatasize*3,1);
        scanner1(:) = {scans{sc}};
        scanner_gm = cat(1,scanner_gm, scanner1);
    end
    
    figure('color','w')
    c = 1;
    for sci = 1:numel(scans)
        for scj = 1:numel(scans)
            ax = subplot(5,5,c);
            %ax.ColorOrder = flip(lines(3),1);
            rho = corr(wm_mp(:,sci),wm_mp(:,scj));
            disp([num2str(rho),' - ',scans{sci},' v ',scans{scj}])
            pl = plot(wm_nomp(:,sci),wm_nomp(:,scj),'r.',wm_mp(:,sci),wm_mp(:,scj),'b.',wm_rmt(:,sci),wm_rmt(:,scj),'g.'); hold on
            set(gca,'fontsize',16)
            pl(1).MarkerSize = 1;
            pl(2).MarkerSize = 1;
            pl(3).MarkerSize = 1;

            r = refline(1,0);
            r.Color = 'k';
            
            if sci == 5
                xlabel(scans{scj})
            end
            if scj == 1
                ylabel(scans{sci})
            end
            if sci == 1 && scj == 1
                legend('no MP','magnitude MP','RMT')
            end
            c = c + 1;
        end
    end
    suptitle([upper(params{p}) ' - White Matter'])
    %saveas(gcf,fullfile(root,'results',['VWcorr_WM_',params{p}]));
    %close all
    
    figure('color','w')
    c = 1;
    for sci = 1:numel(scans)
        for scj = 1:numel(scans)
            ax = subplot(5,5,c);
            %ax.ColorOrder = flip(lines(3),1);
            rho = corr(gm_mp(:,sci),gm_mp(:,scj));
            disp([num2str(rho),' - ',scans{sci},' v ',scans{scj}])
            pl = plot(gm_nomp(:,sci),gm_nomp(:,scj),'r.',gm_mp(:,sci),gm_mp(:,scj),'b.',gm_rmt(:,sci),gm_rmt(:,scj),'g.'); hold on
            set(gca,'fontsize',16)
            pl(1).MarkerSize = 1;
            pl(2).MarkerSize = 1;
            pl(3).MarkerSize = 1;

            r = refline(1,0);
            r.Color = 'k';
            
            if sci == 5
                xlabel(scans{scj})
            end
            if scj == 1
                ylabel(scans{sci})
            end
            if sci == 1 && scj == 1
                legend('no MP','magnitude MP','RMT')
            end
            c = c + 1;
        end
    end
    suptitle([upper(params{p}) ' - Gray Matter'])
    %saveas(gcf,fullfile(root,'results',['VWcorr_GM_',params{p}]));
    %close all
   
figure('color','w')
%col = [0 0 1; 1 0 0; 0 1 0];
h = boxplot(wmdata,{groups_wm, scanner_wm}, 'plotstyle','traditional','colors',lines(5),'factorgap',[5 0]);
set(gca,'fontsize',16)
hLegend = legend(findall(gca,'Tag','Box'), flip({'Magnitude MP','No denoising','RMT'}),'FontSize',16);
set(h,{'linew'},{2})
ylabel(upper(params{p}));
ax = gca;
xticks(2:3.7:128.5)
%xticks(2.5:4.625:128.5)
ax.XGrid = 'on';
set(gca,'XTickLabel',scans,'LineWidth',1,'FontSize',16);
title('White Matter');
%saveas(gcf,fullfile(root,'results',['boxplot_WM_',params{p}]));
%close all

figure('color','w')
%col = [0 0 1; 1 0 0; 0 1 0];
h = boxplot(gmdata,{groups_gm, scanner_gm}, 'plotstyle','traditional','colors',lines(5),'factorgap',[5 0]);
set(gca,'fontsize',16)
hLegend = legend(findall(gca,'Tag','Box'), flip({'Magnitude MP','No denoising','RMT'}),'FontSize',16);
set(h,{'linew'},{2})
ylabel(upper(params{p}));
ax = gca;
xticks(2:3.7:128.5)
%xticks(2.5:4.625:128.5)
ax.XGrid = 'on';
set(gca,'XTickLabel',scans,'LineWidth',1,'FontSize',16);
title('Gray Matter');
%saveas(gcf,fullfile(root,'results',['boxplot_GM_',params{p}]));
%close all


    
end