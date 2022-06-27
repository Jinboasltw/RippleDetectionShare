% SWR detection algorithm.
% Detecting ripples after exclusion of amplitude based artifacts, gradient artifacts, and arousals or movement artifacts.
% The artifacts cleanup algorithm uses one main signals:
% (1) hippocampal raw data
% The detection algorithm use two main input:
% (1) hippocampal raw data
% (2) artifact-free NREM binary mask, datapoints need kept are marked as 1
% Ref. Ngo, et la., 2020, Elife

% Author: Jinbo Zhang, 2022, Liu Lab
% @ State Key Laboratory of Cognitive Neuroscience and Learning, Beijing Normal University
% @ Chinese Institute for Brain Research,Beijing
% Mail to Author:  <a href="jinbozhang@cibr.ac.cn">jinbozhang@cibr.ac.cn</a>
% Date: 2022-06-21
% Last edited: 2022-06-26

run(fullfile(pwd,'..','startup_script.m'));

% close some fieldtrip output:
global ft_default
ft_default.showcallinfo = 'no';

% load neuronal data:
EEG_flist      = getFilelist(fullfile(parentfolder,'RippleDetection','data','pat*_HIPP.mat'));
Stage_flist    = getFilelist(fullfile(parentfolder,'RippleDetection','data','pat*_HIPP_supplement.mat'));
[~,subjects,~] = cellfun(@(x) fileparts(x), EEG_flist, 'UniformOutput', false);
subjects       = cellfun(@(x) x(1:5), subjects, 'UniformOutput', false);

%% Main analysis loop:
for iSub=1:numel(subjects)
    subjid=subjects{iSub};
    s.iSubj=iSub;
    
    %% Artifacts cleaning:
    %%%%% Artifacts PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    artifactDetectionParameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % set output folder:
    tempdir=fullfile(parentfolder,'RippleDetection','temp');
    if ~exist(tempdir,'dir')
        mkdir(tempdir);
        disp('Creating Temp Directory...')
        fprintf('\n\n =============== \n Temp directory: \n %s \n ',tempdir);
    end
    files_out=fullfile(parentfolder,'RippleDetection','temp',...
        sprintf('%s_cleanMask.mat',subjid));
    
    % remove line power in raw data:
    EEG    = f_rmLinePower(EEG_flist{iSub},s);
    Signal = double(cell2mat(EEG.trial));
    Fs     = EEG.fsample;
    %     EEG = load(EEG_flist{iSub});
    %     Signal = double(cell2mat(EEG.trial));
    %     Fs = EEG.fsample;

    % Detect artifacts and get the clean mask (keep=1):
    artifacts=f_artifacts(s,EEG);
    cleanupMask=f_artifacts_aligment(s,artifacts,EEG);
    save(files_out,'cleanupMask')
    
    %% Detection prepare
    %%%%% Ripple detection PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    swrDetectionParameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % method used
    %methodChoose = 'Norman'; % Norman's Method
    methodChoose = 'Ngo'; % user Ngo's Method
    
    % filter-related parameter
    ripple_band = [s.rippleDetection.cfg.hpfreq,s.rippleDetection.cfg.lpfreq]; % in Hz
    wtype       = s.rippleRefine.freqProfile.fft.cfg.taper;
    th          = [s.rippleDetection.bottom_RMS_Threshold,s.rippleDetection.up_RMS_Threshold]; % ripple detection thresholds (in std from the mean) [onset/offset]
    
    % duration threshold
    minRippleDuration = s.rippleDetection.bottom_durationThreshold; % in msec
    maxRippleDuration = s.rippleDetection.up_durationThreshold; % in msec
    
    % set output folder:
    dataoutputfolder=fullfile(parentfolder,'RippleDetection','ripple_times',...
        sprintf('Ripple_times_%s_%g-%dstd_%d-%dHz_%d-%dms_method_%s',...
        wtype,th(1),th(2),ripple_band(1),ripple_band(2),...
        minRippleDuration,maxRippleDuration,methodChoose),subjid); % ADJUST OUTDIR
    if ~exist(dataoutputfolder,'dir')
        mkdir(dataoutputfolder);
        disp('Creating Output Directory...')
    end
    % adjust image output dir:
    outdir=fullfile(parentfolder,'RippleDetection','figures',...
        sprintf('Ripple_images_%s_%g-%dstd_%d-%dHz_%d-%dms_method_%s',...
        wtype,th(1),th(2),ripple_band(1),ripple_band(2),...
        minRippleDuration,maxRippleDuration,methodChoose),subjid);
    fprintf('\n\n =============== \n Image Output directory: \n %s \n ',outdir);
    if ~exist(outdir,'dir')
        mkdir(outdir);
        disp('Creating Image Output Directory...')
        fprintf('\n\n =============== \n Ouput directory: \n %s \n ',outdir);
    end
    
    % get mask for artifact-free and NREM stage:
    temp       = load(s.sleep_stage{iSub});
    bin_msk_N2 = (temp.scoring==2);
    bin_msk_N3 = (temp.scoring==3);
    bin_msk_N4 = (temp.scoring==4);
    
    bin_msk_NREM = double((bin_msk_N2+bin_msk_N3+bin_msk_N4)>0);
    bin_msk = double((bin_msk_N2+bin_msk_N3+bin_msk_N4)>0 & cell2mat(cleanupMask)==1);
    %bin_msk = cell2mat(cleanupMask);
    clear temp

    % band-pass filtered from 80-120 Hz:
    s.rippleDetection.cfg.feedback = 'no';
    bpSignal_temp = ft_preprocessing(s.rippleDetection.cfg,EEG);
    bpSignal = double(cell2mat(bpSignal_temp.trial));
    
    if strcmp(methodChoose,'Norman')
        disp('======> Use Norman''s approach:')
        % Low-pass filter for smoothing the ripple-band envelope:
        LPcutoff = round(mean(ripple_band)/pi); % in Hz (see Stark et al. 2014 for a similar method)
        s.rippleDetection.LPcutoff = LPcutoff;
        % Note: the cutoff frequency is 40 Hz when ripple-band frequency window is 70-180 Hz
        minDistance                   = 0.030; % in sec, between peaks, Norman's method
        s.rippleDetection.minDistance = minDistance;
        
        % find ripples:
        [ripples,squaredSignalNormA] = ripples_detection_excluding_artifacts_Norman(bpSignal,Fs,s,bin_msk);
    elseif strcmp(methodChoose,'Ngo')
        disp('======> Use Ngo''s approach:')
        [ripples,rmsSignal,bpSmooth,minalThreshod,peaks,throughts] = ripples_detection_excluding_artifacts_Ngo(bpSignal,Fs,s,bin_msk_NREM,bin_msk);
    else
        disp('======> Error: Please check input!')
        methodChoose
        return
    end
    outfilename = sprintf('%s_ripples',subjid);
    save(fullfile(dataoutputfolder,outfilename),'ripples')
    
    %% Plot detection figure:
    if strcmp(methodChoose,'Norman')
        plotUseData = squaredSignalNormA;
        close all
        channel = 'HIPP';
        H1=figure('Name',[subjid ' ripples detection Ch ' channel],'units','normalized','outerposition',[0 0 1 0.5],'Color','w');
        hold on;
        title(['Ripples Detection ' ' (ripple band envelope)']);
        %h1 = plot(T,squaredSignalNormB,'linewidth',0.5,'color',[1,0.7,0.7]); hold on;
        rippleUse = [1,2];
        plotStr = ripples.str(rippleUse(1));plotFin = ripples.fin(rippleUse(2));
        T=plotStr:plotFin;
        h1 = plot(T,plotUseData(plotStr:plotFin),'k','Linesmoothing','on');
        ylim([-2,max(plotUseData(plotStr:plotFin))+3]);
        needHighlight=1;
        [x_points,y_points]=f_shadeByMask(bin_msk(plotStr:plotFin),needHighlight,H1.CurrentAxes.YLim(1),H1.CurrentAxes.YLim(2));
        color_keep = [0,1,0];
        for i=1:size(x_points,1)
            h2 = fill(plotStr+x_points(i,:)-1,y_points(i,:),color_keep);h2.FaceAlpha=0.05;h2.EdgeAlpha=0;
        end
        ylim([y_points(i,1),y_points(i,2)]);
        clear x_points y_points
        needHighlight=0;
        [x_points,y_points]=f_shadeByMask(bin_msk(plotStr:plotFin),needHighlight,H1.CurrentAxes.YLim(1),H1.CurrentAxes.YLim(2));
        color_reject = [1,0,0];
        for i=1:size(x_points,1)
            h2 = fill(plotStr+x_points(i,:)-1,y_points(i,:),color_reject);h2.FaceAlpha=0.05;h2.EdgeAlpha=0;
        end
        h0 = plot(T,ones(size(T))*th(2),'b--','Linewidth',0.5);
        ylim([y_points(i,1),y_points(i,2)]);
        xlim([plotStr-100,plotFin+100]);
        if ~isempty(ripples)
            selectedRipples = ripples(rippleUse(1):rippleUse(2),:);
            h3=scatter(selectedRipples.peak,ones(size(selectedRipples,1),1)*-2,30,'bo','fill'); hold on
            for j=1:length(selectedRipples.peak)
                xline(selectedRipples.peak(j),':'); hold on
            end
        end
        ylabel(sprintf('%d-%dHz Amplitude^2 (Z-score)',ripple_band(1),ripple_band(2)))
        xlim([T(1),T(end)])
        ylim([-2.1,y_points(i,2)])
        legend([h1 h2 h3 h0],{'Ripple Band Amplitude (zscore)','Mask Mark (Red = Reject; Green = Keep)','Ripple','Detection thr.'}); legend boxoff;
        set_font_size_and_type;
        drawnow;
    elseif strcmp(methodChoose,'Ngo')
        plotUseData = rmsSignal;
        close all
        channel = 'HIPP';
        H1=figure('Name',[subjid ' ripples detection Ch ' channel],'units','normalized','outerposition',[0 0 1 0.5],'Color','w');
        hold on;
        title(['Ripples Detection ' ' (ripple band envelope)']);
        %h1 = plot(T,squaredSignalNormB,'linewidth',0.5,'color',[1,0.7,0.7]); hold on;
        rippleUse = [2,2];
        plotStr = ripples.str(rippleUse(1));plotFin = ripples.fin(rippleUse(2));
        T=plotStr:plotFin;
        h1 = plot(T,plotUseData(plotStr:plotFin),'k','Linesmoothing','on');hold on
        h4 = plot(T,bpSmooth(plotStr:plotFin)+10,'b','Linesmoothing','on');hold on
        ylim([-2,max(plotUseData(plotStr:plotFin))+3]);
        needHighlight=1;
        [x_points,y_points]=f_shadeByMask(bin_msk(plotStr:plotFin),needHighlight,H1.CurrentAxes.YLim(1),H1.CurrentAxes.YLim(2));
        color_keep = [0,1,0];
        for i=1:size(x_points,1)
            h2 = fill(plotStr+x_points(i,:)-1,y_points(i,:),color_keep);h2.FaceAlpha=0.05;h2.EdgeAlpha=0;
        end
        ylim([y_points(i,1),y_points(i,2)]);
        clear x_points y_points
        needHighlight=0;
        [x_points,y_points]=f_shadeByMask(bin_msk(plotStr:plotFin),needHighlight,H1.CurrentAxes.YLim(1),H1.CurrentAxes.YLim(2));
        color_reject = [1,0,0];
        for i=1:size(x_points,1)
            h2 = fill(plotStr+x_points(i,:)-1,y_points(i,:),color_reject);h2.FaceAlpha=0.05;h2.EdgeAlpha=0;
        end
        h0 = plot(T,ones(size(T))*minalThreshod,'r--','Linewidth',0.5);
        ylim([y_points(i,1),y_points(i,2)]);
        xlim([plotStr-100,plotFin+100]);
        if ~isempty(ripples)
            selectedRipples = ripples(rippleUse(1):rippleUse(2),:);
            h3=scatter(selectedRipples.peak,ones(size(selectedRipples,1),1)*-2,30,'bo','fill'); hold on
            for j=1:length(selectedRipples.peak)
                xline(selectedRipples.peak(j),':'); hold on
            end
        end
        ylabel(sprintf('%d-%dHz Amplitude^2',ripple_band(1),ripple_band(2)))
        xlim([T(1),T(end)])
        ylim([-2.1,y_points(i,2)+20])
        legend([h1 h4 h2 h3 h0],{'Ripple Band RMS','Band-passed Signal 80-120Hz','Mask Mark (Red = Reject; Green = Keep)','Ripple','Detection thr.'}); legend boxoff;
        set_font_size_and_type;
        drawnow;
    else
        methodChoose
        disp('======> Error: Please check input!')
        return
    end
    
    
    % save ripple detection figures:
    saveflag = 1; % change to 1 if you want to save figures
    if saveflag
        figdir = outdir;
        if ~exist(figdir,'dir')
            mkdir(figdir);
            disp('Creating Image Output Directory...')
        end
        set(gcf,'renderer','painters')
        saveas(H1,fullfile(figdir,get(gcf,'name')),'fig');
        %export_fig(fullfile(figdir,get(gcf,'name')),'-nocrop','-png','-r150','-transparent')
        export_fig(fullfile(figdir,get(gcf,'name')),'-nocrop','-png','-r300')
        close all
    end
    
    % frequency profile check
    [results,keptEvent,rmEvent] = f_rippleRefine(s,bpSignal_temp,EEG,Fs,ripples,outdir,subjid);
    
    %% save final results
    outfilename = [outfilename,'_freq_time_profile_checked'];
    save(fullfile(dataoutputfolder,outfilename),'results');
    writetable(results,fullfile(dataoutputfolder,[outfilename,'.csv']));
    writematrix(keptEvent,fullfile(dataoutputfolder,[outfilename,'_keptEv.csv']));
    writematrix(rmEvent,fullfile(dataoutputfolder,[outfilename,'_rmEv.csv']));
    
    %% unit test
    %unittest_rippleDetection
end