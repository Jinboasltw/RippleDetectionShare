function [results] = f_rippleDetection(opts,data,artifactsFree)
results = [];


tfg                     = [];
tfg.hpfilter            = 'yes';
tfg.hpfreq              = 0.3;
tfg.hpfiltord           = 3;
tfg.lpfilter            = 'yes';
tfg.lpfreq              = 150;
tfg.lpfiltord           = 5;
tfg.lpinstabilityfix    = 'reduce';
EEG_initalFiltered = ft_preprocessing(tfg,EEG);
        
%% filter bp 80-120hz
data_filter = ft_preprocessing(opts.rippleDetection.cfg,EEG);

%% only detect in the NREM stage
if ~isempty(opts.sleep_stage)
    temp=load(opts.sleep_stage{opts.iSubj});
    % N2
    bin_msk_N2 = (temp.scoring==2);
    % N3
    bin_msk_N3 = (temp.scoring==3);
    % N4
    bin_msk_N4 = (temp.scoring==4);
    
    bin_msk = (bin_msk_N2+bin_msk_N3+bin_msk_N4)>0;
else
    bin_msk = ones(data_filter.sampleinfo);
end

%% rms signal
rmsWin = round(opts.rippleDetection.windowLength_ms/(1000/data_filter.fsample));
%data_filter = movmean(data_filter,rmsWin); % smoothing
for j=1:numel(opts.rippleDetection.cfg.channel)
    channelIdx = find(strcmp(data_filter.label,opts.rippleDetection.cfg.channel));
    rawSignal = data_filter.trial{1,1}(channelIdx,:);
    %rawSignal = rawSignal - median(rawSignal);
    %rmsSignal = sqrt(movmean(rawSignal .^ 2, rmsWin));
    rmsSignal = f_rms(rawSignal, rmsWin, rmsWin-1, 1);
    rmsSignal_smooth = f_smooth(rmsSignal, rmsWin, rmsWin-1);

    if ~isempty(opts.sleep_stage)
        %% NREM = N2,N3,N4
        % N2
        originalMask_NREM = bin_msk_N2;
        artifactFree_mask = artifactsFree.cleanupMask{1}(channelIdx,:);
        %rmsMean_NREM_ArtFree = median(rmsSignal_smooth(originalMask_NREM==1 & artifactFree_mask==1));rmsStd_NREM_ArtFree=iqr(rmsSignal_smooth(originalMask_NREM==1 & artifactFree_mask==1))/2;
        rmsMean_NREM_ArtFree = mean(rmsSignal_smooth(originalMask_NREM==1 & artifactFree_mask==1));rmsStd_NREM_ArtFree=std(rmsSignal_smooth(originalMask_NREM==1 & artifactFree_mask==1));
        %P99_threshold = prctile(rmsSignal_smooth(originalMask_NREM==1 & artifactFree_mask==1),99);
        %results.rippleDetection_initalRange(channelIdx,:,1)=(rmsSignal > (rmsMean_NREM_ArtFree + opts.rippleDetection.bottom_RMS_Threshold*rmsStd_NREM_ArtFree));
        results.rippleDetection_initalRange(channelIdx,:,1)=(rmsSignal_smooth >= (rmsMean_NREM_ArtFree + opts.rippleDetection.bottom_RMS_Threshold*rmsStd_NREM_ArtFree))&(rmsSignal_smooth <= (rmsMean_NREM_ArtFree + opts.rippleDetection.up_RMS_Threshold*rmsStd_NREM_ArtFree));
        %results.rippleDetection_initalRange(channelIdx,:,1)=(rmsSignal_smooth >= (rmsMean_NREM_ArtFree + opts.rippleDetection.bottom_RMS_Threshold*rmsStd_NREM_ArtFree));
        results.rippleDetection_initalRange_durationFilter(channelIdx,:,1) = f_stateDurationCheck(results.rippleDetection_initalRange(channelIdx,:,1),[opts.rippleDetection.bottom_durationThreshold*round(1000/data_filter.fsample),opts.rippleDetection.up_durationThreshold*round(1000/data_filter.fsample)]);
        [results.rippleDetection_initalRange_durationFilter_cycleFilter(channelIdx,:,1),results.rippleDetection_maxPosPeakIdx{channelIdx,1}] = f_localPeak_detection(EEG.trial{1,1}(channelIdx,:),data_filter.trial{1,1}(channelIdx,:),results.rippleDetection_initalRange_durationFilter(channelIdx,:,1),opts.rippleDetection.minCycles);
        results.rippleDetection_rippleTrlDef{channelIdx} = f_maxPeak2trlDef(data_filter,EEG,channelIdx,results.rippleDetection_maxPosPeakIdx{channelIdx,1},opts.rippleDetection.periRipplePosPeak,originalMask_NREM.*artifactFree_mask);
    else
        originalMask_full = bin_msk;
        artifactFree_mask = artifactsFree.cleanupMask{1}(channelIdx,:);
        rmsMean_full = mean(rmsSignal(originalMask_full==1 & artifactFree_mask == 1));rmsStd_full=std(rmsSignal(originalMask_full==1 & artifactFree_mask==1));
        results.rippleDetection_initalRange(channelIdx,:,1)=(rmsSignal > (rmsMean_NREM_ArtFree + opts.rippleDetection.bottom_RMS_Threshold*rmsStd_NREM_ArtFree))&(rmsSignal < (rmsMean_NREM_ArtFree + opts.rippleDetection.up_RMS_Threshold*rmsStd_NREM_ArtFree));
        results.rippleDetection_initalRange_durationFilter(channelIdx,:,1) = f_stateDurationCheck(results.rippleDetection_initalRange(channelIdx,:,1),[opts.rippleDetection.bottom_durationThreshold*round(1000/data_filter.fsample),opts.rippleDetection.up_durationThreshold*round(1000/data_filter.fsample)]);
        [results.rippleDetection_initalRange_durationFilter_cycleFilter(channelIdx,:,1),results.rippleDetection_initalRange_durationFilter_cycleFilter_record{channelIdx},results.rippleDetection_maxPosPeakIdx{channelIdx,1}] = f_localPeak_detection(EEG.trial{1,1}(channelIdx,:),data_filter.trial{1,1}(channelIdx,:),results.rippleDetection_initalRange_durationFilter(channelIdx,:,1),opts.rippleDetection.minCycles);
        results.rippleDetection_rippleTrlDef{channelIdx} = f_maxPeak2trlDef(data_filter,EEG,channelIdx,results.rippleDetection_maxPosPeakIdx{channelIdx,1},opts.rippleDetection.periRipplePosPeak,originalMask_full.*artifactFree_mask);
    end
end
end

