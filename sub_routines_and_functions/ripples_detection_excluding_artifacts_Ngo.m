function [ripples,rmsSignal_smooth,bp_signal_onePass_filtered,minalThreshod,withinEvent_Pospeak,withinEvent_Negpeak]= ripples_detection_excluding_artifacts_Ngo(bp_signal,Fs,opts,nrem_mask,mask)
% Detecting ripples based on Stark et al. 2014 (Neuron) method.
% input:
% signal = normalized squared bandpassed LFP signal from the hippocampus
%          (assuming sampling rate of 500Hz)
% BP = raw bandpass signal
% t = corresponding time in sec
% th = threshold in stdev [onset/offset peak]
% minDistance = in sec
% minRippleDuration = min ripple duration in Sec
% maxRippleDuration = max ripple duration in Sec
% noise_ch = channel to exclude EMG artefacts identifed as rippples
% IED_ch = channel to exclude IEDs that were identifed as rippples
% Fs = sampling rate
%
% output:
% ripples = output table with ripple-timing and other features (start,
% peak, end, amplitude, etc.)
% ripples_stat = statistics of excluded events
%
% Author: Jinbo Zhang 22/06/21

if size(bp_signal,2)<size(bp_signal,1),bp_signal=bp_signal'; end
th(1) = opts.rippleDetection.bottom_RMS_Threshold;
th(2) = opts.rippleDetection.up_RMS_Threshold;
%minDistance = opts.rippleDetection.minDistance*round(1000/Fs);
minRippleDuration = opts.rippleDetection.bottom_durationThreshold*round(1000/Fs); % in msec
maxRippleDuration = opts.rippleDetection.up_durationThreshold*round(1000/Fs); % in msec
%LPcutoff = opts.rippleDetection.LPcutoff;
% convert bp_signal to RMS signal:
bp_signal = detrend(bp_signal);
% smothing bp_signal
% one-pass moving average filtered the BP-data:
movingavgWin = 2; % 2 data points
bp_signal_onePass_filtered = smooth(bp_signal, movingavgWin, 'moving');

% cal rms signal:
rmsWin = round(opts.rippleDetection.windowLength_ms/(1000/Fs));
if mod(rmsWin,2)==0
    rmsWin=rmsWin-1;
end
[rmsSignal,~] = envelope(bp_signal,rmsWin,'rms');

% smooth rms signal
rmsSignal_smooth = smooth(rmsSignal, rmsWin, 'moving'); % must be odd, this function will correction by -1 if it is even

% avg and std within NREM and actifact-free duration
%avg = mean(rmsSignal_smooth(mask==1));
%stdev = std(rmsSignal_smooth(mask==1));
avg = median(rmsSignal_smooth(mask==1));
stdev = iqr(rmsSignal_smooth(mask==1));

% cutoff and event thresholds
cutoffThr = avg + opts.rippleDetection.up_RMS_Threshold*stdev;
evenThr = avg + opts.rippleDetection.bottom_RMS_Threshold*stdev;

% detect ripple segments:
rippleMask = (rmsSignal_smooth >=  evenThr & rmsSignal_smooth <=  cutoffThr );
rippleMask(mask~=1) = 0;
minalThreshod = evenThr;

% reject event outside duration threshold:
rippleMask = f_stateDurationCheck(rippleMask,[opts.rippleDetection.bottom_durationThreshold,opts.rippleDetection.up_durationThreshold]);

% rject event with cycle condition not fit the threshold:
% segment the signal to event chunks:
eventMarker = rippleMask(find([1,diff(rippleMask)~=0]));
durationInfo = diff(find([1,diff(rippleMask),1]));
numpotentialEvent = sum(eventMarker == 1);
potential_Seg = find(eventMarker);

% find peaks
%bp_signal_onePass_filtered(nrem_mask~=1)=0;
[PcycleCheck_allPeaksInfo_val,PcycleCheck_allPeaksInfo_loc,PcycleCheck_allPeaksInfo_W,PcycleCheck_allPeaksInfo_P]  = findpeaks(bp_signal_onePass_filtered);
[NcycleCheck_allPeaksInfo_val,NcycleCheck_allPeaksInfo_loc,NcycleCheck_allPeaksInfo_W,NcycleCheck_allPeaksInfo_P]  = findpeaks(-1*bp_signal_onePass_filtered);
%PcycleCheck_allPeaksInfo_loc_pass = PcycleCheck_allPeaksInfo_loc(rmsSignal_smooth(PcycleCheck_allPeaksInfo_loc)>=minalThreshod & rmsSignal_smooth(PcycleCheck_allPeaksInfo_loc)<=cutoffThr);
%NcycleCheck_allPeaksInfo_loc_pass = NcycleCheck_allPeaksInfo_loc(rmsSignal_smooth(NcycleCheck_allPeaksInfo_loc)>=minalThreshod & rmsSignal_smooth(NcycleCheck_allPeaksInfo_loc)<=cutoffThr);

% loop each segment for check cycle number
% rmsSignal_smooth_temporalFineStructure = bp_signal'./rmsSignal_smooth;
for i=1:numpotentialEvent
    % extract signal within this event
    temp.segSignal = [];
    temp.targetEvent_beg(i,1) = sum(durationInfo(1:(potential_Seg(i)-1)))+1;
    temp.targetEvent_end(i,1) = sum(durationInfo(1:(potential_Seg(i))));
    temp.segSignal = bp_signal_onePass_filtered(temp.targetEvent_beg(i,1):temp.targetEvent_end(i,1));
%     if temp.targetEvent_beg(i,1)>26403470
%         disp('done');
%     end
    % cheak peaks withiin this event
    withinEvent_Pospeak{i} = PcycleCheck_allPeaksInfo_loc(find(ismember(PcycleCheck_allPeaksInfo_loc,temp.targetEvent_beg(i,1):temp.targetEvent_end(i,1))));
    withinEvent_Negpeak{i} = NcycleCheck_allPeaksInfo_loc(find(ismember(NcycleCheck_allPeaksInfo_loc,temp.targetEvent_beg(i,1):temp.targetEvent_end(i,1))));
    if length(withinEvent_Pospeak)>= opts.rippleDetection.minCycles || length(withinEvent_Negpeak) >= opts.rippleDetection.minCycles
        eventMarker_passCychle(i) = 1;
    else
        eventMarker_passCychle(i) = 0;
    end
end

% loop each segment which pass the cycle check for extracting positive peak
goodEventid = find(eventMarker_passCychle);
countRipple=0;
for i=goodEventid
    % extract signal within this event
%     temp.segSignal = [];
%     temp.targetEvent_beg(i,1) = sum(durationInfo(1:(potential_Seg(i)-1)))+1;
%     temp.targetEvent_end(i,1) = sum(durationInfo(1:(potential_Seg(i))));
%     temp.segSignal = rmsSignal_smooth(temp.targetEvent_beg(i,1):temp.targetEvent_end(i,1));
    
    % find positive peaks
%     [pkval,pkloc,pkw,pkp] = findpeaks(temp.segSignal,'MinPeakHeight',avg + opts.rippleDetection.bottom_RMS_Threshold*stdev);
%     [choosePks,locSel] = max(pkval);
    pksCandidate = withinEvent_Pospeak{i};
    [~,locs] = max(bp_signal_onePass_filtered(pksCandidate));
    peaks = pksCandidate(locs);
    if length(peaks)~=1
        disp('multiple Pks found');
        peaks
    else
        countRipple = countRipple + 1;
        ripple(countRipple,1) = peaks;
    end
end

% record in table format
epochDuration = [-1,1]*1000*(1000/Fs);
str = ripple+epochDuration(1);
peak = ripple;
fin = ripple+epochDuration(2);
amplitude = rmsSignal_smooth(peak);
ripples = array2table([str,peak,fin,amplitude]);
ripples.Properties.VariableNames = {'str','peak','fin','amplitude'};
end
