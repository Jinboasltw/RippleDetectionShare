function [ripples,squaredSignalNormA]= ripples_detection_excluding_artifacts_Norman(bp_signal,Fs,opts,mask)
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
minDistance = opts.rippleDetection.minDistance*round(1000/Fs);
minRippleDuration = opts.rippleDetection.bottom_durationThreshold*round(1000/Fs); % in msec
maxRippleDuration = opts.rippleDetection.up_durationThreshold*round(1000/Fs); % in msec
LPcutoff = opts.rippleDetection.LPcutoff;
% convert bp_signal to RMS signal:
absSignal=double(abs(hilbert(bp_signal))); % hilbert envelope

% Clipping the signal:
topLim=nanmean(absSignal)+th(2)*nanstd(absSignal);
absSignal(absSignal>topLim)=topLim;

% Squaring:
squaredSignal = absSignal.^2;

% Smoothing using a lowpass filter:
% FIR kaiserwin lowpass filter -
lpFilt = designfilt('lowpassfir','PassbandFrequency',LPcutoff, ...
    'StopbandFrequency',LPcutoff+10,'PassbandRipple',0.001, ...
    'StopbandAttenuation',60,'DesignMethod','kaiserwin','SampleRate',Fs);
% fvtool(lpFilt,'OverlayedAnalysis','phase')
squaredSignal = filtfilt(lpFilt,double(squaredSignal));

% Compute means and std:
squaredSignal(mask~=1) = nan;
avg = nanmean(squaredSignal);
stdev = nanstd(squaredSignal,[],2);

% Hilbert envelop (rectification):
absSignalA = double(abs(hilbert(bp_signal)));

% Squaring the signal
squaredSignalA = absSignalA.^2;

% FIR filter
squaredSignalA = filtfilt(lpFilt,double(squaredSignalA));

% norm the signal with ZSCORE method:
squaredSignalNormA = (squaredSignalA-avg)/stdev;

% find peak
[pks,locs] = findpeaks(squaredSignalNormA, 'MINPEAKHEIGHT', th(2));
ENV=abs(hilbert(bp_signal));
ENV(mask~=1)=nan;
ENV=ENV./nanmedian(ENV);
ENV = 10*log10(ENV); % to calculate ripple amplitude in dB relative to median
fprintf('\n *** find %d events based on initial threshold \n',length(locs))

% refine peak detection
counter=1;
ripples=nan(1,4);
ripples=array2table(ripples,'VariableNames',{'str','peak','fin','amplitude'});
ripples(1,:)=[];
for k=locs
    % find the starting point of the peak:
    stop=0;
    str=k;
    while ~stop && ~str==0
        if str==1
            break
        end
        str=str-1;
        if squaredSignalNormA(str)<th(1), stop=1; end
    end
    % find the ending point of the peak:
    stop=0;
    fin=k;
    while ~stop && ~(fin==numel(squaredSignalNormA))
        fin=fin+1;
        if squaredSignalNormA(fin)<th(1), stop=1; end
    end
    % check mask overlap rate
    maskoverlapRate = mean(isnan(ENV(str:fin)));
    
    if maskoverlapRate == 0
        % Detect negative peak position for each ripple (closest to ripple's power peak)
        minIndex = [];
        [~,minpos] = findpeaks(-double(bp_signal(str:fin)));
        [~,maxamp] = max(double(ENV(str:fin)));
        minpos=minpos-1; 
        maxamp=maxamp-1;
        [~,tmp] = min(abs(minpos-maxamp));
        minIndex=minpos(tmp);
        peakPosition = min((str + minIndex),numel(bp_signal));

        if ~isempty(peakPosition)
            try
                ripples(counter,:)=array2table([str,peakPosition, fin, ENV(peakPosition)]);
            catch
                disp(ripples);
                fprintf('Error has occured in event # %d \n',counter);
            end
            counter=counter+1;
        end
    else
        fprintf('\n mask overlap has been found, omit location # %d',k);
    end
end
fprintf('\n After detection by thresholding: %s  events.',num2str(size(ripples,1)));
if isempty(ripples),return; end

% Merge ripples if inter-ripple period is less than minDistance:
ripples_edit=ripples;
rej=zeros(size(ripples,1),1);
for k = 2:size(ripples,1)
    if (ripples.peak(k)-ripples.peak(k-1)) < minDistance
        % Merge
        ripples_edit.fin(k-1) = ripples.fin(k);
        rej(k)=1;
    end
end
if any(rej), ripples_edit(find(rej),:)=[]; end
ripples=ripples_edit;
disp(['After ripple merge: ' num2str(size(ripples,1)) ' events.']);
if isempty(ripples),return; end

% duration test:
duration = ripples.fin-ripples.str;
ripples(duration<minRippleDuration,:) = [];
disp(['After min duration test: ' num2str(size(ripples,1)) ' events.']);
duration = ripples.fin-ripples.str;
ripples(duration>maxRippleDuration,:) = [];
disp(['After max duration test: ' num2str(size(ripples,1)) ' events.']);
end