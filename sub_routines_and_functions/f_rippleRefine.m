function [results,keptEvent,rmEvent] = f_rippleRefine(opts,rawData,EEG,Fs,rippleRawDetection,outdir,subjid)
keptEvent = [];
rmEvent = [];
% segment data as 1s epoch:
peaks = rippleRawDetection.peak;pre_rej_count = size(peaks,1);
durationEpoch = opts.rippleDetection.periRipplePosPeak*round(1000/Fs);

tfg.trl = [peaks-round(durationEpoch/2),peaks+round(durationEpoch/2)-1,ones(size(peaks))*(-500)];
tfg.feedback = 'no';
rippleTrl = ft_redefinetrial(tfg,rawData);

fft.cfg.feedback          = 'no';
tmp_tfr = ft_freqanalysis(opts.rippleRefine.freqProfile.fft.cfg,rippleTrl);
saveflag = 0; % change to 1 if you want to save figures
if saveflag
    for i=1:size(tmp_tfr.powspctrm(:,1,:,:),1)
        H1=figure; imagesc(tmp_tfr.time, tmp_tfr.freq, squeeze(tmp_tfr.powspctrm(i,1,:,:)));
        axis xy % flip vertically
        colormap(brewermap([],"OrRd")); hBar1 = colorbar;
        ylabel(hBar1,'Absolute power (\muV^2)','FontSize',12);
        title(sprintf('Ripple Candidate # %05d',i));
        xlabel('Tims (s)');
        ylabel('Freq (Hz)');
        set_font_size_and_type;
        drawnow;
        
        % save ripple freq profile to figures (optional):
        
        figdir = fullfile(outdir,subjid);
        if ~exist(figdir,'dir')
            mkdir(figdir);
            disp('Creating Output Directory...')
        end
        set(gcf,'renderer','painters')
        saveas(H1,fullfile(figdir,sprintf('Ripple_Candidate-%05d',i)),'fig');
        export_fig(fullfile(figdir,sprintf('Ripple_Candidate-%05d',i)),'-nocrop','-png','-r300')
        close all
    end
end

% get locs for
% time
timeIndex_str = find(opts.rippleRefine.freqProfile.fft.cfg.toi==opts.rippleRefine.freqProfile.avgRange(1));
timeIndex_end = find(opts.rippleRefine.freqProfile.fft.cfg.toi==opts.rippleRefine.freqProfile.avgRange(2));
% freq
freqIndex_str = min(find(opts.rippleRefine.freqProfile.fft.cfg.foi>=opts.rippleRefine.freqProfile.ampDeclineThreshold_targetFreq(1)));
freqIndex_end = max(find(opts.rippleRefine.freqProfile.fft.cfg.foi<=opts.rippleRefine.freqProfile.ampDeclineThreshold_targetFreq(2)));
freqUnerScale_str = opts.rippleRefine.freqProfile.fft.cfg.foi(freqIndex_str);
freqUnerScale_end = opts.rippleRefine.freqProfile.fft.cfg.foi(freqIndex_end);
freqInterval = unique(diff(opts.rippleRefine.freqProfile.fft.cfg.foi));
if freqUnerScale_str > opts.rippleRefine.freqProfile.ampDeclineThreshold_targetFreq(1)
    freqIndex_str=freqIndex_str-freqInterval;
    freqIndex_end=freqIndex_end+freqInterval;
end
% main loop for calc freq profile:
rej_count = 0;
for i = 1:size(tmp_tfr.powspctrm,1)
        if peaks(i) == 11328457
            disp('done');
        end
    %freq check
    extractedEventTrial_freqProfile{i} = mean(squeeze(tmp_tfr.powspctrm(i,:,:,timeIndex_str:timeIndex_end)),2);
    extractedEventTrial_timeProfile{i} = mean(squeeze(tmp_tfr.powspctrm(i,:,freqIndex_str:freqIndex_end,:)),1);
    [PKS_freq,loc_freq,~,Prominence_freq] = findpeaks(extractedEventTrial_freqProfile{i},'MinPeakHeight',max(extractedEventTrial_freqProfile{i})/2);
    passPeakLoc = loc_freq(find(Prominence_freq./PKS_freq>opts.rippleRefine.freqProfile.ampDeclineThreshold));
    if ~isempty(passPeakLoc)
        %checkProdPeak = loc_freq(Prominence_freq./PKS_freq==max(Prominence_freq./PKS_freq));
        if any(passPeakLoc >= freqIndex_str) & any(passPeakLoc <= freqIndex_end) & (mean(diff(passPeakLoc))<=2 | length(passPeakLoc)==1)
            [~,timeLoc]=max(extractedEventTrial_timeProfile{i});
            if timeLoc>timeIndex_str+13 & timeLoc<timeIndex_end-13 % 1/3cycle, center to the midline
                rippleRawDetection.rippleDetection_rippleTrlDef(i) = 1;
            else
                rippleRawDetection.rippleDetection_rippleTrlDef(i) = 0;
            end
        else
            % edge condition check
            [~,index]=max(mean(squeeze(tmp_tfr.powspctrm(i,:,freqIndex_str:freqIndex_end,timeIndex_str:timeIndex_end)),2));
            if (abs(passPeakLoc-freqIndex_str)==1 | abs(passPeakLoc-freqIndex_end)==1) &(index==1 | index==length(mean(squeeze(tmp_tfr.powspctrm(i,:,freqIndex_str:freqIndex_end,timeIndex_str:timeIndex_end)),2)))
                [~,timeLoc]=max(extractedEventTrial_timeProfile{i});
                if timeLoc>timeIndex_str+13 & timeLoc<timeIndex_end-13
                    rippleRawDetection.rippleDetection_rippleTrlDef(i) = 1;
                else
                    rippleRawDetection.rippleDetection_rippleTrlDef(i) = 0;
                end
            else
                rippleRawDetection.rippleDetection_rippleTrlDef(i) = 0;
            end
        end
    else
        rippleRawDetection.rippleDetection_rippleTrlDef(i) = 0;
    end
    if rippleRawDetection.rippleDetection_rippleTrlDef(i) == 0
        rej_count = rej_count + 1;
        rmEvent = [rmEvent;i];
    else
        keptEvent = [keptEvent;i];
    end
    %end
end
% frequency profile check
fprintf('\n *** rejected %d / %d events based on freq and time profile check \n',rej_count,pre_rej_count)

%% Results Check
results = rippleRawDetection;
H1=figure;
subplot(2,2,1)
temp = table2array(results(results.rippleDetection_rippleTrlDef==1,2));
tfg.trl = [temp-500,temp+500,ones(size(temp))*(-500)];
clear temp
temp = tfg.trl;
tfg.feedback = 'no';
tstring = ['80-120 Hz; N = ',sprintf('%d',size(tfg.trl,1))];
clear temp
rippleTrl_ripple = ft_redefinetrial(tfg,rawData);rippleAvg = ft_timelockanalysis([],rippleTrl_ripple);
tfg             = [];
tfg.baseline    = [-0.5,-0.4];
rippleAvg       = ft_timelockbaseline(tfg,rippleAvg);
%figure('Name','Ripple average');
plot(rippleAvg.time,rippleAvg.avg,'-b')
title(tstring);
xlim([-0.5,0.5]);
%ylim([-80,70]);
xlabel('Time (s)');
ylabel('Amplitude (\muV)');

subplot(2,2,3)
tfg                     = [];
tfg.feedback            = 'no';
tfg.hpfilter            = 'yes';
tfg.hpfreq              = 0.3;
tfg.hpfiltord           = 3;
tfg.lpfilter            = 'yes';
tfg.lpfreq              = 150;
tfg.lpfiltord           = 5;
tfg.lpinstabilityfix    = 'reduce';
%tfg.bsfilter            = 'yes';
%tfg.bsfreq              = [48 52]; % line clean implemented as a
%independent step
EEG = ft_preprocessing(tfg,EEG);
temp = table2array(results(results.rippleDetection_rippleTrlDef==1,2));
tfg.trl = [temp-500,temp+499,ones(size(temp))*(-500)];
tstring = ['0.3-150 Hz; N = ',sprintf('%d',size(tfg.trl,1))];
rippleTrl_SW = ft_redefinetrial(tfg,EEG);swrAvg = ft_timelockanalysis([],rippleTrl_SW);
tfg          = [];
tfg.baseline = [-0.5,-0.4];
swrAvg       = ft_timelockbaseline(tfg,swrAvg);
%figure('Name','Sharp wave - Ripple average');
plot(swrAvg.time,swrAvg.avg,'-r')
title(tstring);
xlim([-0.5,0.5]);
%ylim([-80,70]);
xlabel('Time (s)');
ylabel('Amplitude (\muV)');

subplot(2,2,[2,4])
avgTFR = squeeze(mean(squeeze(tmp_tfr.powspctrm(keptEvent,1,:,:)),1));
%imagesc(tmp_tfr.time, tmp_tfr.freq, pow2db(avgTFR./median(avgTFR)));
imagesc(tmp_tfr.time, tmp_tfr.freq, avgTFR);
axis xy % flip vertically
colormap(brewermap([],"OrRd")); hBar1 = colorbar('southoutside');
%caxis([-2,0]);
%ylabel(hBar1,'Power (dB)','FontSize',8);
ylabel(hBar1,'Absolute power (\muV^2)','FontSize',8);
title('TFR');
xlabel('Tims (s)');
ylabel('Freq (Hz)');
set_font_size_and_type;
drawnow;

% save ripple detection figures (optional):
saveflag = 1; % change to 1 if you want to save figures
if saveflag
    figdir = outdir;
    if ~exist(figdir,'dir')
        mkdir(figdir);
        disp('Creating Output Directory...')
    end
    set(gcf,'renderer','painters')
    saveas(H1,fullfile(figdir,'rippleAvg'),'fig');
    export_fig(fullfile(figdir,'rippleAvg'),'-nocrop','-png','-r300')
    close all
end
end

