s.channel_used = {'HIPP'}; % n by 1 cell-array

%%% ================== IF HAVE SLEEP STAGING INFO ===================== %%%
s.sleep_stage = Stage_flist;

%%% ========== MAKE SURE THE FILTER PARAMETER ARE CORRECT ============= %%%
s.rippleDetection.cfg.channel          = s.channel_used;
s.rippleDetection.cfg.type             = 'fir';
s.rippleDetection.cfg.hpfilter         = 'yes';
s.rippleDetection.cfg.hpfreq           = 80;
s.rippleDetection.cfg.hpfiltord        = 3;
s.rippleDetection.cfg.lpfilter         = 'yes';
s.rippleDetection.cfg.lpfreq           = 120;
s.rippleDetection.cfg.lpfiltord        = 5;
s.rippleDetection.cfg.lpinstabilityfix = 'reduce';
%s.rippleDetection.cfg.bsfilter         = 'yes';
%s.rippleDetection.cfg.bsfreq           = [48 52];

%%% ==== MAKE SURE THE INITIAL RIPPLE REGION PARAMETER ARE CORRECT ============= %%%
s.rippleDetection.windowLength_ms                         = 20;
s.rippleDetection.bottom_RMS_Threshold                    = 2.5; % std of RMS (Ngo,2.5) or ENV (Norman,2)
s.rippleDetection.up_RMS_Threshold                        = 9; % std of RMS (Ngo,9) or ENV (Norman,4)
s.rippleDetection.bottom_durationThreshold                = round(3*(1000/s.rippleDetection.cfg.hpfreq)); % ms
s.rippleDetection.up_durationThreshold                    = 500; % ms
s.rippleDetection.minCycles                               = 3; % 3 peaks or 3 troughts for the over-threshold rms signal
s.rippleDetection.periRipplePosPeak                       = 1000; % ms, total is 2000ms

%%% ==== MAKE SURE THE FR RIPPLE REFINE PARAMETER ARE CORRECT ============= %%%
s.rippleRefine.freqProfile.filter.cfg.hpfreq              = 65;
s.rippleRefine.freqProfile.filter.cfg.hpfiltord           = 3;
s.rippleRefine.freqProfile.filter.cfg.lpfreq              = 135;
s.rippleRefine.freqProfile.filter.cfg.lpfiltord           = 5;
s.rippleRefine.freqProfile.filter.cfg.lpinstabilityfix    = 'reduce';
s.rippleRefine.freqProfile.fft.cfg.output                 = 'pow';
%s.rippleRefine.freqProfile.fft.cfg.method                 = 'wavelet';
s.rippleRefine.freqProfile.fft.cfg.method                 = 'mtmconvol'; % multitaper time-frequency transformation based on multiplication in the frequency domain
s.rippleRefine.freqProfile.fft.cfg.taper                  = 'hanning';
%s.rippleRefine.freqProfile.fft.cfg.pad                    = 'nextpow2';
%mwCycl                                                    = ceil((65:2:135) * 0.5);     %% number of cycles per frequency for morlet wavelets
%mwCycl(mwCycl < 5)                                        = 5;
%s.rippleRefine.freqProfile.fft.cfg.width                  = mwCycl;
%s.rippleRefine.freqProfile.fft.cfg.polyremoval            = 1;
s.rippleRefine.freqProfile.fft.cfg.toi                    = -0.1:0.002:0.1;
s.rippleRefine.freqProfile.fft.cfg.foi                    = 65:2:135;
s.rippleRefine.freqProfile.fft.cfg.t_ftimwin              = ones(length(s.rippleRefine.freqProfile.fft.cfg.foi),1).*0.2;
s.rippleRefine.freqProfile.fft.cfg.keeptrials             = 'yes';
s.rippleRefine.freqProfile.avgRange                       = [-0.05,0.05];% s
s.rippleRefine.freqProfile.ampDeclineThreshold            = 0.2;
s.rippleRefine.freqProfile.ampDeclineThreshold_targetFreq = [80,120];
