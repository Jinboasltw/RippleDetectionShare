s.channel_used = {'HIPP'}; % n by 1 cell-array

%%% ================== IF HAVE SLEEP STAGING INFO ===================== %%%
s.sleep_stage = Stage_flist;

%%% ========== MAKE SURE THE FILTER PARAMETER ARE CORRECT ============= %%%
% line power
s.line_power                                          = [50,100,150,200];

% artifacts: amplitude
s.artifacts.amplitude.cfg.channel                     = s.channel_used;
s.artifacts.amplitude.cfg.type                        = 'fir';
s.artifacts.amplitude.cfg.bpfilter                    = 'yes';
s.artifacts.amplitude.cfg.bpfreq                      = [0.3,150];
s.artifacts.amplitude.cfg.bpfiltord                   = 4;
s.artifacts.amplitude.threshold                       = 750; % muV, >

% artifacts: gradient
s.artifacts.gradient.cfg.channel                      = s.channel_used;
s.artifacts.gradient.cfg.type                         = 'fir';
s.artifacts.gradient.cfg.bpfilter                     = 'yes';
s.artifacts.gradient.cfg.bpfreq                       = [0.3,150];
s.artifacts.gradient.cfg.bpfiltord                    = 4;
s.artifacts.gradient.threshold                        = 6; % IQR (Interquartile Range) of changes of amplitude of t and t-1, outof this range were marked as artifacts

% artifacts: arousals or movement
s.artifacts.arousals_movement.cfg.channel             = s.channel_used;
s.artifacts.arousals_movement.type                    = 'fir';
s.artifacts.arousals_movement.cfg.hpfilter            = 'yes';
s.artifacts.arousals_movement.cfg.hpfreq              = 150;
s.artifacts.arousals_movement.cfg.hpfiltord           = 4;
s.artifacts.arousals_movement.windowLength_ms         = 100; % 100ms
s.artifacts.arousals_movement.threshold               = 4; % IQR (Interquartile Range) of changes of amplitude of t and t-1, outof this range were marked as artifacts
s.artifacts.arousals_movement.over_threshold_duration = 100; % 100ms, exceeded the corresponding threshold for at least 100 ms

%%% ========= MAKE SURE THE PADDING PARAMETER ARE CORRECT ============= %%%
s.pad_duration  = 250; % unit: ms
s.minWithoutArt = 3000; % unit: ms