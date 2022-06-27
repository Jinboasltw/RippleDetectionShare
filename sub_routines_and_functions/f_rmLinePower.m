function EEG = f_rmLinePower(data,opts)
EEG_raw               = load(data);
cfg                   = [];
cfg.feedback          = 'no';
cfg.dftfilter         = 'yes';
cfg.dftfreq           = opts.line_power;
cfg.dftreplace        = 'zero';
cfg.dftneighbourwidth = [2,2,2,2];
EEG                   = ft_preprocessing(cfg,EEG_raw);
end

