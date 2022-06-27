function [results] = f_artifacts(opts,EEG)
results = [];

fdNames = fieldnames(opts.artifacts);
for i=1:numel(fdNames)
    switch fdNames{i}
        case 'amplitude'
            fprintf('artifacts detection: amplitude.doing!\n')
            try
                opts.artifacts.amplitude.cfg.feedback          = 'no';
                data_filter = ft_preprocessing(opts.artifacts.amplitude.cfg,EEG);
                for j=1:numel(opts.artifacts.amplitude.cfg.channel)
                    channelIdx = find(strcmp(data_filter.label,opts.artifacts.amplitude.cfg.channel));
                    results.artifacts_amplitude_bin_mask(channelIdx,:)=abs(data_filter.trial{1,1}(channelIdx,:))>opts.artifacts.amplitude.threshold;
                    results.artifacts_amplitude_report{channelIdx,1}=tabulate(results.artifacts_amplitude_bin_mask(channelIdx,:));
                end
                fprintf('artifacts detection: amplitude.success!\n')
                clear data_filter
            catch
                fprintf('artifacts detection: amplitude.fail!\n')
            end
        case 'gradient'
            fprintf('artifacts detection: gradient.doing!\n')
            try
                opts.artifacts.amplitude.cfg.feedback = 'no';
                data_filter = ft_preprocessing(opts.artifacts.amplitude.cfg,EEG);
                if ~isempty(opts.sleep_stage)
                    temp=load(opts.sleep_stage{opts.iSubj});
                    % N2
                    bin_msk_N2 = (temp.scoring==2);
                    % N3
                    bin_msk_N3 = (temp.scoring==3);
                    % N4
                    bin_msk_N4 = (temp.scoring==4);
                end
                for j=1:numel(opts.artifacts.gradient.cfg.channel)
                    channelIdx = find(strcmp(data_filter.label,opts.artifacts.gradient.cfg.channel));
                    deltaArray = diff(data_filter.trial{1,1}(channelIdx,:));
                    if ~isempty(opts.sleep_stage)
                        % N2
                        deltaArray_median_N2 = median(deltaArray(1,bin_msk_N2(2:end)));deltaArray_iqr_N2=iqr(deltaArray(1,bin_msk_N2(2:end))); % use 2:end because, diff function is calculated by T_(n+1) - T_(n)
                        results.artifacts_gradient_bin_mask(channelIdx,:,1)=(deltaArray > (deltaArray_median_N2 + opts.artifacts.gradient.threshold*deltaArray_iqr_N2))|(deltaArray < (deltaArray_median_N2 - opts.artifacts.gradient.threshold*deltaArray_iqr_N2));
                        results.artifacts_gradient_bin_mask(channelIdx,:,1) = results.artifacts_gradient_bin_mask(channelIdx,:,1) .* bin_msk_N2(2:end);
                        % N3
                        deltaArray_median_N3 = median(deltaArray(1,bin_msk_N3(2:end)));deltaArray_iqr_N3=iqr(deltaArray(1,bin_msk_N3(2:end)));
                        results.artifacts_gradient_bin_mask(channelIdx,:,2)=(deltaArray > (deltaArray_median_N3 + opts.artifacts.gradient.threshold*deltaArray_iqr_N3))|(deltaArray < (deltaArray_median_N3 - opts.artifacts.gradient.threshold*deltaArray_iqr_N3));
                        results.artifacts_gradient_bin_mask(channelIdx,:,2)= results.artifacts_gradient_bin_mask(channelIdx,:,2) .* bin_msk_N3(2:end);
                        % N4
                        deltaArray_median_N4 = median(deltaArray(1,bin_msk_N4(2:end)));deltaArray_iqr_N4=iqr(deltaArray(1,bin_msk_N4(2:end)));
                        results.artifacts_gradient_bin_mask(channelIdx,:,3)=(deltaArray > (deltaArray_median_N4 + opts.artifacts.gradient.threshold*deltaArray_iqr_N4))|(deltaArray < (deltaArray_median_N4 - opts.artifacts.gradient.threshold*deltaArray_iqr_N4));
                        results.artifacts_gradient_bin_mask(channelIdx,:,3)= results.artifacts_gradient_bin_mask(channelIdx,:,3) .* bin_msk_N4(2:end);
                        % combine
                        results.artifacts_gradient_bin_mask_fusion(channelIdx,:) = (sum(results.artifacts_gradient_bin_mask,3)>0);
                        results.artifacts_gradient_report{channelIdx,1} = tabulate(sum(results.artifacts_gradient_bin_mask,3)>0);
                        results.artifacts_gradient_report_N2{channelIdx,1} = tabulate(results.artifacts_gradient_bin_mask(channelIdx,:,1));
                        results.artifacts_gradient_report_N3{channelIdx,1} = tabulate(results.artifacts_gradient_bin_mask(channelIdx,:,2));
                        results.artifacts_gradient_report_N4{channelIdx,1} = tabulate(results.artifacts_gradient_bin_mask(channelIdx,:,3));
                    else
                        deltaArray_median_all = median(deltaArray(1,:));deltaArray_iqr_all=iqr(deltaArray(1,:));
                        results.artifacts_gradient_bin_mask_fusion(channelIdx,:,1)=(deltaArray > (deltaArray_median_all + opts.artifacts.gradient.threshold*deltaArray_iqr_all))|(deltaArray < (deltaArray_median_all - opts.artifacts.gradient.threshold*deltaArray_iqr_all));
                        results.artifacts_gradient_report{channelIdx,1} = tabulate(results.artifacts_gradient_bin_mask_fusion(channelIdx,:,1));
                    end
                end
                fprintf('artifacts detection: gradient.success!\n')
                clear data_filter
            catch
                fprintf('artifacts detection: gradient.fail!\n')
            end
        case 'arousals_movement'
            fprintf('artifacts detection: arousals movements.doing!\n')
            try
                opts.artifacts.arousals_movement.cfg.feedback          = 'no';
                data_filter = ft_preprocessing(opts.artifacts.arousals_movement.cfg,EEG);
                if ~isempty(opts.sleep_stage)
                    temp=load(opts.sleep_stage{opts.iSubj});
                    % N2
                    bin_msk_N2 = (temp.scoring==2);
                    % N3
                    bin_msk_N3 = (temp.scoring==3);
                    % N4
                    bin_msk_N4 = (temp.scoring==4);
                end
                rmsWin = round(opts.artifacts.arousals_movement.windowLength_ms/(1000/EEG.fsample));
                for j=1:numel(opts.artifacts.arousals_movement.cfg.channel)
                    channelIdx = find(strcmp(data_filter.label,opts.artifacts.arousals_movement.cfg.channel));
                    rawSignal = data_filter.trial{1,1}(channelIdx,:);
                    [rmsSignal,~] = envelope(rawSignal,rmsWin,'rms');
                    if ~isempty(opts.sleep_stage)
                        % N2
                        originalMask_N2 = bin_msk_N2;
                        rmsMedian_N2 = median(rmsSignal(originalMask_N2==1));rmsMedian_irq_N2=iqr(rmsSignal(originalMask_N2==1));
                        results.artifacts_arousals_movement_bin_mask(channelIdx,:,1)=(rmsSignal > (rmsMedian_N2 + opts.artifacts.arousals_movement.threshold*rmsMedian_irq_N2))|(rmsSignal < (rmsMedian_N2 - opts.artifacts.arousals_movement.threshold*rmsMedian_irq_N2));
                        results.artifacts_arousals_movement_bin_mask(channelIdx,:,1) = f_stateDurationCheck(results.artifacts_arousals_movement_bin_mask(channelIdx,:,1),opts.artifacts.arousals_movement.over_threshold_duration*round(1000/data_filter.fsample));
                        results.artifacts_arousals_movement_bin_mask(channelIdx,:,1)= results.artifacts_arousals_movement_bin_mask(channelIdx,:,1) .* originalMask_N2;
                        % N3
                        originalMask_N3 = bin_msk_N3;
                        rmsMedian_N3 = median(rmsSignal(originalMask_N3==1));rmsMedian_irq_N3=iqr(rmsSignal(originalMask_N3==1));
                        results.artifacts_arousals_movement_bin_mask(channelIdx,:,2)=(rmsSignal > (rmsMedian_N3 + opts.artifacts.arousals_movement.threshold*rmsMedian_irq_N3))|(rmsSignal < (rmsMedian_N3 - opts.artifacts.arousals_movement.threshold*rmsMedian_irq_N3));
                        results.artifacts_arousals_movement_bin_mask(channelIdx,:,2) = f_stateDurationCheck(results.artifacts_arousals_movement_bin_mask(channelIdx,:,2),opts.artifacts.arousals_movement.over_threshold_duration*round(1000/data_filter.fsample));
                        results.artifacts_arousals_movement_bin_mask(channelIdx,:,2)= results.artifacts_arousals_movement_bin_mask(channelIdx,:,2) .* originalMask_N3;
                        % N4
                        originalMask_N4 = bin_msk_N4;
                        rmsMedian_N4 = median(rmsSignal(originalMask_N4==1));rmsMedian_irq_N4=iqr(rmsSignal(originalMask_N4==1));
                        results.artifacts_arousals_movement_bin_mask(channelIdx,:,3)=(rmsSignal > (rmsMedian_N4 + opts.artifacts.arousals_movement.threshold*rmsMedian_irq_N4))|(rmsSignal < (rmsMedian_N4 - opts.artifacts.arousals_movement.threshold*rmsMedian_irq_N4));
                        results.artifacts_arousals_movement_bin_mask(channelIdx,:,3)= f_stateDurationCheck(results.artifacts_arousals_movement_bin_mask(channelIdx,:,3),opts.artifacts.arousals_movement.over_threshold_duration*round(1000/data_filter.fsample));
                        results.artifacts_arousals_movement_bin_mask(channelIdx,:,3)= results.artifacts_arousals_movement_bin_mask(channelIdx,:,3) .* originalMask_N4;
                        
                        % combine
                        results.artifacts_arousals_movement_bin_mask_fusion(channelIdx,:) = (sum(results.artifacts_arousals_movement_bin_mask,3)>0);
                        results.artifacts_arousals_movement_report{channelIdx,1} = tabulate(sum(results.artifacts_arousals_movement_bin_mask,3)>0);
                        results.artifacts_arousals_movement_report_N2{channelIdx,1} = tabulate(results.artifacts_arousals_movement_bin_mask(channelIdx,:,1));
                        results.artifacts_arousals_movement_report_N3{channelIdx,1} = tabulate(results.artifacts_arousals_movement_bin_mask(channelIdx,:,2));
                        results.artifacts_arousals_movement_report_N4{channelIdx,1} = tabulate(results.artifacts_arousals_movement_bin_mask(channelIdx,:,3));
                    else
                        rmsMedian_all = median(rmsSignal(:));rmsMedian_irq_all=iqr(rmsSignal(:));
                        results.artifacts_arousals_movement_bin_mask_fusion(channelIdx,:,1)=(rmsSignal > (rmsMedian_all + opts.artifacts.arousals_movement.threshold*rmsMedian_irq_all))|(rmsSignal < (rmsMedian_all - opts.artifacts.arousals_movement.threshold*rmsMedian_irq_all));
                        results.artifacts_arousals_movement_bin_mask_fusion(channelIdx,:,1) = f_stateDurationCheck(results.artifacts_arousals_movement_bin_mask_fusion(channelIdx,:,1),opts.artifacts.arousals_movement.over_threshold_duration*round(1000/data_filter.fsample));
                        results.artifacts_arousals_movement_report{channelIdx,1}= tabulate(results.artifacts_arousals_movement_bin_mask_fusion(channelIdx,:,1));
                    end
                end
                fprintf('artifacts detection: arousals movement.success!\n')
                clear data_filter
            catch
                fprintf('artifacts detection: arousals movement.fail!\n')
            end
        otherwise
            fprintf('Error, I do not know to do what! Try again!\n')
    end
end
end

