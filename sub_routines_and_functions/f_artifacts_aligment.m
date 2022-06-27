function [results] = f_artifacts_aligment(opts,previousResults,EEG)
results = [];
for j=1:numel(opts.artifacts.amplitude.cfg.channel)
    channelIdx = find(strcmp(EEG.label,opts.artifacts.amplitude.cfg.channel));
    art_amplitude = previousResults.artifacts_amplitude_bin_mask(channelIdx,:);
    art_gradient = previousResults.artifacts_gradient_bin_mask_fusion(channelIdx,:);
    art_arousal_movement = previousResults.artifacts_arousals_movement_bin_mask_fusion(channelIdx,:);
    
    fullMask{channelIdx,:} = [art_amplitude;[art_gradient(1),art_gradient];art_arousal_movement];
    % padded
    combineMask(channelIdx,:) = (sum(fullMask{channelIdx,:},1)>0);
    artpadding = round(opts.pad_duration/round(1000/EEG.fsample));
    
    artbeg = find(diff([0 combineMask(channelIdx,:)])== 1);
    artend = find(diff([0 combineMask(channelIdx,:)])== -1);
    artbeg = artbeg - artpadding;
    artend = artend + artpadding;
    artbeg(artbeg<1) = 1;
    artend(artend>EEG.sampleinfo(2)) = EEG.sampleinfo(2);
    temp = combineMask(channelIdx,:);
    if ~isempty(artbeg)
        for artlop=1:length(artbeg)
            temp(artbeg(artlop):artend(artlop)) = 1;
        end
    end
    % < xs mark as artifacts
    temp_withoutArt = ~temp;
    noartifactsMsk_duration = diff(find([1,diff(temp_withoutArt),1]));
    tooShort_msk = noartifactsMsk_duration < round(opts.minWithoutArt/round(1000/EEG.fsample));
    tooShort_msk_fulldatalen=all([temp_withoutArt;repelem(tooShort_msk,noartifactsMsk_duration)]);
    holeFill_artifacts = (temp+tooShort_msk_fulldatalen)>0;
    %results{j,1} = (~holeFill_artifacts) .* bin_msk_NREM;
    results{j,1} = (~holeFill_artifacts);
end
end

