function [msk_update] = f_stateDurationCheck(msk,durationThreshold)
if size(msk,1)>size(msk,2), msk = msk'; end
msk_duration = diff(find([1,diff(msk),1]));
if size(durationThreshold,2) ~= 1
    pass_msk = (msk_duration > durationThreshold(1,1)) & (msk_duration < durationThreshold(1,2));
    %artifactsMsk_updated = repelem(pass_msk,artifactsMsk_duration);
    msk_update = msk.*repelem(pass_msk,msk_duration);
%     sum(msk(1682093:1682143))
%     imagesc(msk(1682093:1682143));
%     figure
%     imagesc(msk_update(1682093:1682143));
    
else
    pass_msk = msk_duration > durationThreshold;
    %artifactsMsk_updated = repelem(pass_msk,artifactsMsk_duration);
    msk_update = msk.*repelem(pass_msk,msk_duration);
end



