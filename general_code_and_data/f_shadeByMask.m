function [x,y] = f_shadeByMask(mask,highlight,str,fin)
areaDetect = (mask==highlight);
changePoint = find(diff(areaDetect));
if areaDetect(1)==1
    changePoint=[1,changePoint];
end
if areaDetect(end)==1
    changePoint=[changePoint,length(areaDetect)];
end

if ~isempty(changePoint)
    if mod(length(changePoint),2)==0
        temp=reshape(changePoint,2,[])';
        for i=1:size(temp,1)
            x(i,1)= temp(i,1);x(i,2)= temp(i,1);
            x(i,3)= temp(i,2);x(i,4)= temp(i,2);
        end
    end
else
    i=1;
    x(i,1)= 1;x(i,2)= 1;
    x(i,3)= length(mask);x(i,4)= length(mask);
end
y = repmat([str,fin,fin,str],[size(x,1),1]);
end
