obj = findall(gcf,'-property','FontSize');
if numel(obj)>1
    tmp = obj(cellfun(@(x)x>8,get(obj,'fontsize')));
    set(tmp,'FontSize',8)
else
    set(findall(gcf,'-property','FontSize'),'FontSize',8)
end
set(findall(gcf,'-property','FontType'),'FontType','Arial')
tmp = get(gca,'YColor');
if any(tmp > 0.9) % if yaxis is colored
    set(gca,'Xcolor','k','Zcolor','k')
else
    set(gca,'Xcolor','k','YColor','k','Zcolor','k')
end