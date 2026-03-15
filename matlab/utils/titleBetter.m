function tt = titleBetter(cax,total_title,varargin)

fig = cax.Parent;

axes_ch = findobj(fig,'Type','axes');
MaxNumberOfSubPlots = numel(axes_ch);
if MaxNumberOfSubPlots == 1
    rown = 1; % the same Y position means the same row
    coln = 1; % the same X position means the same column
    count = 1;
else
    pos = cell2mat(get(axes_ch,'position'));
    rown = numel(unique(pos(:,2))); % the same Y position means the same row
    coln = numel(unique(pos(:,1))); % the same X position means the same column
    count = 1;
end

cax.TitleFontSizeMultiplier=0.8;
NoCharPerLine=round(1.5*fig.Position(3)/coln/(cax.FontSize*cax.TitleFontSizeMultiplier));
retitle=ChopPlotTitles2FixedLength(total_title,NoCharPerLine);
%retitle = strrep(retitle,'\','\\'); % cope with latex problem
tt = title(cax,retitle,varargin{:});

end