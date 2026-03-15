function [dlines_ , mlines_, clines_] = CMDA_getlines(cax,varargin)
% Return data lines, lines with markers, and ConstantLine objects.
%   [dlines, mlines, clines] = CMDA_getlines(cax)
%   [...] = CMDA_getlines(cax, flipd, flipm, flipc)
%   [...] = CMDA_getlines(..., 'VisibleOnly', true)

visibleOnly = false;
nvIdx = [];
for k = 1:numel(varargin)
    if (ischar(varargin{k}) || isstring(varargin{k})) ...
            && strcmpi(varargin{k}, 'VisibleOnly') && k < numel(varargin)
        visibleOnly = logical(varargin{k+1});
        nvIdx = [k, k+1];
    end
end
varargin(nvIdx) = [];

flipdlines = false;
flipmlines = false;
flipclines = false;
if numel(varargin) > 0
    if ~isempty(varargin{1})
        flipdlines = varargin{1};
    end
    if numel(varargin) > 1
        if ~isempty(varargin{2})
            flipmlines = varargin{2};
        end
        if numel(varargin) > 2
            if ~isempty(varargin{3})
                flipclines = varargin{3};
            end
        end
    end
end
lines_ = cax.Children;

dlines_ = [];
mlines_ = [];
clines_ = [];

for ls = 1:numel(lines_)
    line_ = lines_(ls);
    if visibleOnly && strcmp(line_.Visible, 'off')
        continue;
    end
    if strcmp(line_.Type,'line')
        if contains(line_.Tag,'Markers')
            mlines_ = [mlines_ line_];
        else
            dlines_ = [dlines_ line_];
        end
    elseif strcmp(line_.Type,'constantline')
        clines_ = [clines_ line_];
    end
end

if flipdlines
    dlines_ = fliplr(dlines_);
end
if flipmlines
    mlines_ = fliplr(mlines_);
end
if flipclines
    clines_ = fliplr(clines_);
end


end