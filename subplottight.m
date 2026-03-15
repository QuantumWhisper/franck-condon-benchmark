function h = subplottight(n,m,i,varargin)
% Plot subplot more tightly with exact space for title and labels

[c,r] = ind2sub([m n], i);
widthper=0.13;
heightper=0.2;
bottomper=0.13;
if numel(varargin) > 0
    if ~isempty(varargin{1})
        widthper = varargin{1};        
    end
    if numel(varargin) > 1
        if ~isempty(varargin{2})
            heightper = varargin{2};
        end
        if numel(varargin) > 2
            if ~isempty(varargin{3})
                bottomper = varargin{3};
            end
        end
    end
end
ax = subplot('Position', [(c-1)/m+widthper*1/m, 1-(r)/n+(1+heightper)*bottomper*1/n, (1-2*widthper)*1/m, (1-2*bottomper)*1/n]);
if(nargout > 0)
    h = ax;
end
end