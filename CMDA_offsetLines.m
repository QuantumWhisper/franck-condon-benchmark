function CMDA_offsetLines(handle,varargin)

UserDataFN = 'CMDA_offsetLines';
ORDERNAME = 'Temperature';%'Temperature';%'SweepNumber';
if isa(handle,'matlab.ui.Figure')
    %{
    caxall = findobj(handle,'Type','axes');
    for ii = 1:numel(caxall)
        cax = caxall(ii);
        CMDA_offsetLines(cax)
    end
    %}
    cax = handle.CurrentAxes;
    CMDA_offsetLines(cax,varargin{:})
    
elseif isa(handle,'matlab.graphics.axis.Axes')
    if ~isempty(handle)
        if numel(varargin) == 2
            if ~isempty(varargin{2})
                ORDERNAME = varargin{2};
            end
        end
        cax = handle;
        [lines_ , mlines_, clines_] = CMDA_getlines(cax);
        nlines = numel(lines_);
        nmlines = numel(mlines_);
        nclines = numel(clines_);
        % data lines
        xdata_lists = cell([1 nlines]);
        ydata_lists = cell([1 nlines]);
        userdata_lists = cell([1 nlines]);
        linespec_lists = cell([1 nlines]);
        linespec_names = {'LineStyle','LineWidth','Marker','MarkerEdgeColor','MarkerFaceColor','MarkerSize','Color'};
        Plotdata = cell([1 nlines]);
        SweepNumberList = zeros(1,nlines);
        for ll = 1:nlines
            line_ = lines_(ll);
            temp_x = line_.XData;
            temp_y = line_.YData;
            xdata_lists{ll} = temp_x;
            ydata_lists{ll} = temp_y;
            userdata_lists{ll} = line_.UserData;
            linespec_lists{ll} = get(line_,linespec_names);
            if ~iscolumn(temp_x)
                temp_x = temp_x';
            end
            if ~iscolumn(temp_y)
                temp_y = temp_y';
            end
            Plotdata{ll}(:,1) = temp_x;
            Plotdata{ll}(:,2) = temp_y;
            if isstruct(line_.UserData)
                if isfield(line_.UserData,ORDERNAME)%'SweepNumber')
                    %SweepNumberList = [SweepNumberList line_.UserData.SweepNumber];
                    SweepNumberList(ll) = line_.UserData.(ORDERNAME);
                end
            end
        end

        current_offset = 0;
        try
            line_offset = cax.UserData.offsetLines.line_offset + 0.01;
            cax.UserData.offsetLines.line_offset = line_offset;
        catch
            line_offset = 0.025;
            cax.UserData.offsetLines.line_offset = line_offset;
        end
        if numel(varargin) == 1
            if ~isempty(varargin{1})
                line_offset = varargin{1};
            end
        elseif numel(varargin) == 2
            if ~isempty(varargin{1})
                line_offset = varargin{1};
            end
        end
        offset = (max(cax.YLim)-min(cax.YLim))*line_offset;

        if ~isempty(cax.XLabel)
            caxXLabel = cax.XLabel.String;
        end
        if ~isempty(cax.YLabel)
            caxYLabel = cax.YLabel.String;
        end

        tagAxes = 'CMDA_offsetLinesAxis';
        try
            prevAxis = findobj('Tag',tagAxes);
            prevFig = prevAxis.Parent;
            axisPos = prevAxis.Position;
            figPos = prevFig.Position;
        catch
            axisPos = [];
            figPos = [];
        end

        tagFigure = 'CMDA_offsetLinesFigure';
        fighs = CMDA_createTagedFigure(sprintf(': %s',tagFigure));
        fighs.UserData.(UserDataFN) = handle;
        axhs = subplottight(1,1,1);
        
        axhs.Tag = tagAxes;
        
        legend(axhs,'AutoUpdate','off')
        if ~isempty(figPos)
            fighs.Position = figPos;
        end
        if ~isempty(axisPos)
            %axhs.Position = axisPos;
        end

        hold(axhs,'on');
        % try animate in order of sweep numer
        [~,SwpIndex] = sort(SweepNumberList); % ascending e.g 1 2 3 4
        if isempty(SwpIndex)
            SwpIndex = 1:nlines;
        end

        im = {};
        line_color = rand(1,3);
        hs = [];
        LegendStringCell = cell([1 numel(SwpIndex)]);

        for jj = 1:numel(SwpIndex)
            ii = SwpIndex(jj);
            swp = SweepNumberList(ii);
            xdata = xdata_lists{ii};
            ydata = ydata_lists{ii};
            Specs = linespec_lists{ii};
            userdata = userdata_lists{ii};

            if swp ~= 0
                legstr = sprintf('sweep = %g',swp); 
                if strcmp(ORDERNAME,'Temperature')
                    legstr = sprintf('T \\approx %.0f K',swp); 
                end
            else
                legstr = sprintf('Line number = %g',swp); 
            end

            h = plot(axhs,xdata,ydata + current_offset);

            hs = [hs, h];
            set(h,linespec_names,Specs)
            h.UserData = userdata;
            hold(axhs,'on')
            current_offset = current_offset + offset;
            LegendStringCell{jj} = legstr;
        end
        hold(axhs,'off');
        cax.UserData.offsetLines.line_offset_abs = offset;
        axhs.UserData.offsetLines.line_offset_abs = offset;
        legend(axhs,LegendStringCell,'Location','bestoutside');
        if ~isempty(caxXLabel)
            xlabel(axhs,cax.XLabel.String);
        end
        if ~isempty(caxYLabel)
            ylabel(axhs,cax.YLabel.String);
        end

    end
elseif iscell(handle)
    nfig = numel(handle);
    for ii = 1:nfig
        CMDA_offsetLines(handle{ii},varargin);
    end
else
    warning('Please check your inputs for CMDA_offsetLines: figure, axes or cell!')
    return
end

end