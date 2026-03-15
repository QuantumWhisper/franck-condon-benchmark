function [fig,tag,pre_ran] = CMDA_createTagedFigure(varargin)
%%% @title A function to create a figure with tag tag, and close the other
%%% figures with the same tag. Now supports standard figure() syntax.
%%% @author Shanglong Ning (sn538)
%%% @notice 
%{
%%% Usage examples:
%%%   fig = CMDA_createTagedFigure('myTag')  % Original style
%%%   fig = CMDA_createTagedFigure('myTag', true, false)  % Original style with options
%%%   fig = CMDA_createTagedFigure('Position', [100, 100, 800, 600])  % figure() style
%%%   fig = CMDA_createTagedFigure('Position', [100, 100, 800, 600], 'Tag', 'myTag')  % Mixed style
%%%   fig = CMDA_createTagedFigure()  % No arguments, auto-generated tag
%%% 
%%% Original parameters (when first argument is not a property name):
%%% @param tag (string) - figure tag
%%% @param varargin{1} (logical) - attach a 64 random string to given tag with prefix $r$
%%% @param varargin{2} (logical) - deep removal 
%%% @param varargin{3} (string) - delete figures if contain this
%%%
%%% Standard figure properties (name-value pairs):
%%% Any property supported by figure() function
%%%
%%% @return fig figure handle
%%% @return tag string, final tag used
%%% @return pre_ran string, prefix for random string
%%% @dev version v1 @ 2021/July/29
%%%   @dev version v1 @ 2021/Aug/13 add randomness to tag; ability to turn off closure function
%%%   @dev version v1 @ 2021/Oct/15 return pre_ran
%%%   @dev version v2 @ 2025/May/24 add compatibility with figure() syntax
%}

% Initialize default values
tag_input = '';
random_string = '';
pre_ran = '$r$';
random_string_on = false;
closePreviousFigure = true;
del_str_exception = '';
figure_props = {};

% Determine if we're using original style or figure() style
if nargin == 0
    % No arguments - generate automatic tag
    tag_input = ['AutoTag_' datestr(now, 'yyyymmdd_HHMMSS_FFF')];
    
elseif nargin >= 1 && ischar(varargin{1}) && ~any(strcmpi(varargin{1}, {'Position', 'Name', 'Color', 'MenuBar', 'ToolBar', 'NumberTitle', 'Visible', 'Units', 'OuterPosition', 'PaperPosition', 'PaperSize', 'PaperType', 'PaperUnits', 'Resize', 'WindowStyle', 'CloseRequestFcn', 'CreateFcn', 'DeleteFcn', 'ButtonDownFcn', 'KeyPressFcn', 'KeyReleaseFcn', 'WindowButtonDownFcn', 'WindowButtonMotionFcn', 'WindowButtonUpFcn', 'WindowKeyPressFcn', 'WindowKeyReleaseFcn', 'WindowScrollWheelFcn'}))
    % Original style - first argument is tag, not a property name
    tag_input = varargin{1};
    
    % Process remaining arguments in original style
    remaining_args = varargin(2:end);
    num_v = numel(remaining_args);
    
    if num_v >= 1
        random_string_on = remaining_args{1};
    end
    if num_v >= 2
        closePreviousFigure = remaining_args{2};
    end
    if num_v >= 3
        del_str_exception = remaining_args{3};
    end
    
else
    % Figure style - name-value pairs or just properties
    % Look for 'Tag' property to use as our tag
    tag_idx = find(strcmpi(varargin(1:2:end), 'Tag'), 1);
    if ~isempty(tag_idx)
        tag_input = varargin{2*tag_idx};
        % Remove Tag from the properties we'll pass to figure
        prop_indices = 1:nargin;
        prop_indices([2*tag_idx-1, 2*tag_idx]) = [];
        figure_props = varargin(prop_indices);
    else
        % No Tag specified, generate automatic tag
        tag_input = ['AutoTag_' datestr(now, 'yyyymmdd_HHMMSS_FFF')];
        figure_props = varargin;
    end
end

% Handle random string generation
if random_string_on
    random_string = char(randi([33 126],1,64)); % size 1 64
    random_string = [pre_ran random_string];
end

% Create final tag
tag = ['TagedFigure' tag_input random_string];

% Close previous figures with same tag
if closePreviousFigure
    fighs = findobj('Tag',tag);
    for ii = 1:numel(fighs)
        figToClose = fighs(ii);
        close(figToClose);
        delete(figToClose);
        clear figToClose;
    end
end

% Handle deep removal with exception string
if closePreviousFigure && (~isempty(del_str_exception))
    fighs2 = findobj('Type','figure');
    for ii = 1:numel(fighs2)
        figToClose = fighs2(ii);
        if contains(figToClose.Tag,'TagedFigure')
            if ~contains(figToClose.Tag,del_str_exception)
                close(figToClose);
                delete(figToClose);
                clear figToClose;
            end
        end
    end
end

% Create figure with properties
if isempty(figure_props)
    fig = figure;
else
    fig = figure(figure_props{:});
end

% Set custom properties
fig.Tag = tag;
% fig.WindowKeyPressFcn = @CMDA_WindowKeyPressFcn_tools; % Removed for standalone use
%fig.WindowButtonMotionFcn = @CMDA_WindowButtonDownFcn_cuttingProfile_General;
%fig.WindowButtonDownFcn = @CMDA_WindowButtonDownFcn_cuttingProfile_General;

% Only set Name if it wasn't already set by figure properties
if ~any(strcmpi(figure_props(1:2:end), 'Name'))
    fig.Name = tag;
end

end