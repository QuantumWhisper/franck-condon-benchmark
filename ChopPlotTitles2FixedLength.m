function new_title = ChopPlotTitles2FixedLength(old_title, NumberOfCharsPerLine, varargin)
%%% @title A function to resize plot titles
%%% @author Shanglong Ning (sn538)
%%% @notice 
%{
%%% @param old_title string or cell
%%% @param NumberOfCharsPerLine number of characters for each line in title
%%% @param varargin
    varargin{1}, true - reshape title strings based on fixed lines
%%% @return new_title cell, title string
%%% @dev version v1 @ 2021/Aug/13
    old comments:
    This supportive function works with sn538 databse to reshape title into
    fixed length initially to cope with wordy comments
%}



new_title = {};
keep_lines = false; 
nvar = numel(varargin);
if nvar == 1
    keep_lines = varargin{1};
end

if iscell(old_title)
    t_string = strjoin(reshape(string(old_title),1,[]));
    if keep_lines
        index1s=find((cellfun(@length,old_title)>NumberOfCharsPerLine)==1);
        
        for nt=1:numel(old_title)
            
            if ~ismember(nt,index1s)
                tempind=numel(new_title);
                new_title{tempind+1}=old_title{nt};
            else
                currenttitle=regexp(old_title{nt}, sprintf('.{1,%d}', NumberOfCharsPerLine), 'match');
                tempind=numel(new_title);
                for cti=1:numel(currenttitle)
                    new_title{tempind+cti}=currenttitle{cti};
                end
            end
            
        end
    else
        new_title = regexp(t_string, sprintf('.{1,%d}',NumberOfCharsPerLine), 'match');%\\w
    end

    
else
    new_title = regexp(old_title, sprintf('.{1,%d}',NumberOfCharsPerLine), 'match');%\\w
end
end