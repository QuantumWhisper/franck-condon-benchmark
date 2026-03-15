function CMDA_savefig(varargin)
%%% @title A function to save figures or axes with scitific printing
%%% quality
%%% @author Shanglong Ning (sn538)
%%% @notice
%{
%%% @param fig, figure,axes,figure cell, axes cell, handles
%%% @param varargin
varargin{1}: figure handle, default gcf
varargin{2}: figure name, string, default untitled or figure's tag or
figure's name
varargin{3}: save directory, -d, folder path, default pwd
varargin{4}: save everything as one pdf file, true/false, default false
varargin{5}: save figure .fig, true/false, default false
%%% @return 
%%% @dev version v1 @ 2021/Oct./19
example: CMDA_savefig(gcf)
@dev version v2 @ 2022/Apr./27
updates on plot pdf as default for publication qualities
updates on saving figures now all have same default size unless specified
otherwise
% change ext = '.pdf';%'.pdf'; %'.png'; from .eps tp .pdf
%}

tic;
fprintf('Start saving ...\n')
fig = gcf;
if isa(fig,'matlab.ui.Figure')
    fname = sprintf('%s',strrep(erase(fig.Name,'TagedFigure:'),'/','-'));
elseif isa(fig,'matlab.graphics.axis.Axes')
    fname = fig.Tag;
else
    fname = 'untitled';
end
folderpath = pwd;
ext = '.pdf';%'.pdf'; %'.png';
resolution = 600; % Specifying the resolution has no effect when the ContentType is 'vector'.

PDFCollection = false;
isSaveFig = false;

defaultSize = false;

if numel(varargin) > 0
    if ~isempty(varargin{1})
        fig = varargin{1};
    end 
    if numel(varargin) > 1
        if ~isempty(varargin{2})
            fname = varargin{2};
        end
    end
    if numel(varargin) > 2
        if exist(varargin{3},'dir') == 7
            folderpath = varargin{3};
        end
    end
    if numel(varargin) > 3
        if ~isempty(varargin{4})
            PDFCollection = varargin{4};
        end
    end
    if numel(varargin) > 4
        if ~isempty(varargin{5})
            isSaveFig = varargin{5};
        end
    end
    if numel(varargin) > 5
        if ~isempty(varargin{6})
            ext = varargin{6};
        end
    end
end

if isa(fig,'matlab.graphics.axis.Axes')
    if defaultSize
        h = fig.Parent;
        h.Position = get(0, 'Screensize');
        hunits = get(h,'Units');
        punits = get(h,'PaperUnits');
        psize = get(h,'PaperSize');
        set(h,'PaperUnits','centimeters');
        set(h,'Units','centimeters');
        pos = get(h,'Position');
        set(h,'PaperSize', [pos(3) pos(4)]);
        set(h,'PaperPositionMode', 'manual');
        set(h,'PaperPosition',[0 0 pos(3) pos(4)]);
    end
    if strcmp(fig.Tag,fname)
        fname = sprintf('%s - %s',fig.Tag, fname);
    end    
    fpath = [fullfile(folderpath,fname) ext];
    count = 1;
    while (exist(fpath,'file')==2)
        fname0 = sprintf('%s - %g',fname,count);
        fpath = [fullfile(folderpath,fname0) ext];
        count = count + 1;
    end
    exportgraphics(fig,fpath,'Resolution',resolution,'ContentType','auto','BackgroundColor','none')
    if isSaveFig
        savefig(fig.Parent,[erase(fpath,ext) '.fig'],'compact');
    end
    if defaultSize
        set(h,'PaperSize',psize);% restoring this seems to prevent plotting pngs from clearing existing images!
        set(h,'Units',hunits);% restore old units
        set(h,'PaperUnits',punits);% restore old units
    end
elseif isa(fig,'matlab.ui.Figure')
    fpath = [fullfile(folderpath,fname) ext];
    count = 1;
    while (exist(fpath,'file')==2)
        fname0 = sprintf('%s - %g',fname,count);
        fpath = [fullfile(folderpath,fname0) ext];
        count = count + 1;
    end
    exportgraphics(fig,fpath,'Resolution',resolution,'ContentType','auto','BackgroundColor','none')
    if isSaveFig
        savefig(fig,[erase(fpath,ext) '.fig'],'compact');
    end
elseif iscell(fig)
    nfig = numel(fig);
    if ~PDFCollection
        for ii = 1:nfig
            CMDA_savefig(fig{ii},fname,folderpath,PDFCollection,isSaveFig);
        end
    else
        pnf = nfig;
        tolpdf = cell([1 pnf]);
        for pn = 1:pnf
            if isa(fig{pn},'matlab.ui.Figure')
                fig{pn}.Units = 'inches';
                figwidth = fig{pn}.Position(3);
                figheight = fig{pn}.Position(4);
                fig{pn}.PaperUnits = 'inches';
                fig{pn}.PaperSize = [ceil(1.1*figwidth) ceil(1.1*figheight)];
                filecount = 1;
                savefilename = sprintf('%s-%g.pdf',fig{pn}.Name,filecount);
                savefilename = strrep(savefilename,'/','_');
                singlepath1 = fullfile(folderpath,savefilename);
                while (exist(singlepath1,'file')==2)
                    filecount = filecount+1;
                    savefilename = sprintf('%s-%g.pdf',fig{pn}.Name,filecount);
                    savefilename = strrep(savefilename,'/','_');
                    singlepath1 = fullfile(folderpath,savefilename);
                end
                tolpdf{pn} = singlepath1;
                if pn == 1
                    foldername = strrep(fig{pn}.Name,'/','_');
                end
            end
        end
        parfor newpn = 1:pnf
            if isa(fig{newpn},'matlab.ui.Figure')
                print(fig{newpn},tolpdf{newpn},'-painters','-dpdf','-r600','-noui')
            end            
        end        
        sumcount = 1;
        SummaryFileName = fullfile(folderpath,sprintf('0-%s-%g.pdf',foldername,sumcount));
        while (exist(SummaryFileName,'file')==2)
            sumcount=sumcount+1;
            SummaryFileName=fullfile(folderpath,sprintf('0-%s-%g.pdf',foldername,sumcount));
        end
        if numel(tolpdf) > 2
            if isa(fig{1},'matlab.ui.Figure')
                try
                    append_pdfs(SummaryFileName,tolpdf{:})
                catch
                    warning('append_pdfs not available. Individual PDFs saved but not merged.')
                end
            end            
        end        
    end
end

rtime = toc;
fprintf('%s\nsaved @%s\nRunning time: %g s.\n',fname,folderpath,rtime);
end