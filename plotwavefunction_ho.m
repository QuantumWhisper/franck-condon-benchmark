fig = CMDA_createTagedFigure('harmonic oscillator');
fig.Position = [ 372   470   836   547];
nrows = 1;
ncols = 2;
cax1 = subplottight(nrows,ncols,1);
cax2 = subplottight(nrows,ncols,2);
x = linspace(-10,10,200);
line_offset = 1.115;
for n = 0:5
    psi = wavefunction_ho(x,n);
    line_ = plot(cax1,x,psi);
    line_.UserData.SweepNumber = n;
    hold(cax1,'on')

    line2_ = plot(cax2,x,psi.^2,'LineWidth',2);
    line2_.UserData.SweepNumber = n;
    hold(cax2,'on')

end
CMDA_offsetLines(cax1,line_offset)
CMDA_offsetLines(cax2,line_offset)
%% 0.16*x.^2
fig = CMDA_createTagedFigure('harmonic oscillator on hyperbola');
fig.Position = [ 372   470   836   547];
nrows = 1;
ncols = 2;
cax1 = subplottight(nrows,ncols,1);
cax2 = subplottight(nrows,ncols,2);
OSLineWidth = 4; %2
ParaLineWidth = 4;
xl = 5;
x = linspace(-xl,xl,200);
ye = 0.16*x.^2;
offset = 0.669;
current_offset = 0.5;
c1 = [170 28 11];
c2 = [62 40 174];
depth =6;
[plotTsColors,~] = colorGradient(c2,c1,depth);
for n = 0:3%5
    plotTsColor = plotTsColors(n+1,:);
    xr = sqrt(current_offset/0.16);
    x_cutoff = linspace(-xr,xr,200);

    psi = wavefunction_ho(x_cutoff,n);

    line_base = plot(cax2,x_cutoff,current_offset*ones(size(x_cutoff)),'Color','k',LineStyle='-.');
    hold(cax2,'on')

    line2_ = plot(cax2,x_cutoff,psi.^2 + current_offset,'LineWidth',OSLineWidth,'Color',plotTsColor);
    line2_.UserData.SweepNumber = n;
    hold(cax2,'on')
    current_offset = current_offset + offset;

end
plot(cax2,x,ye,'LineWidth',ParaLineWidth,'Color','k')
set(cax2,'XTick',[], 'YTick', [])

lambda = 7;
yoffset = 0.5*(lambda/7)^(1/30);
xl = 5;
x = linspace(-xl,xl,200);
ye = 0.16*x.^2;
offset = 0.669;
current_offset = 0.5;

for n = 0:5
    plotTsColor = plotTsColors(n+1,:);
    xr = sqrt(current_offset/0.16);
    x_cutoff = linspace(-xr,xr,200);

    psi = wavefunction_ho(x_cutoff,n);

    line_base = plot(cax2,lambda + x_cutoff,yoffset + current_offset*ones(size(x_cutoff)),'Color','k',LineStyle='-.');
    hold(cax2,'on')

    line2_ = plot(cax2,lambda + x_cutoff,yoffset + psi.^2 + current_offset,'LineWidth',OSLineWidth,'Color',plotTsColor);
    line2_.UserData.SweepNumber = n;
    hold(cax2,'on')
    current_offset = current_offset + offset;

end
plot(cax2,lambda+x,yoffset + ye,'LineWidth',ParaLineWidth,'Color','k')
set(cax2,'XTick',[], 'YTick', [])

fname = sprintf('demo-harmonic oscillator-lambda%g',lambda);
fdir = pwd; % Save to current directory
CMDA_savefig(cax2,fname,fdir,[],true)
%%

