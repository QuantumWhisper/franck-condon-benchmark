
xlists = [];
for ii = 1:numel(cursor_info)
    cursor0 = cursor_info(ii);
    pos = cursor0.Position;
    xlists = [xlists pos(2)];
end
%%
xlists = fittingResultsTable.Prominence;
xlists = fliplr(xlists);
%%
xlists = [2.6735e-7 1.8136e-07 1.7867e-07 1.1288e-07 9.1673e-08 6.6922e-08 5.8633e-08];

formula = 'FranckCondonElementRatio(n,lambda)';
coefficientscell = {'lambda'};

s = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[0],...
    'Upper',[10],...
    'Startpoint',[1]);

f = fittype(formula,'dependent',{'y'},'independent',{'n'},'coefficients',coefficientscell,'options', s);
x2fit = 1:numel(xlists);
y2fit = xlists/xlists(1);
if ~iscolumn(x2fit)
    x2fit = x2fit';
end
if ~iscolumn(y2fit)
    y2fit = y2fit';
end
[gfit,gof] = fit(x2fit,y2fit,f);

lambda = gfit.lambda;
figure
cax = subplottight(1,1,1);
plot(cax,x2fit,y2fit,'ok')
hold(cax,"on")
n = x2fit;
plot(cax,n,eval(formula),'r')

