function [Outputs,OutputsLinCheck] = CombineFitPublish
clear all
clc
warning off
load('Surveys.mat')
Outputs = [];
OutputsLinCheck = [];
%% Waterladder
x(:,1) = Surveys.PopDensity;
x(:,2) = Surveys.ToiletsRoad;
x(:,3) = Surveys.BinsEuc;
x(:,4) = Surveys.NalasRoad;
x(:,5) =  Surveys.ParksEuc;
x(:,6) =  Surveys.PlayGroundRoad;
x(:,7) =  Surveys.OpenSpaceRoads;
x(:,8) =  Surveys.ForestEuc;
x(:,9) =  Surveys.WateRoads;
x(:,10) =  Surveys.TownCentreEuc;

Trans.logs = [1,4];
Trans.Quads = [2,3;4,5;6,7;8,NaN;9,10;11,12;13,14;15,16;17,18;19,NaN];
Trans.transY = 'None';
modelFun = @(b,xD) b(1) +...
    (b(2).*xD(:,1)) + (b(3).*(xD(:,1).^2)) + ...
    (b(4).*xD(:,2)) + (b(5).*(xD(:,2).^2)) + ....
    (b(6).*xD(:,3)) + (b(7).*(xD(:,3).^2)) + ....
    (b(8).*xD(:,4)) + ....
    (b(9).*xD(:,5)) + (b(10).*(xD(:,5).^2)) + ....
    (b(11).*xD(:,6)) + (b(12).*(xD(:,6).^2)) + ....
    (b(13).*xD(:,7)) + (b(14).*(xD(:,7).^2)) + ....
    (b(15).*xD(:,8)) + (b(16).*(xD(:,8).^2)) + ....
    (b(17).*xD(:,9)) + (b(18).*(xD(:,9).^2)) + ....
    (b(19).*xD(:,10));
Trans.prior = ones(1,19).*(1/19);
Outputs = CreateOutputs(Outputs,1,x,modelFun,List,Surveys,Values,Trans);
%Linear Check
Trans.logs = [1,2,3,4,5,6,7,8,9,10];
Trans.Quads = [2,NaN;3,NaN;5,NaN;6,NaN;7,NaN;8,NaN;9,NaN;10,NaN;11,NaN];
Trans.transY = 'None';
modelFun = @(b,xD) b(1) +...
    (b(2).*xD(:,1)) + ...
    (b(3).*xD(:,2)) + ....
    (b(4).*xD(:,3)) + ....
    (b(5).*xD(:,4)) + ....
    (b(6).*xD(:,5)) + ....
    (b(7).*xD(:,6))+ ....
    (b(8).*xD(:,7)) + ....
    (b(9).*xD(:,8)) + ....
    (b(10).*xD(:,9)) + ....
    (b(11).*xD(:,10));
Trans.prior = ones(1,11).*(1/11);
OutputsLinCheck = CreateOutputs(OutputsLinCheck,1,x,modelFun,List,Surveys,Values,Trans);
clearvars -except Outputs List Surveys Values OutputsLinCheck 
%%
% Other indices are done identically
end

function Outputs = CreateOutputs(Outputs,nr,xIn,modelFun,List,Surveys,Values,Trans)
y = Surveys.(genvarname([char(List.Dependent(nr))]));
outs = regfunc(xIn,y,modelFun,Trans);
outs = SensiFunc(outs,modelFun,xIn,Trans);
Outputs.(genvarname([char(List.Dependent(nr))])).Orginal = outs;
y = Values.(genvarname([char(List.Dependent(nr))])).marginal;
outs = regfunc(xIn,y,modelFun,Trans);
outs = SensiFunc(outs,modelFun,xIn,Trans);
Outputs.(genvarname([char(List.Dependent(nr))])).Marginal = outs;
display('finished')
end

function outs = regfunc(xIn,y,modelFun,Trans)
xD = logfunctions(xIn,Trans.logs);
[y,Trans] = yTransfer(y,Trans);
Options = statset('FunValCheck','off','Display','off','MaxIter',200,...
    'TolFun',1.0000e-4, 'TolX',1.0e-4, 'DerivStep', 6.0555e-06,'Robust',...
    'off', 'WgtFun', 'bisquare');
[beta,residuals,J,~,~]  = nlinfit(xD,y,modelFun,Trans.prior,'Options',Options);
ci = nlparci(beta,residuals,'Jacobian',J);
t = tinv(1-(0.05/2),length(y)-length(beta));
sebeta= (ci(:,2)-ci(:,1)) ./ (2*t);  % Standard Error
z = beta./sebeta';
pvalu = 2*(1 - normcdf(abs(z)));
%Rsquare
sstot = nansum((y-nanmean(y)).^2);
ssres = residuals.^2;
outs.Rsquare(1) = max(0,(1- (nansum(ssres)/sstot))); % Rsquare
if isnan(outs.Rsquare(1,1))==1 || isinf(outs.Rsquare(1,1))==1
    outs.Rsquare(1,1) = 0;
end
Beta = dataset(reshape(beta,length(beta),1),'Varnames',{'Beta'});
Beta.SeBeta = reshape(sebeta,length(sebeta),1);
Beta.PValue =reshape(pvalu,length(pvalu),1);
outs.Beta = Beta;
outs.Values = y-residuals;
outs = CurveType(outs,Trans,Beta,xD);
% Ranges overview
ranges_txt = {'Range';'Min';'Max';'Mean';'STD'};
Ranges(1,1) = range(y);
Ranges(2,1) = min(y);
Ranges(3,1) = max(y);
Ranges(4,1) = nanmean(y);
Ranges(5,1) = nanstd(y);
outs.Ranges = dataset(Ranges,'ObsNames', ranges_txt,'Varnames',{'Observed'});
clear Ranges
Ranges(1,1) = range(outs.Values);
Ranges(2,1) = min(outs.Values);
Ranges(3,1) = max(outs.Values);
Ranges(4,1) = nanmean(outs.Values);
Ranges(5,1) = nanstd(outs.Values);
outs.Ranges.ThisModel = Ranges;
end

function outs = SensiFunc(outs,modelFun,xIn,Trans)
% sensitivity
save('all')
xOrg = xIn;
xD = logfunctions(xIn,Trans.logs);
b = outs.Beta.Beta;
RefSens = nanmean(modelFun(b,xD));
clear xD
Sens = zeros(size(xIn,2),1);
for i = 1:1:size(xIn,2)
    xTrans = xOrg;
    xTrans(:,i) = xTrans(:,i)*1.25;
    xD = logfunctions(xTrans,Trans.logs);
    tester = modelFun(b,xD);
    value = nanmean(tester);
    Sens(i,1) = abs(value./RefSens);
    if Sens(i,1) < 1
        Sens(i,1) = 1/Sens(i,1);
    end
    clear xTrans xD
end
Sens(:,1) = Sens(:,1) -1;
Sens = Sens./sum(Sens);
outs.Sens = dataset(Sens,'Varnames','Sensitivity');
outs.Sens.Direction = outs.TextDirection;
outs = rmfield(outs,'TextDirection');
end

function outs = CurveType(outs,Trans,Beta,xD)
for i = 1:size(Trans.Quads,1)
    test = find((isnan(Trans.Quads(i,:))==1));
    if isempty(test) == 1
        xRange = max(xD(:,i));
        save('all')
        factor = abs(Beta.Beta(Trans.Quads(i,1))./Beta.Beta(Trans.Quads(i,2)));
        Txt = 'Concave(firstUp)';
        Txt2 = '_but_Decrease';
        if  Beta.Beta(Trans.Quads(i,2)) > 0 && Beta.Beta(Trans.Quads(i,1)) < 0
            Txt = 'Convex(firstDown)';
            if factor < xRange;
                Txt2= '_but_Increase';
            else
                Txt2= '_remaining_Decrease';
                if factor < 2*(xRange);
                    Txt2= [Txt2,'_True_Humb'];
                end
            end
        elseif Beta.Beta(Trans.Quads(i,2)) < 0 && Beta.Beta(Trans.Quads(i,1)) > 0
            if factor > xRange;
                Txt2= '_remaining_Increase';
                if factor < 2*(xRange);
                    Txt2= [Txt2,'_True_Humb'];
                end
            end
        else
            Txt = 'Full';
            Txt2 = '_Decreasing';
            if Beta.Beta(Trans.Quads(i,2)) > 0
                Txt2 = '_Increasing';
            end
        end
        Txt = [Txt,Txt2]; %#ok<*AGROW>
    else
        Txt = 'Increasing';
        if Beta.Beta(Trans.Quads(i,1)) < 0
            Txt = 'Decreasing';
        end
    end
    outs.TextDirection(i,1) = {Txt};
end
end

function xD = logfunctions(x,logs)
xD = x;
for i = logs
    xD(:,i) = log10(x(:,i)+1);
end
end

function [Yval,Trans] = yTransfer(Yval,Trans)
if isfield(Trans,'transY') == 1
    if strcmp('Log',Trans.transY) % reverse = 10^b
    Yval = log10(Yval+1);
    display('Log y transfer')
    elseif strcmp('SquareRoot',Trans.transY)% reverse = b^2
     Yval = sqrt(Yval);
    display('Sqrt y transfer')
    elseif strcmp('Arcsin',Trans.transY) % reverse = (sin(b))^2
    Yval = asin(sqrt(Yval));
    display('Arcsn y transfer')  
    end
end
end
