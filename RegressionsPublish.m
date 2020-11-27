clear all
clc
warning off
load('Surveys.mat')
maxDepen = length(List.Dependent);
maxIndepen = length(List.Independent) ;
for yDepen = 1:maxDepen
    clc
    str = sprintf('      Running dependent = %d %s ',yDepen, char(List.Dependent(yDepen)));
    disp(str)
    y = Surveys.(genvarname([char(List.Dependent(yDepen))]));
    %% Create marginal values
    ds = dataset(y);
    [~,~,~,wij] = Morans(Surveys.Longitude_Residence,Surveys.Lattitude_Residence,1,y,1);
    Auto = zeros(1280,1);
    for s = 1:1280
        Autot = 0;
        for t = 1:1280
            if s ~= t && (isnan(y(t)))~=1
                Autot = Autot + (wij(s,t).*y(t));
            end
        end
        Auto(s,1) = Autot./nansum(wij(s,:)); %#ok<*AGROW>
    end
    ds.Auto = Auto;
    [~,outsAuto,statsAuto] = anovan(ds.y,{ds.Auto},'sstype',1,...
        'model','linear','continuous',[1],'display', 'off',...
        'varnames', {'Autocorrelation'});
    statsAuto.resid = reshape(statsAuto.resid,length(statsAuto.resid),1);
    ds.InBetween = statsAuto.resid;
    ds.Sites = Surveys.(genvarname([char(List.CoVar(1))]));
    ds.Age = Surveys.(genvarname([char(List.CoVar(2))]));
    ds.Embedding = ones(length(y),1);
    ds.Income = ones(length(y),1);
    ds.Education = ones(length(y),1);
    if yDepen ~= 5
        ds.Income = log10((Surveys.(genvarname([char(List.CoVar(4))])))+1);
    end
    if yDepen ~= 6
        ds.Education = Surveys.(genvarname([char(List.CoVar(3))]));
    end
    if yDepen ~= 4
        ds.Embedding =  (Surveys.(genvarname([char(List.CoVar(5))])));
    end
    [~,outsModel,statsModel] = anovan(ds.InBetween,{ds.Age,ds.Income,ds.Education},'sstype',3,...
        'model','linear','continuous',[1,2,3],'display', 'off',...
        'varnames', {'Age','Income','Education'}); % STEP 1a
    outs = outsModel(1,1:7);
    outs(2,1:7) = outsAuto(2,1:7);
    outs(3:7,1:7) = outsModel(2:6,1:7);
    outs(1,8) = {'Direction'};
    outs(3,8) = {statsModel.coeffs(length(statsModel.coeffs)-2)};
    outs(4,8) = {statsModel.coeffs(length(statsModel.coeffs)-1)};
    outs(5,8) = {statsModel.coeffs(length(statsModel.coeffs))};
    statsModel.resid = reshape(statsModel.resid,length(statsModel.resid),1);
    %% clear Autocorrelations and reshape residuals
    yRes = statsModel.resid;%unexplained variation by the model above
    %Scale y
    if min(yRes) <min(y)
        yRes = yRes + (min(y)-(min(yRes)));
    end
    yRes =  ((yRes-min(yRes))./(range(yRes./range(y))))+min(yRes);
    CoVar.(genvarname([char(List.Dependent(yDepen))])) = outs;
    %% regressions
    RegressionOuts = dataset({'Dummy';'Dummy'},'Varnames',char('Independent'));
    RegressionOuts.Function = {'Dummy';'Dummy'};
    RegressionOuts.RSquare =  [NaN;NaN];
    RegressionOuts.AIC =  [NaN;NaN];
    RegressionOuts.Constant = [NaN;NaN];
    RegressionOuts.B2 = [NaN;NaN];
    RegressionOuts.B3 = [NaN;NaN];
    RegressionOuts.P_Constant = [NaN;NaN];
    RegressionOuts.P_B2 = [NaN;NaN];
    RegressionOuts.P_B3 = [NaN;NaN];
    for xIndepen = 20:maxIndepen
        display(List.Independent(xIndepen))
        x = Surveys.(genvarname([char(List.Independent(xIndepen))]));
        for func = 1:1:8
            [Outs,text] = curvefitPublish(func,x,yRes);  % STEP 1b
            outfac = ((xIndepen-1).*10)+(func+1);
            RegressionOuts.Independent(outfac,1) = {List.Independent(xIndepen)};
            RegressionOuts.Function(outfac,1) = {text};
            RegressionOuts.RSquare(outfac,1) = Outs.Rsquare(1);
            RegressionOuts.AIC (outfac,1) = Outs.AIC;
            %RegressionOuts.P_Rsquare(outfac,1) = Outs.Rsquare(2);
            RegressionOuts.Constant(outfac,1) = Outs.Beta(1);
            RegressionOuts.B2(outfac,1) = Outs.Beta(2);
            RegressionOuts.B3(outfac,1) = Outs.Beta(3);
            RegressionOuts.P_Constant(outfac,1) = Outs.Beta(7);
            RegressionOuts.P_B2(outfac,1) =Outs.Beta(8);
            RegressionOuts.P_B3(outfac,1) = Outs.Beta(9);
        end
        %Single Variable Log Lin check for factor significance and distance
        %choice
        [~,outsCheck,statsCheck] = anovan(yRes,{log10(x+1)},'sstype',1,...
            'model','linear','continuous',[1],'display', 'off',...
            'varnames', {'CrossCheck'});
        facout = ((xIndepen-1).*10)+1;
        RegressionOuts.Independent(facout,1) = {List.Independent(xIndepen)};
        RegressionOuts.Function(facout,1) = {'Crosscheck'};
        RegressionOuts.RSquare(facout,1) =  cell2mat(outsCheck(2,6));
        RegressionOuts.AIC (facout,1) = -999;
        RegressionOuts.Constant(facout,1) = statsCheck.coeffs(1);
        RegressionOuts.B2(facout,1) = statsCheck.coeffs(length(statsCheck.coeffs));
        RegressionOuts.B3(facout,1) = -999;
        RegressionOuts.P_Constant(facout,1) = -999;
        RegressionOuts.P_B2(facout,1) = cell2mat(outsCheck(2,7));
    end
    RegressionOutputs.(genvarname([char(List.Dependent(yDepen))]))= RegressionOuts;
    clear RegressionOuts
    Values.(genvarname([char(List.Dependent(yDepen))])).marginal = yRes;
    Values.(genvarname([char(List.Dependent(yDepen))])).observed = y;
    Values.(genvarname([char(List.Dependent(yDepen))])).residuals = statsModel.resid;
end
clearvars -except List Surveys RegressionOutputs CoVar Values Statscheck

% Surveys Content:
% List.CoVar = {'SiteLocation','Age','HighestEducation','MonthlyIncome','SocialEmbedding'};
% List.Dependent = {'WatterLadderInverse','SanitationLadder','WasteRiskInverse','SocialEmbedding','MonthlyIncome','HighestEducation','RecreatioNal','RegligiousActN'};
% List.Independent = {'PopDensity','SocialEmbedding','PopGrowth','BuildingDensity','ToiletDisEuc','ToiletsRoad',...
%     'BinsEuc','BinsRoad','HealthEuc','HealthRoad','NalasEuc','NalasRoad',...
%     'SlumEuc','SlumRoad','ParksEuc','ParksRoad','PlayGroundEuc','PlayGroundRoad',...
%     'OpenSpaceEuc','OpenSpaceRoads','ForestEuc','ForestRoads','WaterEuc','WateRoads',...
%     'TownCentreEuc','TownCentreRoads'};
