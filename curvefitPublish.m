function [Outs,text] = curvefit(func,Inx,Iny)
Inlogx(:,1) = log10(Inx+1);
Inx(:,2) = ones(length(Iny),1).*max(double(Iny));
Inlogx(:,2) = Inx(:,2);
[text,modelFun, prior] = funcchoice(func);

MaxItR = 1;
MaxItP = 25;
PerToTake = 250;

%% Full run for Rsquare
x = Inx;
logx = Inlogx;
y = Iny;
outs.Rsquare = -999;
outs.Beta = zeros(1,9)-999;
outs = Regfunc(outs,func,x,logx,y,modelFun,prior,1);
Outs.Rsquare(1) = outs.Rsquare(1);
Outs.AIC = NaN;
if Outs.Rsquare(1) > 0
    job = createJob();
    prior = outs.Beta(1,1:3);
    if prior(3) == -999
        prior(3) = [];
    end
    createTask(job, @curvefitAIC, 1,{func,x,logx,y,modelFun,prior});
    submit(job);
    time_out_time = 30;
    waitForState(job, 'finished',time_out_time);
    name_file = 'aic.mat';
    test = exist(name_file); %#ok<*EXIST>
    if test~= 0
        load(name_file)
        Outs.AIC = Aic;
        delete(name_file)
    end
    destroy(job)
end
%% p values beta as bootstrap
clear outs
outs.Rsquare = zeros(MaxItP,1)-999;
outs.Beta = zeros(MaxItP,9)-999;
for iter = 1:MaxItP
    %display(iter)
    PerNr = (randperm(length(Inx)));
    TkPerNr = PerNr(1:PerToTake);
    x = Inx(TkPerNr,:);
    logx = Inlogx(TkPerNr,:);
    y = Iny(TkPerNr,:);
    outs = Regfunc(outs,func,x,logx,y,modelFun,prior,iter);
end
Outs.Beta = nanmean(outs.Beta,1);
clear outs
end
%% Regression functions
function outs = Regfunc(outs,func,x,logx,y,modelFun,prior,iter)
Options = statset('FunValCheck','off','Display','off','MaxIter',200,...
    'TolFun',1.0000e-4, 'TolX',1.0e-4, 'DerivStep', 6.0555e-06,'Robust',...
    'off', 'WgtFun', 'bisquare');
if func == 1 || func == 3 ||  func == 5 || func == 7
    [beta,residuals,J,~,~]  = nlinfit(x,y,modelFun,prior,'Options',Options);
else
    [beta,residuals,J,~,~]  = nlinfit(logx,y,modelFun,prior,'Options',Options);
end
% standard error & p values
ci = nlparci(beta,residuals,'Jacobian',J);
t = tinv(1-(0.05/2),length(y)-length(beta));
sebeta= (ci(:,2)-ci(:,1)) ./ (2*t);  % Standard Error
z = beta./sebeta';
pvalu = 2*(1 - normcdf(abs(z)));
%Rsquare
sstot = nansum((y-nanmean(y)).^2);
ssres = residuals.^2;
outs.Rsquare(iter,1) = max(0,(1- (nansum(ssres)/sstot))); % Rsquare
if isnan(outs.Rsquare(iter,1))==1 || isinf(outs.Rsquare(iter,1))==1
    outs.Rsquare(iter,1) = 0;
end
outs.Beta(iter,1:length(beta)) = beta';
outs.Beta(iter,4:(3+length(sebeta))) = sebeta';
outs.Beta(iter,7:(6+length(sebeta))) = pvalu;
end

function Aic = curvefitAIC(func,x,logx,y,modelFun,prior)
if func == 1 || func == 3 ||  func == 5 || func == 7
    xtmp = x;
else
    xtmp = logx;
end
ytmp = y;
test1 = find(isnan(xtmp)==1);
test2 = find(isnan(ytmp)==1);
test = unique([test1;test2]);
xtmp(test,:) = [];
ytmp(test,:) = [];

    Options = statset('FunValCheck','off','Display','off','MaxIter',100,...
    'TolFun',1.0000e-4, 'TolX',1.0e-4, 'Jacobian','off', 'DerivStep', 6.0555e-05, 'OutputFcn',' ');
    Group = grp2idx(ytmp);
    save('all')
    [~,~,stats] = nlmefit(xtmp,ytmp,Group,[],modelFun,prior,'RefineBeta0','off','Options',Options);
    Aic = stats.aic;
    save('aic.mat','Aic')
end

%% Functions to fit
function [text,modelFun, prior]  = funcchoice(func)
if func == 1
    %Linear
    prior = [0.5,0.5];
    modelFun = @(b,x) (b(2).*x(:,1)) + b(1);
    text = 'Linear';
elseif func == 2
    %Linear logaritmic
    prior = [0.5,0.5];
    modelFun = @(b,logx) (b(2).*logx(:,1)) + b(1);
    text = 'LogLinear';
elseif func == 3
    %Quadratic
    prior = [0.5,0.5,0.5];
    modelFun = @(b,x) (b(2).*x(:,1)) + ((b(3).*x(:,1).^2)) + b(1);
    text = 'Quadratic';
elseif func == 4
    %QuadraticLogaritmic
    prior = [0.5,0.5,0.5];
    modelFun = @(b,logx) (b(2).*logx(:,1)) + (b(3).*(logx(:,1).^2)) + b(1);
     text = 'LogQuadratic';
elseif func == 5
    %logistic
    prior = [0.5,0.5];
    modelFun = @(b,x) (x(:,2)./(1+exp(-b(2).*(x(:,1)-b(1)))));
     text = 'Logistic';
elseif func == 6
    % Logaritmic logistic
    prior = [0.5,0.5];
    modelFun = @(b,logx) (logx(:,2)./(1+exp(-b(2).*(logx(:,1)-b(1)))));
     text = 'LogLogistic';
elseif func == 7
     prior = [0.5,0.5];                    
     modelFun = @(b,x) b(1).*(exp((b(2).*(x(:,1)))));
    text = 'Exponential';                      
elseif func == 8
     prior = [0.5,0.5];                    
     modelFun = @(b,logx) b(1).*(exp((b(2).*(logx(:,1)))));
    text = 'LogExponential'; 
end
end
