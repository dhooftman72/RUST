% needed is a list of results (datapoint,value) and a list of coordinates
% (datapoint,x,y);
clc
clear list Output
load('coordinates_Points.mat')
coordinatestmp = round(coordinates(:,(2:3)).*10000)./10000;
% load('coordinates_Roads.mat')
%  coordinatestmp = round(coordinates(:,(2:3))./2.5).*2.5;
Output=zeros(1280,3);
NrList =  coordinates(:,1);
countstore =1 ;
for count = 1:1:length(NrList)
    Output(count,1) = NrList(count);
    y = find(results(:,1)==NrList(count));
    if isempty(y) ~=1
        Output(count,2) = results(y,2);
    else
        storelistID(countstore) = NrList(count);
        storelistCountNum(countstore) = count;
        countstore = countstore+1;
    end
    clear y
end
clear count
for count = 1:1:length(storelistID)
    testing = coordinatestmp(storelistCountNum(count),1);
    if isempty(testing) ~=1
        i1 = find(coordinatestmp(:,1)==coordinatestmp(storelistCountNum(count),1));
        i2 = find(coordinatestmp(:,2)==coordinatestmp(storelistCountNum(count),2));
        equallist = intersect(i1,i2);
        resulttest =  Output(equallist,2);
        tester = find(resulttest);
        if isempty(tester)~=1
            Output(storelistCountNum(count),2) = resulttest(tester(1));
            Output(storelistCountNum(count),3) = equallist(tester(1));
        else
            testd = find(coordinates(:,4)== coordinates(storelistCountNum(count),4));
            Tester = Output(testd,2);
            Output(storelistCountNum(count),2) = nanmean(Tester(find(Tester))); %#ok<FNDSB>
            Output(storelistCountNum(count),3) = -9999;
        end
        clear i1 i2 equallist resulttest tester 
    else
        Output(storelistCountNum(count),:) = NaN;
    end
end
clearvars -except coordinates Output results
results = 1;