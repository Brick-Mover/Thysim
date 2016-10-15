function InitializeDataPoints()
%initializes data points as globals

global tspan time_interval
global dataPointsTime dataPointsT3 dataPointsT4 dataTimeForPlots
global dataPointsT3Red errorT3Red
global dataPointsT3Blue errorT3Blue
global dataPointsT3Pooled errorT3Pooled
global dataPointsT4Red errorT4Red
global dataPointsT4Blue errorT4Blue
global dataPointsT4Pooled errorT4Pooled
global dataPointsWeightT3 dataPointsWeightT4
global fitIndex

%%calculate tspan
time_start = 0;
time_end = 7*24; %168 hours = 7 days
time_interval = 0.01;
tspan = [time_start:time_interval:time_end];

%% Initialize Red and Blue data points
dataPointsTime=[0,1,2,4,6,24,48,96,168];
dataTimeForPlots = dataPointsTime;
dataPointsWeightT3=[1,1,1,1,1,1,1,1,1];
dataPointsWeightT4=[1,1,1,1,1,1,1,1,1];
dataPointsWeightT3Pooled=dataPointsWeightT3;
dataPointsWeightT4Pooled=dataPointsWeightT4;

% 2 is blue, 1 is red, 3 is pooled
fitWeight = 1;
dataPointsT3Red=[1.1049908, 1.0561098, 1.2216583, 1.2962135, 1.4297527, 1.3874069, 1.9845822, 1.9310332, 1.4789005]; %T3 Red
dataPointsT4Red=[83.1675, 149.50249, 253.15092, 249.00497, 259.36981, 241.54228, 193.86402, 152.40464, 103.06799]; %T4 Red
dataPointsT3Blue=[1.3403265, 1.3629332, 1.5703606, 1.7374763, 1.79526, 1.9000906, 2.3585252, 2.3050095, 1.7436115]; %T3 Blue
dataPointsT4Blue=[95.19072, 198.42455, 243.20066, 273.8806, 265.58872, 235.3234, 200.08291, 163.59867, 117.57877]; %T4 Blue
dataPointsT3Pooled=(dataPointsT3Red + dataPointsT3Blue)/2;
dataPointsT4Pooled=(dataPointsT4Red + dataPointsT4Blue)/2;

errorT4Red=[10.15884 41.69082 38.71584 45.49379 24.18903 39.60051 17.88969 15.12839 10.1579];
errorT4Blue=[23.83826 58.6422 50.83728 11.72012 17.47732 21.34463 19.42265 16.0 12.64433]; % 16.0
errorT3Red=[0.13770046 0.0627368 0.12051824 0.0612537 0.0295975 0.14185385 0.1395091 0.1242382 0.1221517];
errorT3Blue=[0.0757315 0.072874 0.1107138 0.1710828 0.1026163 0.0377608 0.3403078 0.25 0.1109819]; % 0.25

for i = 1:9
   dataPointsWeightT3Pooled(i) = sqrt((3+3-2)/(2*errorT3Red(i)^2 + 2*errorT3Blue(i)^2));
   dataPointsWeightT4Pooled(i) = sqrt((3+3-2)/(2*errorT4Red(i)^2 + 2*errorT4Blue(i)^2));
   errorT3Pooled(i) = sqrt((2*errorT3Red(i)^2 + 2*errorT3Blue(i)^2)/(3+3-2));
   errorT4Pooled(i) = sqrt((2*errorT4Red(i)^2 + 2*errorT4Blue(i)^2)/(3+3-2));
end

%% Plot Wenzel data set that Simon sent to us
global wenzelDataTime wenzelDataT3 wenzelDataT4
wenzelDataTimeT4 = [-0.0230670663333333, 1.65233373333333, 3.02302613333333, 4.0184445, 5.30914633333333, 6.46450436666667, 7.38013933333333, 8.70870766666667, 10.272803, 24.351073, 48.475298, 72.490616, 95.95985, 124.807616666667, 148.69012];
wenzelDataTimeT3 = [0.327546983333333, 1.69815353333333, 2.35567173333333, 3.6829224, 4.99693313333333, 6.32604856666667, 7.4641182, 8.79080966666667, 9.815035, 24.1086896666667, 48.3195686666667, 72.42201, 96.7957766666667, 125.261225, 148.679396666667];
wenzelDataTime = (wenzelDataTimeT4 + wenzelDataTimeT3)/2;

wenzelDataT4 = [73.5095353333333, 129.40033, 195.278063333333, 254.221873333333, 239.932876666667, 230.637306666667, 210.801353333333, 222.309916666667, 216.20337, 170.675846666667, 188.924123333333, 161.26392, 121.954496666667, 103.297638, 97.5519053333333];
wenzelDataT3 = [0.685065516666667, 0.888772786666667, 1.04152513333333, 1.06919696666667, 1.39803123333333, 1.38328583333333, 1.06524166666667, 1.10563866666667, 1.0526945, 1.34645713333333, 1.53068573333333, 1.8570044, 1.82279576666667, 1.69557886666667, 1.31559626666667];

%% Select the data points we care to fit our data to
if fitIndex == 1
    dataPointsT3 = dataPointsT3Red;
    dataPointsT4 = dataPointsT4Red;
    errorT3Norm = norm(errorT3Red);
    errorT4Norm = norm(errorT4Red);
elseif fitIndex == 2
    dataPointsT3 = dataPointsT3Blue;
    dataPointsT4 = dataPointsT4Blue;
    errorT3Norm = norm(errorT3Blue);
    errorT4Norm = norm(errorT4Blue);
elseif fitIndex == 3
    dataPointsT3 = dataPointsT3Pooled;
    dataPointsT4 = dataPointsT4Pooled;
    dataPointsWeightT3 = dataPointsWeightT3Pooled;
    dataPointsWeightT4 = dataPointsWeightT4Pooled;
elseif fitIndex == 4
    dataPointsTime = wenzelDataTime;
    dataPointsT3 = wenzelDataT3;
    dataPointsT4 = wenzelDataT4;
    dataPointsWeightT3 = ones(1,length(wenzelDataT3));
    dataPointsWeightT4 = ones(1,length(wenzelDataT4));
end

if fitWeight == 1
    if fitIndex == 1
        dataPointsWeightT3 = errorT3Norm * dataPointsWeightT3 ./ errorT3Red;
        dataPointsWeightT4 = errorT4Norm * dataPointsWeightT4 ./ errorT4Red;
    elseif fitIndex == 2
        dataPointsWeightT3 = errorT3Norm * dataPointsWeightT3 ./ errorT3Blue;
        dataPointsWeightT4 = errorT4Norm * dataPointsWeightT4 ./ errorT4Blue;
    end
end

end

