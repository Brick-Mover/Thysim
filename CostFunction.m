function retval = CostFunction(input)
%------------------------------------------------------------
% function CostFunction
% takes:
%     V,K
%
% returns: 
%     retval
%
%------------------------------------------------------------

global options tspan time_interval;
global dataPointsTime dataPointsT3 dataPointsT4
global dataPointsWeightT3 dataPointsWeightT4
global T3conv T4conv
global p n
global y0
global iteration
global searchPoints

global newK2 D1stim D1deg newK3;
global k45Mult k54Mult k46Mult k64Mult k05Mult k06Mult;
global k12Mult k21Mult k13Mult k31Mult k02Mult k03Mult;
global k52Mult k63MultD1 k63MultD2 D2inhibit;


%% update search parameters

% STEP 5: set parameters equal to inputs
if searchPoints == 1
elseif searchPoints == 5
    k12Mult = abs(input(1));
    k21Mult = abs(input(2));
    k13Mult = abs(input(3));
    k31Mult = abs(input(4));
    k02Mult = abs(input(5));
    D2inhibit = abs(input(6));
    D1stim = abs(input(7));
    k45Mult = abs(input(8));
    k54Mult = abs(input(9));
    k46Mult = abs(input(10));
    k64Mult = abs(input(11));
    k05Mult = abs(input(12));
elseif searchPoints == 9
    k12Mult = abs(input(1));
    k21Mult = abs(input(2));
    k13Mult = abs(input(3));
    k31Mult = abs(input(4));
    k02Mult = abs(input(5));
    k03Mult = abs(input(6));
    D2inhibit = abs(input(7));
    D1stim = abs(input(8));
    k45Mult = abs(input(9));
    k54Mult = abs(input(10));
    k46Mult = abs(input(11));
    k64Mult = abs(input(12));
    k05Mult = abs(input(13));
elseif searchPoints == 10
    k12Mult = abs(input(1));
    k21Mult = abs(input(2));
    k13Mult = abs(input(3));
    k31Mult = abs(input(4));
    k02Mult = abs(input(5));
    D2inhibit = abs(input(6));
    D1stim = abs(input(7));
    k45Mult = abs(input(8));
    k54Mult = abs(input(9));
    k46Mult = abs(input(10));
    k64Mult = abs(input(11));
    k05Mult = abs(input(12));
    k06Mult = abs(input(13));
elseif searchPoints == 11
    k12Mult = abs(input(1));
    k21Mult = abs(input(2));
    k13Mult = abs(input(3));
    k31Mult = abs(input(4));
    k02Mult = abs(input(5));
    k03Mult = abs(input(6));
    D2inhibit = abs(input(7));
    D1stim = abs(input(8));
    k45Mult = abs(input(9));
    k54Mult = abs(input(10));
    k46Mult = abs(input(11));
    k64Mult = abs(input(12));
    k05Mult = abs(input(13));
    k06Mult = abs(input(14));
end

k63MultD2 = ( (6.6781e-4) * ( p(18) + 0.639 + 0.639^2 * D2inhibit ) ) / ( p(17) * 0.639 );

%% Display Iteration Counter
iteration = iteration + 1;
display(iteration)
display(input);

%% Calculate our curve
[x, y]=ode15s(@ODEs, tspan, y0, options);

%% Calculate Curve/Data Residuals
retval = [0 0];
NUMPTS = size(dataPointsTime);
t = round(dataPointsTime/time_interval+1);
for i = 1:NUMPTS(2)
   m = dataPointsWeightT3(i)*(abs(y(t(i),4)*T3conv - dataPointsT3(i)))^2;
   retval(1) = retval(1) + m^2;
   
   m = dataPointsWeightT4(i)*(abs(y(t(i),1)*T4conv - dataPointsT4(i)))^2;
   retval(2) = retval(2) + m^2;
end

% TODO: figure out how to add weights for slow compartment data
%retval(1) = retval(1) + (mean(y((end - 24/time_interval):end,6)) - );
%retval(2) = retval(2) + (mean(y((end - 24/time_interval):end,3)) - );


%divide cost by the weighted variance
T3Mean = sum(dataPointsWeightT3 .* dataPointsT3) / sum(dataPointsWeightT3);
T4Mean = sum(dataPointsWeightT4 .* dataPointsT4) / sum(dataPointsWeightT4);
SStotT3 = sum(dataPointsWeightT3 .* (dataPointsT3 - T3Mean).^2 );
SStotT4 = sum(dataPointsWeightT4 .* (dataPointsT4 - T4Mean).^2 );

%% Marquardt Levenberg can use multidimensional residuals
retval(1) = retval(1)/SStotT3;
retval(2) = retval(2)/SStotT4;


if sum(dataPointsWeightT3) == 0
    retval = retval(2);
elseif sum(dataPointsWeightT4) == 0
    retval(2) = retval(1);
    retval = retval(2);
end



%% Nelder Mead can only use single value residuals
%if sum(dataPointsWeightT3) ~= 0 && sum(dataPointsWeightT4) ~= 0
%weight = 1000;
%display([retval(2),retval(1)*weight])
%retval = 0 + retval(2) + weight * retval(1);
%end