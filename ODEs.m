function dqdt = ODEs(t, y)
%------------------------------------------------------------
% function ODEs 
% takes:
%     t - elapsed time since the beginning of simulation
%     y - vector of current concentrations
% returns: 
%     dqdt - vector of rates of changes of concentrations
%
% ODEs can be used with matlabs odeN function: 
% [t,x] = ode23(@ODEs, tspan, y0, options);
%     tspan - vector of sampling points during time course
%     y0 - initial concentrations
%------------------------------------------------------------

% t is in hours?

%kinetic parameters
global p d1 d2 d3 d4 u1 u4  kdelay n y0;
global T3conv T4conv;

dqdt = zeros(20,1); 

% maping species concentrations to y
q = y;

%************************************************
%   ODE Functions List
%*************************************************/

%% Auxillary equations
q4F = (p(24)+p(25)*q(1)+p(26)*q(1)^2+p(27)*q(1)^3)*q(4);            %FT3p
q1F = (p(7) +p(8) *q(1)+p(9) *q(1)^2+p(10)*q(1)^3)*q(1);            %FT4p

SR3 = (p(19)*q(19))*d3; % Brain delay
SR4 = (p(1) *q(19))*d1; % Brain delay

fCIRC = 1+(p(32)/(p(31)*exp(-q(9)))-1)*(1/(1+exp(10*q(9)-55)));

SRTSH = (p(30)+p(31)*fCIRC*sin(pi/12*t-p(33)))*exp(-q(9));

fdegTSH = p(34)+p(35)/(p(36)+q(7));
fLAG = p(41)+2*q(8)^11/(p(42)^11+q(8)^11);
f4 = p(37)+5*p(37)/(1+exp(2*q(8)-7));
%NL = p(13)/(p(14)+q(2));


%% Hill Functions:
H1Fast = p(13) * q(2)^n / (p(14)^n + q(2)^n);                               
H1Slow = p(15) * q(3)^n / (p(16)^n + q(3)^n);
H2Slow = p(17) * q(3)^1 / (p(18)^1 + q(3)^1);

%% New Stuff for May
global newK2 D1stim D1deg newK3
global k45Mult k54Mult k46Mult k64Mult k05Mult k06Mult;
global k12Mult k21Mult k13Mult k31Mult k02Mult k03Mult;
global k52Mult k63MultD1 k63MultD2;

H1Fast = q(2) * q(20) / p(14);
qdot(20) = D1deg * p(13) + D1stim * (q(5) - 0.0112585458130547) - D1deg * q(20);

H1Slow = q(3) * q(20) / p(16);

global D2inhibit
qdot(21) = 0;
H2Slow = ( p(17) * q(3) ) / ( p(18) + q(3) + q(3)^2 * D2inhibit );

%% ODEs
qdot12 = p(3)*k12Mult*q(2);         %p(3) = k12
qdot21 = p(6)*k21Mult*q1F;          %p(6) = k21
qdot13 = p(4)*k13Mult*q(3);         %p(4) = k13
qdot31 = p(5)*k31Mult*q1F;          %p(5) = k31
qdot02 = p(12)*k02Mult*q(2);
qdot03 = k03Mult * q(3);

qdot52 = H1Fast*k52Mult;
qdot63 = H1Slow*k63MultD1 + H2Slow*k63MultD2;

qdot45 = p(20)*k45Mult*q(5);
qdot54 = p(23)*k54Mult*q4F;
qdot46 = p(21)*k46Mult*q(6);
qdot64 = p(22)*k64Mult*q4F;
qdot05 = p(29)*k05Mult*q(5);
qdot06 = k06Mult*q(6);

qdot(1) = qdot12 + qdot13 - qdot21 - qdot31 + SR4 + p(11)*q(11) + u1;      %T4dot
qdot(2) = qdot21 - qdot12 - qdot02 - qdot52;                               %T4fast
qdot(3) = qdot31 - qdot13 - qdot03 - qdot63;                               %T4slow
qdot(4) = qdot45 + qdot46 - qdot54 - qdot64 + SR3 + p(28)*q(13) + u4;      %T3pdot
qdot(5) = qdot54 + qdot52 - qdot45 - qdot05;                               %T3fast
qdot(6) = qdot64 + qdot63 - qdot46 - qdot06;                               %T3slow
qdot(7) = SRTSH-fdegTSH*q(7);                                              %TSHp
qdot(8) = f4/p(38)*q(1)+p(37)/p(39)*q(4)-p(40)*q(8);                       %T3B
qdot(9) = fLAG*(q(8)-q(9));                                                %T3B LAG
qdot(10) = -p(43)*q(10);                                                   %T4PILLdot
qdot(11) =  p(43)*q(10)-(p(44)+p(11))*q(11);                               %T4GUTdot
qdot(12) = -p(45)*q(12);                                                   %T3PILLdot
qdot(13) =  p(45)*q(12)-(p(46)+p(28))*q(13);                               %T3GUTdot

%% Delay ODEs
qdot(14) = -kdelay*q(14) +q(7);                                             %delay1
qdot(15) = kdelay*(q(14) -q(15));                                           %delay2
qdot(16) = kdelay*(q(15) -q(16));                                           %delay3
qdot(17) = kdelay*(q(16) -q(17));                                           %delay4
qdot(18) = kdelay*(q(17) -q(18));                                           %delay5
qdot(19) = kdelay*(q(18) -q(19));                                           %delay6

% ODE vector
dqdt = qdot';