%function dqdt = ODEs(t, y)
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


syms k12 k21 k13 k31 k02 k03 k45 k54 k46 k64 k05 k06
x = sym('x', [1 21]);
p = sym('p', [1 51]);
syms t

%kinetic parameters
global d1 d3 u1 u4  kdelay n;

% maping species concentrations to y
% x = y;

%************************************************
%   ODE Functions List
%*************************************************/

% Auxillary equations
q4F = (p(24)+p(25)*x(1)+p(26)*x(1)^2+p(27)*x(1)^3)*x(4);            %FT3p
q1F = (p(7) +p(8) *x(1)+p(9) *x(1)^2+p(10)*x(1)^3)*x(1);            %FT4p

SR3 = (p(19)*x(19))*d3; % Brain delay
SR4 = (p(1) *x(19))*d1; % Brain delay

fCIRC = 1+(p(32)/(p(31)*exp(-x(9)))-1)*(1/(1+exp(10*x(9)-55)));

SRTSH = (p(30)+p(31)*fCIRC*sin(pi/12*t-p(33)))*exp(-x(9));

fdegTSH = p(34)+p(35)/(p(36)+x(7));
fLAG = p(41)+2*x(8)^11/(p(42)^11+x(8)^11);
f4 = p(37)+5*p(37)/(1+exp(2*x(8)-7));


% Hill Functions:
H1Fast = p(13) * x(2)^n / (p(14)^n + x(2)^n);                               
H1Slow = p(15) * x(3)^n / (p(16)^n + x(3)^n);
H2Slow = p(17) * x(3)^1 / (p(18)^1 + x(3)^1);

% New Stuff for May
global D1stim D1deg
global k45Mult k54Mult k46Mult k64Mult k05Mult k06Mult;
global k12Mult k21Mult k13Mult k31Mult k03Mult;
global k52Mult k63MultD1 k63MultD2;
global D2inhibit

%% For sensitivity ODEs
k12value = p(3)*k12Mult;
k21value = p(6)*k21Mult;
k13value = p(4)*k13Mult;
k31value = p(5)*k31Mult;
k02value = p(12)*k02Mult;
k03value = k03Mult;
k45value = p(20)*k45Mult;
k54value = p(23)*k54Mult;
k46value = p(21)*k46Mult;
k64value = p(22)*k64Mult;
k05value = p(29)*k05Mult;
k06value = k06Mult;


%% ODEs
qdot12 = k12*x(2);         %p(3) = k12
qdot21 = k21*q1F;          %p(6) = k21
qdot13 = k13*x(3);         %p(4) = k13
qdot31 = k31*q1F;          %p(5) = k31
qdot02 = k02*x(2);         %p(12) = k02
qdot03 = k03*x(3);            

qdot52 = H1Fast*k52Mult;
qdot63 = H1Slow*k63MultD1 + H2Slow*k63MultD2;

qdot45 = k45*x(5);
qdot54 = k54*q4F;
qdot46 = k46*x(6);
qdot64 = k64*q4F;
qdot05 = k05*x(5);
qdot06 = k06*x(6);

qdot(1) = qdot12 + qdot13 - qdot21 - qdot31 + SR4 + p(11)*x(11) + u1;      %T4p        q1
qdot(2) = qdot21 - qdot12 - qdot02 - qdot52;                               %T4fast     q2
qdot(3) = qdot31 - qdot13 - qdot03 - qdot63;                               %T4slow     q3
qdot(4) = qdot45 + qdot46 - qdot54 - qdot64 + SR3 + p(28)*x(13) + u4;      %T3p        q4
qdot(5) = qdot54 + qdot52 - qdot45 - qdot05;                               %T3fast     q5
qdot(6) = qdot64 + qdot63 - qdot46 - qdot06;                               %T3slow     q6
qdot(7) = SRTSH-fdegTSH*x(7);                                              %TSHp       q7
qdot(8) = f4/p(38)*x(1)+p(37)/p(39)*x(4)-p(40)*x(8);                       %T3B        q8  
qdot(9) = fLAG*(x(8)-x(9));                                                %T3B LAG    q9
qdot(10) = -p(43)*x(10);                                                   %T4PILL     q10
qdot(11) =  p(43)*x(10)-(p(44)+p(11))*x(11);                               %T4GUT      q11
qdot(12) = -p(45)*x(12);                                                   %T3PILL     q12
qdot(13) =  p(45)*x(12)-(p(46)+p(28))*x(13);                               %T3GUT      q13

% Delay ODEs
qdot(14) = -kdelay*x(14) +x(7);                                             %delay1
qdot(15) = kdelay*(x(14) -x(15));                                           %delay2
qdot(16) = kdelay*(x(15) -x(16));                                           %delay3
qdot(17) = kdelay*(x(16) -x(17));                                           %delay4
qdot(18) = kdelay*(x(17) -x(18));                                           %delay5
qdot(19) = kdelay*(x(18) -x(19));                                           %delay6
H1Fast = x(2) * x(20) / p(14);
qdot(20) = D1deg * p(13) + D1stim * (x(5) - 0.0112585458130547) - D1deg * x(20); %D1 ODE 

H1Slow = x(3) * x(20) / p(16);

qdot(21) = 0; %D2 assume constant
H2Slow = ( p(17) * x(3) ) / ( p(18) + x(3) + x(3)^2 * D2inhibit );
 
% sensitivity ODEs

% x(1) = q1(T4)
% x(?) = q2(T3)       ?????????

% x(22) = dq1dk12
% x(23) = dq1dk21 
% x(24) = dq1dk13 
% x(25) = dq1dk31 
% x(26) = dq1dk02 
% x(27) = dq1dk03 
% x(28) = dq1dk45 
% x(29) = dq1dk54 
% x(30) = dq1dk46
% x(31) = dq1dk64 
% x(32) = dq1dk05
% x(33) = dq1dk06

% x(34) = dq2dk12
% x(35) = dq2dk21 
% x(36) = dq2dk13 
% x(37) = dq2dk31 
% x(38) = dq2dk02 
% x(39) = dq2dk03 
% x(40) = dq2dk45 
% x(41) = dq2dk54 
% x(42) = dq2dk46
% x(43) = dq2dk64 
% x(44) = dq2dk05
% x(45) = dq2dk06

dq1dotdq1 = diff(qdot(1), x(1));
dq1dotdq2 = diff(qdot(1), x(2));
dq1dotdq3 = diff(qdot(1), x(3));
dq1dotdq4 = diff(qdot(1), x(4));
dq1dotdq5 = diff(qdot(1), x(5));


% ODE vector
dqdt = qdot';