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
ps = sym('p', [1 51]);
syms t

%kinetic parameters
global d1 d3 u1 u4 kdelay n p;

% maping species concentrations to y
% x = y;

%************************************************
%   ODE Functions List
%*************************************************/

% Auxillary equations
q4F = (ps(24)+ps(25)*x(1)+ps(26)*x(1)^2+ps(27)*x(1)^3)*x(4);            %FT3p
q1F = (ps(7) +ps(8) *x(1)+ps(9) *x(1)^2+ps(10)*x(1)^3)*x(1);            %FT4p

SR3 = (ps(19)*x(19))*d3; % Brain delay
SR4 = (ps(1) *x(19))*d1; % Brain delay

fCIRC = 1+(ps(32)/(ps(31)*exp(-x(9)))-1)*(1/(1+exp(10*x(9)-55)));

SRTSH = (ps(30)+ps(31)*fCIRC*sin(pi/12*t-ps(33)))*exp(-x(9));

fdegTSH = ps(34)+ps(35)/(ps(36)+x(7));
fLAG = ps(41)+2*x(8)^11/(ps(42)^11+x(8)^11);
f4 = ps(37)+5*ps(37)/(1+exp(2*x(8)-7));


% Hill Functions:
H1Fast = ps(13) * x(2)^n / (ps(14)^n + x(2)^n);                               
H1Slow = ps(15) * x(3)^n / (ps(16)^n + x(3)^n);
H2Slow = ps(17) * x(3)^1 / (ps(18)^1 + x(3)^1);

% New Stuff for May
global D1stim D1deg
global k45Mult k54Mult k46Mult k64Mult k05Mult k06Mult;
global k12Mult k21Mult k13Mult k31Mult k03Mult;
global k52Mult k63MultD1 k63MultD2;
global D2inhibit

%% For sensitivity ODEs
k12value = ps(3)*k12Mult;
k21value = ps(6)*k21Mult;
k13value = ps(4)*k13Mult;
k31value = ps(5)*k31Mult;
k02value = ps(12)*k02Mult;
k03value = k03Mult;
k45value = ps(20)*k45Mult;
k54value = ps(23)*k54Mult;
k46value = ps(21)*k46Mult;
k64value = ps(22)*k64Mult;
k05value = ps(29)*k05Mult;
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

qdot(1) = qdot12 + qdot13 - qdot21 - qdot31 + SR4 + ps(11)*x(11) + u1;      %T4p        q1
qdot(2) = qdot21 - qdot12 - qdot02 - qdot52;                               %T4fast     q2
qdot(3) = qdot31 - qdot13 - qdot03 - qdot63;                               %T4slow     q3
qdot(4) = qdot45 + qdot46 - qdot54 - qdot64 + SR3 + ps(28)*x(13) + u4;      %T3p        q4
qdot(5) = qdot54 + qdot52 - qdot45 - qdot05;                               %T3fast     q5
qdot(6) = qdot64 + qdot63 - qdot46 - qdot06;                               %T3slow     q6
qdot(7) = SRTSH-fdegTSH*x(7);                                              %TSHp       q7
qdot(8) = f4/ps(38)*x(1)+ps(37)/ps(39)*x(4)-ps(40)*x(8);                       %T3B        q8  
qdot(9) = fLAG*(x(8)-x(9));                                                %T3B LAG    q9
qdot(10) = -ps(43)*x(10);                                                   %T4PILL     q10
qdot(11) =  ps(43)*x(10)-(ps(44)+ps(11))*x(11);                               %T4GUT      q11
qdot(12) = -ps(45)*x(12);                                                   %T3PILL     q12
qdot(13) =  ps(45)*x(12)-(ps(46)+ps(28))*x(13);                               %T3GUT      q13

% Delay ODEs
qdot(14) = -kdelay*x(14) +x(7);                                             %delay1
qdot(15) = kdelay*(x(14) -x(15));                                           %delay2
qdot(16) = kdelay*(x(15) -x(16));                                           %delay3
qdot(17) = kdelay*(x(16) -x(17));                                           %delay4
qdot(18) = kdelay*(x(17) -x(18));                                           %delay5
qdot(19) = kdelay*(x(18) -x(19));                                           %delay6
H1Fast = x(2) * x(20) / ps(14);
qdot(20) = D1deg * ps(13) + D1stim * (x(5) - 0.0112585458130547) - D1deg * x(20); %D1 ODE 

H1Slow = x(3) * x(20) / ps(16);

qdot(21) = 0; %D2 assume constant
H2Slow = ( ps(17) * x(3) ) / ( ps(18) + x(3) + x(3)^2 * D2inhibit );
 
% sensitivity ODEs

% x(1) = q1(T4)
% x(?) = q2(T3)       ?????????

% x(22) = dq1dk12
% x(23) = dq2dk12
% x(24) = dq3dk12 
% x(25) = dq11dk12
% x(26) = dq19dk12 
% x(27) = dq1dk12 
% x(28) = dq1dk45 


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


% dq1dotdq = sym('dq1dotdq', [1 21]);
% dq2dotdq = sym('dq2dotdq', [1 21]);

for i = 1 : 21
    dq1dotdq(i) = diff(qdot(1), x(i));  % terms with index 1, 2, 3, 11, 19 are not 0
    dq2dotdq(i) = diff(qdot(2), x(i));  % terms with index 1, 2 are not zero
    dq3dotdq(i) = diff(qdot(3), x(i));  % terms with index 1, 3 are not zero
    dq11dotdq(i) = diff(qdot(11), x(i));  % terms with index 1, 3 are not zero
    dq19dotdq(i) = diff(qdot(19), x(i));  % deplay ODEs, 19 -> (18, 19), etc.
    

    dq7dotdq(i) = diff(qdot(7), x(i));  % terms with index 8, 9 are not zero
    dq9dotdq(i) = diff(qdot(9), x(i));  % terms with index 8, 9 are not zero
    dq8dotdq(i) = diff(qdot(8), x(i));  % terms with index 1, 4, 8 are not zero
    dq4dotdq(i) = diff(qdot(4), x(i));  % terms with index 1, 4, 5, 6, 13, 19 are not zero
    dq13dotdq(i) = diff(qdot(13), x(i));  % terms with index 12, 13 are not zero
    dq12dotdq(i) = diff(qdot(12), x(i));  % terms with index 12 is not zero
end

 
  display(diff(qdot(12), );
%  display(dq2dotdq);


% ODE vector
dqdt = qdot';



