function dqdt = ODEsClassic(t, y)
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

dqdt = zeros(21,1); 

% maping species concentrations to y
q = y;

%************************************************
%   ODE Functions List
%*************************************************/

%% Auxillary equations
q4F = (p(24)+p(25)*q(1)+p(26)*q(1)^2+p(27)*q(1)^3)*q(4);                    %FT3p
q1F = (p(7) +p(8) *q(1)+p(9) *q(1)^2+p(10)*q(1)^3)*q(1);                    %FT4p

SR3 = (p(19)*q(19))*d3; % Brain delay
SR4 = (p(1) *q(19))*d1; % Brain delay
fCIRC = 1+(p(32)/(p(31)*exp(-q(9)))-1)*(1/(1+exp(10*q(9)-55)));
SRTSH = (p(30)+p(31)*fCIRC*sin(pi/12*t-p(33)))*exp(-q(9));
fdegTSH = p(34)+p(35)/(p(36)+q(7));
fLAG = p(41)+2*q(8)^11/(p(42)^11+q(8)^11);
f4 = p(37)+5*p(37)/(1+exp(2*q(8)-7));
%NL = p(13)/(p(14)+q(2));


%% Hill Functions:
%H1Fast = p(13) * q(2)^n / (p(14)^n + q(2)^n);                              %without D1 compartment
H1Fast = p(20) * q(2)^n / (p(14)^n + q(2)^n);                               %with D1 compartment
H1Slow = p(15) * q(3)^n / (p(16)^n + q(3)^n);
H2Slow = p(17) * p(18)^1 / (p(18)^1 + q(3)^1);
%Original Functions:
H1Fast = p(13) * q(2) / (p(14) + q(2));
H1Slow = p(15) * q(3) / (p(16) + q(3));
H2Slow = p(17) * q(3) / (p(18) + q(3));

%% Noninhibitory:
%qdot(3) = p(5)*q1F-(p(4)+p(17)/(p(18)+q(3)))*q(3)+H;                       %T4slow
%qdot(6) = p(22)*q4F+H+p(17)*q(3)/(p(18)+q(3))-p(21)*q(6);                  %T3slow

%% ODEs
qdot(1) = SR4+p(3)*q(2)+p(4)*q(3)-(p(5)+p(6))*q1F+p(11)*q(11)+u1;           %T4dot
qdot(2) = p(6)*q1F - (p(3) + p(12))*q(2) - H1Fast;                          %T4fast
qdot(3) = p(5)*q1F - p(4)*q(3) - H1Slow - H2Slow;                           %T4slow
qdot(4) = SR3+p(20)*q(5)+p(21)*q(6)-(p(22)+p(23))*q4F+p(28)*q(13)+u4;       %T3pdot
qdot(5) = p(23)*q4F + H1Fast - (p(20)+p(29))*q(5);                          %T3fast
qdot(6) = p(22)*q4F + H1Slow + H2Slow - p(21)*q(6);                         %T3slow
qdot(7) = SRTSH-fdegTSH*q(7);                                               %TSHp
qdot(8) = f4/p(38)*q(1)+p(37)/p(39)*q(4)-p(40)*q(8);                        %T3B
qdot(9) = fLAG*(q(8)-q(9));                                                 %T3B LAG
qdot(10) = -p(43)*q(10);                                                    %T4PILLdot
qdot(11) =  p(43)*q(10)-(p(44)+p(11))*q(11);                                %T4GUTdot
qdot(12) = -p(45)*q(12);                                                    %T3PILLdot
qdot(13) =  p(45)*q(12)-(p(46)+p(28))*q(13);                                %T3GUTdot

%% do we use a variable d1 or not?
qdot(20) = p(50)*(p(13) - q(20)) ; %qdot20 = K(IntoD1) - D1 * K(leavingD1). This SHOULD maintain an equilibrium.
if q(4) > y0(4)                 % if T3 levels are higher than normal
    qdot(20) = qdot(20) + q(4)*p(51); %then allow T3 Plasma to Stimulate D1 production
end
qdot(20) = 0;

%% Delay ODEs
qdot(14) = -kdelay*q(14) +q(7);                                             %delay1
qdot(15) = kdelay*(q(14) -q(15));                                           %delay2
qdot(16) = kdelay*(q(15) -q(16));                                           %delay3
qdot(17) = kdelay*(q(16) -q(17));                                           %delay4
qdot(18) = kdelay*(q(17) -q(18));                                           %delay5
qdot(19) = kdelay*(q(18) -q(19));                                           %delay6
qdot(21) = 0;

% ODE vector
dqdt = qdot';