function InitializeClassicParameters()
%Initializes various parameters as globals

%% used in ThysimModel.m
global p y0 kdelay 
global n
global T3conv T4conv TSHconv 

%% used in ODEs.m
global d1 d2 d3 d4 u1 u4 

d1 = 1;
d2 = 0.88;
d3 = 1;
d4 = 0.88;

u1 = 0;                     % Infusion into plasma T4
u4 = 0;                     % Infusion into plasma T3
kdelay = 5/8;               %(n-1)/k = t; n comps, t = 8hr

%% Original kinetic parameter values
p(1) = 0.00174155;        %S4
p(2) = 8;                 %tau
p(3) = 0.868;             %k12
p(4) = 0.108;             %k13
p(5) = 584;               %k31free
p(6) = 1503;              %k21free
p(7) = 0.000289;          %A
p(8) = 0.000214;          %B
p(9) = 0.000128;          %C
p(10) = -8.83*10^-6;      %D
p(11) = 0.88;             %k4absorb; originally 0.881
p(12) = 0.0189;           %k02
p(13) = 0.00998996;       %VmaxD1fast
p(14) = 2.85;             %KmD1fast
p(15) = 6.63*10^-4;       %VmaxD1slow
p(16) = 95;               %KmD1slow
p(17) = 0.00074619;       %VmaxD2slow
p(18) = 0.075;            %KmD2slow
p(19) = 3.3572*10^-4;     %S3
p(20) = 5.37;             %k45
p(21) = 0.0689;           %k46
p(22) = 127;              %k64free
p(23) = 2043;             %k54free
p(24) = 0.00395;          %a
p(25) = 0.00185;          %b
p(26) = 0.00061;          %c
p(27) = -0.000505;        %d
p(28) = 0.88;             %k3absorb; originally 0.882
p(29) = 0.207;            %k05
p(30) = 1166;             %Bzero
p(31) = 581;              %Azero
p(32) = 2.37;             %Amax
p(33) = -3.71;            %phi
p(34) = 0.53;             %kdegTSH-HYPO
p(35) = 0.037;            %VmaxTSH
p(36) = 23;               %K50TSH
p(37) = 0.118;            %k3
p(38) = 0.29;             %T4P-EU
p(39) = 0.006;            %T3P-EU
p(40) = 0.037;            %KdegT3B
p(41) = 0.0034;           %KLAG-HYPO
p(42) = 5;                %KLAG
p(43) = 1.3;              %k4dissolve
% p(44) = 0.12*d2;          %k4excrete; originally 0.119
p(44) = 0.119;
p(45) = 1.78;             %k3dissolve
% p(46) = 0.12*d4;          %k3excrete; originally 0.118
p(46) = 0.118;
p(47) = 3.2;              %Vp
p(48) = 4.3;              %VTSH

%% Conversion factors (to match thyrosim plots):
T4conv  = 777/3.2; %T4 (q1) = 777/3.2 (ug/umol)/l
T3conv  = 651/p(47); %T3 (q2) = 651/3.2 (ug/umol)/l
TSHconv = 5.6/p(48); %TSH(q7) = 5.6/4.3 (uU)/l

%% New D1 Parameters
p(49) = 1 - 0.5^(1/12); %kD1 enzyme. halflife every 12 hours
p(50) = p(49) * 4;        %K0D
p(51) = p(13) * T3conv / 4 / 1.3;   %KD4
p(52) = 1; %dummy parameter

%% Hill parameter
n = 1;

%% Initial compartment values
y0(1) = 0.322113909290887;
y0(2) = 0.201296716637979;
y0(3) = 0.638966576237812;
y0(4) = 0.00663210649038208;
y0(5) = 0.0112585458130547;
y0(6) = 0.0652960546216066;
y0(7) = 1.78829810045363;
y0(8) = 7.0572729618671;
y0(9) = 7.05714369002947;
y0(10) = 0.0965250965250965*40; %3000ug T4: IC10=3.8610?;
y0(11) = 0;
y0(12) = 0;
y0(13) = 0;
y0(14) = 3.34290132669233;
y0(15) = 3.69277705107427;
y0(16) = 3.87942611068103;
y0(17) = 3.90062380721701;
y0(18) = 3.77876195311949;
y0(19) = 3.55364904493231;
y0(20) = p(13);
y0(21) = 0;

end

