%------------------------------------------------------------
% Thysim Model
%------------------------------------------------------------
% @author Rukan Shao
clear;
clear global;   % @Jiuru Shao important!
%kinetic parameters

%from InitializeParameters.m
global p y0
global options tspan

%counter for CostFunction.m
global iteration
iteration = 2;

%used by Plotdata
global time y yOld



%% newer model

% STEP 0: select searchMode by changing the variable below
%searchMode 0 is no search
%searchMode 1 is Marquardt Levenberg
%searchMode 2 is Nelder Mead
%searchMode 3 is simulated annealing
global searchMode searchPoints fitIndex
searchMode = 1;
searchPoints = 11;
fitIndex = 2; % 1 is red, 2 is blue, 3 is pooled, 4 is Wenzel, 5 is sum(rgw)

InitializeDataPoints()
InitializeParameters()
% set different ICs depending on whether red or blue data is being used
if fitIndex == 1
    y0(4) = 0.00543159840245776;
end
if fitIndex == 2
    y0(4) = 0.00663210649038208;
end

global newK2 D1stim D1deg newK3
global k45Mult k54Mult k46Mult k64Mult k05Mult k06Mult;
global k12Mult k21Mult k13Mult k31Mult k02Mult k03Mult;
global k52Mult k63MultD1 k63MultD2 D2inhibit;


if fitIndex == 2
    % blue data
    k12Mult = 1.86983471842468; % only remaining insensitive parameter
    k21Mult = 3.75022695555200;
    k13Mult = 0.438791628591551;
    k31Mult = 0.171970978802829;
    k02Mult = 0.428505706431140;
    k03Mult = 0;

    k52Mult = 0.01; % 0.00104672858111850;
    k63MultD1 = 6.96e-04;
	D2inhibit = 0;
    D1deg = 0.0460289892016915; 
    D1stim = 0.00638279818977864;

    k45Mult = 0.0177414391657301;
    k54Mult = 0.0911312077890824;
    k46Mult = 9.52E-02;
    k64Mult = 1.38201556543528;
    k05Mult = 0.526679372979334;
    k64Mult = 0.5;
    k06Mult = 0;
    x = 0;
    % STEP 1: Insert new startpoint "x" for the run you would like to perform here
    % If the run ends due to max iterations, you may continue the search by
    % setting "x" to the parameter set returned at the end of the run
    % k12 k21 k13 k31 k02 D2inhibit D1stim k45 k54 k46 k64 k05
    if searchPoints == 5
        % set1
        x = [1.83060699863826,3.67710388947459,0.407529831273974,0.161368050187716,0.446360455226742,0,0.000586824593491215,0.00866883371765757,0.0236071191939354,0.102128953953089,0.00114292152313195,1.14375894746932];
        % set2, 9.8914
        x = [3.32566514029077,5.78088913508297,1.24114544775457,0.542186486193182,0.595054247018666,-0.0268313774398741,-0.0409139894465678,0.105793140812043,0.0622793829523557,0.386435906286003,0.647136698406839,1.81250554039012];
        % set3, 2.0423
        x = [2.98322171664175,4.85265926079353,1.56954903990335,0.780679781658607,0.547532893361813,-0.0586452202009170,0.0387126605756299,0.0168192373233941,0.0172865152636419,0.0975831423052269,0.179439908503138,1.43161859442455];
        % set4, 1.9165
        x = [2.89251602537961,4.89443326909524,1.56477983569084,0.731685571941464,0.518990586949911,0.0669612368798001,0.0369907566852221,0.0158798656076155,0.0159638357473210,0.0941824411586663,0.185786570082403,1.45280875587848];
        % set5, 1.7947
        x = [2.88262168760867,4.95009298450830,1.44514613423136,0.660753562856723,0.517145812206416,0.134419965591523,0.0316761056293973,0.0142104389379391,0.0135439490872514,0.0864263867705712,0.204468064263803,1.33247294203970];
        % set6, 1.7280
        %x = [2.89294588667972,4.92299890685160,1.32856221130365,0.621013919790329,0.501087890752480,-0.0556878071793347,0.0326121102939624,0.0132682955712718,0.0137411827335557,0.0800411756761469,0.183308370161494,1.33053124987434];
    end  
    if searchPoints == 9
        % 0.9513
        x = [2.12979096167914,4.15424421793373,0.659080161848407,0.251172149576821,0.478226030861308,0.00251236991533257,0.0958288352517887,6.11865349282402e-09,0.00190798205113309,0.00224984775139380,-4.25634515379251e-08,0.160409684509318,0.103059025039266];
    end
    if searchPoints == 10
        % 0.9513
        x = [2.12979096167914,4.15424421793373,0.659080161848407,0.251172149576821,0.478226030861308,0.0958288352517887,6.11865349282402e-09,0.00190798205113309,0.00224984775139380,-4.25634515379251e-08,0.160409684509318,0.103059025039266,0.00827662944669750];
    end
    if searchPoints == 11
        % set1
        x = [1.73463528833347,3.50910340337092,0.352400377858623,0.161000734767735,0.257862711523098,0.000831246632411831,0,0.0114059571559135,0.00366914408580661,0.00510620616906094,0.0293254758624637,0.162807552092413,0.640816999050954,3.46446501763967e-05];
        % set2, 20.6201
        x = [1.06048788111024,2.46442691614501,0.867040549656642,0.112791520106930,0.477275742675120,0.000132678178045928,-0.000110085951236657,-0.00372586799023299,0.0343050527874047,0.0221146137602505,0.390863364971375,0.532164797986662,2.04846509517296,0.00148629172688373];
        % set3, 5.2803
        x = [1.81045792535250,3.80449726360714,0.237331561473304,0.0815795792814736,0.353295696300684,0.000693933757400121,-0.000767499895774245,0.0434373416697466,0.0197652641604250,0.0161045950738769,0.343311486267736,0.422422484122474,2.21561562369772,0.00223663355352062];
        % set4, 1.2318
        x = [2.55885015417048,4.83544356264460,0.752227340941029,0.321928706064792,0.312213845232757,-0.00788816613156630,0.00795103508410109,7.98353149045693e-10,0.00206673881315248,-0.00280912168754769,0.0276147097080319,0.178671128841340,0.106907602191113,0.0191748809020306];
        % set5, 0.9708
        x = [2.14203235161276,4.15533563243389,0.695098184465431,0.268475775667620,0.464237797300158,0.00309946129196619,-0.0721889011840092,3.39276252469079e-09,0.00193562892424693,0.00244077416688521,0.00192258492976752,0.160378200423773,0.104189485092606,0.00947837209339469];
        % set6, 0.9513
        %x = [2.12979096167914,4.15424421793373,0.659080161848407,0.251172149576821,0.478226030861308,0.00251236991533257,0.0958288352517887,6.11865349282402e-09,0.00190798205113309,0.00224984775139380,-4.25634515379251e-08,0.160409684509318,0.103059025039266,0.00827662944669750];
    end
end

if fitIndex == 1
    % red data
    k12Mult = 7.20368263299630; % only remaining insensitive parameter
    k21Mult = 14.0601031989419;
    k13Mult = 1.48634860287821;
    k31Mult = 0.512809258393033;
    k02Mult = 0.589470927298263;

    k52Mult = 0.01; % 0.00104672858111850;
    k63MultD1 = 6.96e-04;
    D2inhibit = 0;
    D1deg = 0.0460289892016915; 
    D1stim = 0.00499739790232132;

    k45Mult = 0.108811058771968;
    k54Mult = 0.587522854722568;
    k46Mult = 5.83028090728259;
    k64Mult = 0.769148876184718;
    k05Mult = 0.493175409763661;
    k64Mult = 1;
    k12Mult = 5;
    
    k03Mult = 0;
    k06Mult = 0;
    x = 0;
    % STEP 1: Insert new startpoint "x" for the run you would like to perform here
    % If the run ends due to max iterations, you may continue the search by
    % setting "x" to the parameter set returned at the end of the run
    if searchPoints == 5
        x = [20.7075326835196,36.2531624514889,2.34412423017598,1.08991310706337,0.449761697657861,0,0.0179971430012946,0.00679728390400401,0.0184509474638053,0.0745947992807474,0.168290559533564,0.782803110819559];
        % 116.1860
        x = [6.33360570995515,11.3678746008716,1.07880442332103,0.634264587985298,0.501986615785278,0.000426651484682506,0.00452504410346933,0.105907360249835,0.634354225287014,6.99216908773713,1.18546519708551,0.346708968532009];
        % 90.7107
        x = [6.55862749165728,11.1395596730341,1.49385697841929,0.955879924185005,0.465899835505641,0.000149409690044896,0.00643732999947426,0.102954997299594,0.672071063409997,7.03533553751637,0.000331458285914925,0.360470770260557];
        % just T4
        x = [7.92636363240778,15.2817505953953,1.81773225170251,0.639718035077871,0.573616408205106,0.000156311744572021,-0.00503801018883265,0.0221066015738825,0.486514559648968,1.84621295930491,0.00102966586696387,0.363169371547584];
        % 24.8750
        x = [12.6038234025158,22.5001237265460,2.10561322954398,1.00882768364537,0.503745486125755,7.79934560009886e-05,0.00816092936591971,0.00369159087494673,0.0252677790484185,0.0760118359878367,0.00279425740261994,0.442520313778177];
        % 23.8940
        x = [12.3366641299636,22.4662317260539,1.67716203140636,0.743968889597854,0.540270092865755,-0.000170606921305440,0.00686976957629031,0.00340720577472084,0.0245225493782555,0.0740692996033673,0.00768499337971234,0.397358521730076];
        % 23.7997 - anneal
        % x = [12.3148061080790,22.4932045503122,1.74576873506947,0.756054236628289,0.546758707151960,0.000877599456657782,-0.00686521272154174,0.00340720577472084,0.0248274480680594,0.0831362501021965,0.0133562302639337,0.398157060094457];
    end
    if searchPoints == 9
        % 21.9272
        x = [14.3946485449118,26.0752729489569,1.89077467366323,0.823995942694748,0.514318023025704,2.25199733414223e-06,0.0901233797939732,0.0254978661160917,0.00915806154262962,0.0165611121955703,0.0747307194003001,0.235972103101233,1.05535174487625];
    end
    if searchPoints == 10
        % 21.9272
        x = [14.3946485449118,26.0752729489569,1.89077467366323,0.823995942694748,0.514318023025704,0.0901233797939732,0.0254978661160917,0.00915806154262962,0.0165611121955703,0.0747307194003001,0.235972103101233,1.05535174487625,4.57959682445615e-06];
    end
    if searchPoints == 11
        x = [20.7075326835196,36.2531624514889,2.34412423017598,1.08991310706337,0.449761697657861,0.00222160228051881,0,0.0179971430012946,0.00679728390400401,0.0184509474638053,0.0745947992807474,0.168290559533564,0.782803110819559,3.36943139447099e-09];
        x = [12.3366641299636,22.4662317260539,1.67716203140636,0.743968889597854,0.540270092865755,0,-0.000170606921305440,0.00686976957629031,0.00340720577472084,0.0245225493782555,0.0740692996033673,0.00768499337971234,0.397358521730076,0];
        % 32.5316
        x = [6.93307099677820,13.5254150613337,1.78565669888597,0.623562961676262,0.551597335310323,0.00236928648206026,0.00394406866000398,0.00978077602282805,0.0148065790632696,0.0217934827239224,0.400421717068690,0.654208910759310,0.950700524396928,-0.00378043275212842];
        % 22.9349
        x = [10.8761568489517,20.0002933763930,1.94097913334035,0.833373796627343,0.381560425918036,-0.00423265561054392,0.0277357977274353,0.0407930877901192,0.0172031362781757,0.0181215765738485,0.0924299078459315,0.344106965141866,1.48901541003240,2.88923089691061e-06];
        % 21.9272
        x = [14.3946485449118,26.0752729489569,1.89077467366323,0.823995942694748,0.514318023025704,2.25199733414223e-06,0.0901233797939732,0.0254978661160917,0.00915806154262962,0.0165611121955703,0.0747307194003001,0.235972103101233,1.05535174487625,4.57959682445615e-06];
        % 22.078 - anneal
        % x = [14.5113773235641,26.3772167307743,1.89253382493767,0.812772128389465,0.508577292976405,-2.56600835341829e-05,1.47864232213096e-06,0.0263130121823276,0.00874567338535943,0.0156622130194128,0.0802200842007336,0.243460116215446,1.07809518621165,-0.000378970249922317];
        % 21.8338 - anneal
        % x = [14.7332677322975,26.5545886041297,1.97890571615354,0.866768217228309,0.525255668211095,3.65743118536249e-05,2.07940194054839e-06,0.0258753195806174,0.00946380362021786,0.0156880176798872,0.0707442090506123,0.249888176415776,1.05653907424810,-4.23165334881073e-07];
    end
end


% STEP 2: If you are creating a new run with a different set of parameters,
% set all parameters here
if x ~= 0
    if searchPoints == 1
    elseif searchPoints == 5 % all including 64
        k12Mult = abs(x(1));
        k21Mult = abs(x(2));
        k13Mult = abs(x(3));
        k31Mult = abs(x(4));
        k02Mult = abs(x(5));
        D2inhibit = abs(x(6));
        D1stim = abs(x(7));
        k45Mult = abs(x(8));
        k54Mult = abs(x(9));
        k46Mult = abs(x(10));
        k64Mult = abs(x(11));
        k05Mult = abs(x(12));
    elseif searchPoints == 9 % all + 03
        k12Mult = abs(x(1));
        k21Mult = abs(x(2));
        k13Mult = abs(x(3));
        k31Mult = abs(x(4));
        k02Mult = abs(x(5));
        k03Mult = abs(x(6));
        D2inhibit = abs(x(7));
        D1stim = abs(x(8));
        k45Mult = abs(x(9));
        k54Mult = abs(x(10));
        k46Mult = abs(x(11));
        k64Mult = abs(x(12));
        k05Mult = abs(x(13));
    elseif searchPoints == 10 % all + 06
        k12Mult = abs(x(1));
        k21Mult = abs(x(2));
        k13Mult = abs(x(3));
        k31Mult = abs(x(4));
        k02Mult = abs(x(5));
        D2inhibit = abs(x(6));
        D1stim = abs(x(7));
        k45Mult = abs(x(8));
        k54Mult = abs(x(9));
        k46Mult = abs(x(10));
        k64Mult = abs(x(11));
        k05Mult = abs(x(12));
        k06Mult = abs(x(13));
    elseif searchPoints == 11 % all + 03 + 06
        k12Mult = abs(x(1));
        k21Mult = abs(x(2));
        k13Mult = abs(x(3));
        k31Mult = abs(x(4));
        k02Mult = abs(x(5));
        k03Mult = abs(x(6));
        D2inhibit = abs(x(7));
        D1stim = abs(x(8));
        k45Mult = abs(x(9));
        k54Mult = abs(x(10));
        k46Mult = abs(x(11));
        k64Mult = abs(x(12));
        k05Mult = abs(x(13));
        k06Mult = abs(x(14));
    end
end

k63MultD2 = ( (6.6781e-4) * ( p(18) + 0.639 + 0.639^2 * D2inhibit ) ) / ( p(17) * 0.639 );

% STEP 3: If you are creating a new run with a new set of parameters, set
% up the startpoint array here
if searchPoints == 1
elseif searchPoints == 5
    startPoint = [k12Mult k21Mult k13Mult k31Mult k02Mult D2inhibit D1stim k45Mult k54Mult k46Mult k64Mult k05Mult];
elseif searchPoints == 9
    startPoint = [k12Mult k21Mult k13Mult k31Mult k02Mult k03Mult D2inhibit D1stim k45Mult k54Mult k46Mult k64Mult k05Mult];
elseif searchPoints == 10
    startPoint = [k12Mult k21Mult k13Mult k31Mult k02Mult D2inhibit D1stim k45Mult k54Mult k46Mult k64Mult k05Mult k06Mult];
elseif searchPoints == 11
    startPoint = [k12Mult k21Mult k13Mult k31Mult k02Mult k03Mult D2inhibit D1stim k45Mult k54Mult k46Mult k64Mult k05Mult k06Mult];
end

% run simulation using Matlab ODE Solver
if searchMode == 1
	LMoptions = optimset('Algorithm','levenberg-marquardt','TolFun',1e-6,'MaxFunEvals',30);
    % Solves nonlinear least-squares curve fitting problems
	[x,resnorm,res,eflag,output1,lambda,Jacobian] = lsqnonlin(@CostFunction,startPoint,[],[],LMoptions);
    
    Jacobian = full(Jacobian);
    varp = inv(Jacobian'*Jacobian);
    CORRMat = corrcov(varp);
end

if searchMode == 2
    fminOptions = optimset('TolFun',1e-4,'MaxIter',10000,'MaxFunEvals',30);
	[x, returnCost, exitflag, options] = fminsearch(@CostFunction,startPoint,fminOptions);
	display(options)
end

% @author Alan Chen
if searchMode == 3
    [x, returnCost] = anneal(@CostFunction, startPoint);
    % added because the final run of cost function is not the 
    % actual values that anneal decides is the best input
    CostFunction(x);
end

% y -> q
[time, y]=ode15s(@ODEs, tspan, y0, options);




%% Display...
% display('T3 muscle/serum ratio');
% display((mean(y((end - 2400):end,6)))/(mean(y((end - 2400):end,4))));
% display('T4 muscle/serum ratio');
% display((mean(y((end - 2400):end,3)))/(mean(y((end - 2400):end,1))));

PlotSetup()

plotOld = 1;
if fitIndex == 1
    PlotData('red')
    PlotExperiment('k')
end
if fitIndex == 2
    PlotData('blue')
    PlotExperiment('k')
end