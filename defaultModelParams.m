%% Load Default Model Parameters

% Total Simulation Time
tSim = 10;

% Reference Input
stepTime = 0;
thetaStep = 0;

% Angle of Attack Disturbance
alphaDisturbanceValue = 0;
alphaStepTime = 0;

% Wind Disturbance
windStepTime = 0;
windDisturbanceValue =0;

% Ramp Wind Disturbance Initialization
rampDistSlope = 0;
rampDistStartTime = 0;

% Initial Conditions
theta0 = 0;
thetadot0 = 0;

% Ramp Input Initialization
rampInputStartTime = 0;
rampInputSlope = 0;

% Sensor Dynamics Included
thetaPGnum = 25;
thetaPGden = [1 25];
thetaRGnum = 21199.36;
thetaRGden = [1 116.48 21199.36];

% Set Perturbation to none
num1 = 0 ;
den1 = 1 ;
num2 = 0 ;
den2 = 1 ;

% Set Wp and WI to be static Gain of 1
Wp = tf(1,1);
WI = tf(1,1);

%% Plant Dynamics
% ----------------------
% 5 states
% ----------------------
A = [0 1 0 0 0;...
    3.9848 -0.00794 0 0 0;...
    25 0 -25 0 0;...
    0 0 0 0 1;...
    0 21199.36 0 -21199.36 -116.48];
B = [0 0 0;...
    -7.8235 -9.413 0.00245;...
    0 0 0;...
    0 0 0;...
    0 0 0];
C = [0 0 1 0 0;...
    0 0 0 1 0];
D = [0 0 0;...
    0 0 0];

% ----------------------------------
% Reduced Order Model with 2 states
% ----------------------------------
% Approximate system by first two states and first input
A2 = A(1:2,1:2);
B2 = B(1:2,1);
C2 = eye(2);
D2 = zeros(2,1);

% For simulink model
A_ol = A2;
B_ol = B2;
C_ol = C2;
