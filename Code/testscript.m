%% Spacecraft Launch Vehicle Attitude Control System Design

% Robust Multivariable Control, Spring 2018
% Authors: Jyot Buch, Alex Hayes, Sepehr Seyedi
% Group 14

%% Clear workspace
clear;clc;close all;
bdclose all;
format short g;

%% Setup Model
% Stop for updating the bode plot
updatePlotsManually = true;

% Simulink Model to be used
model = 'LaunchVehicle17b';
load_system(model);

% Default: Set the Input to Step
set_param(sprintf([model '/InputSwitch']), 'sw', '1');

% Default: Set the Wind Disturbance to Step
set_param(sprintf([model '/DistSwitch']), 'sw', '1');

%% Overall Nominal LTI Plant Linearized Dynamics

% State Space
sys = ss(A,B,C,D);
G = tf(sys);

%% Proportional + Integral Feedback Control

% Loop through various gains
kp1 = -linspace(0.1,10,25);
kp2 = -linspace(0.1,10,25);
ki1 = -linspace(0.1,10,25);
ww = logspace(-4,4,1000);
minTotalRP = inf;

% Set the weights
Wp = tf(1,1);
WI = tf(1,1);

% OutputStructure
id = 1;

tic;

for i = 1:numel(kp1)
    
    for j = 1:numel(kp2)
        
        for k = 1:numel(ki1)
            
            % Set Gains
            Kp2 = [kp1(i) kp2(j)];
            Ki2 = [ki1(k) 0];
            Kp5 = [Kp2;zeros(2)];
            Ki5 = [Ki2;zeros(2)];
            
            % fprintf('\nkp1=%f\tkp2=%f\tki1=%f',kp1(i),kp2(j),ki1(k));
            % Assign it to structure
            out(id).kp1 = kp1(i); %#ok<*SAGROW>
            out(id).kp2 = kp2(j);
            out(id).ki1 = ki1(k);
            [A_cl,B_cl,C_cl,D_cl] = linmod(model);
            eigenvalues = eig(A_cl);
            
            if real(eigenvalues)<0
                % fprintf('\tstable\n');
                %% Nominal Performance
                cltf = tf(ss(A_cl,B_cl,C_cl,D_cl));
                
                %% Single loop at a time analysis for breaking loop at input
                
                % Loop Transfer Function
                L_SLAT_U = -inv(1+cltf(1,1))*cltf(1,1);
                [Gm,Pm,Wgm,Wpm] = margin(L_SLAT_U);
                out(id).marginL_SLAT_U.Gm = Gm;
                out(id).marginL_SLAT_U.Pm = Pm;
                out(id).marginL_SLAT_U.Wgm = Wgm;
                out(id).marginL_SLAT_U.Wpm = Wpm;
                
                %% Single loop at a time analysis for breaking loop at Sensor
                
                % Loop Transfer Function
                L_SLAT_Y1 = -inv(1+cltf(2,2))*cltf(2,2);
                L_SLAT_Y2 = -inv(1+cltf(3,3))*cltf(3,3);
                [Gm,Pm,Wgm,Wpm] = margin(L_SLAT_Y1);
                out(id).marginL_SLAT_Y1.Gm = Gm;
                out(id).marginL_SLAT_Y1.Pm = Pm;
                out(id).marginL_SLAT_Y1.Wgm = Wgm;
                out(id).marginL_SLAT_Y1.Wpm = Wpm;
                
                [Gm,Pm,Wgm,Wpm] = margin(L_SLAT_Y2);
                out(id).marginL_SLAT_Y2.Gm = Gm;
                out(id).marginL_SLAT_Y2.Pm = Pm;
                out(id).marginL_SLAT_Y2.Wgm = Wgm;
                out(id).marginL_SLAT_Y2.Wpm = Wpm;
                
                %% Check Robust Performance
                % Set the bandwidth frequency as per the specification
                WbStar = 1.05;
                
                % Weighted Sensitivity Transfer Function
                As = 0.1; % Based on the low frequency asymptote
                M = 2;
                Wp = tf([1/M WbStar],[1 As*WbStar]);
                
                % Define Weighting Transfer Function WI(s) and Robust Stability
                tau = 1/12;
                r0 = 0.1;
                rInf = 2;
                WI = tf([tau r0],[tau/rInf 1]);
                
                % Get M from simulink model
                [A_cl_M,B_cl_M,C_cl_M,D_cl_M] = linmod(model);
                cltf_M = minreal(tf(ss(A_cl_M,B_cl_M,C_cl_M,D_cl_M)));
                
                % Define M matrixc
                M = [cltf_M(1,1) cltf_M(4,1); cltf_M(1,4) cltf_M(4,4)];
                
                % Consider SISO insted
                sigmaWpS = sigma(M(1,1),ww);
                sigmaWIT = sigma(M(2,2),ww);
                out(id).maxTotalRP = max(sigmaWpS + sigmaWIT);
                
                % Increement ID
                clc;
                
                if out(id).maxTotalRP < minTotalRP
                    minTotalRP = out(id).maxTotalRP;
                end
                
                if (out(id).maxTotalRP<1)
                    out(id).flag = true;
                end
                id = id + 1
                minTotalRP
            end
        end
    end
end

toc;

%% Close Simulink Model Without Saving
close_system(model,0);