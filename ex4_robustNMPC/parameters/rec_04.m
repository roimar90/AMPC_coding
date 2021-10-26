%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2021, ETH Zurich, {adidier, jsieber}@ethz.ch
%
% This code is only made available for students taking the advanced MPC class
% in the fall semester of 2021 (151-0371-00L) and is NOT to be distributed.
%
% Authors: Alexandre Didier, Jérôme Sieber
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function params = rec_04()
%% Control Parameters
params.ctrl.name = 'non-linear RMPC';
params.ctrl.N = 10;
params.ctrl.Q = eye(2)*100;
params.ctrl.R = 10;
params.ctrl.tightening = 'simple';

%% True System Parameters
% system paramters
params.sys.linear = false; %[bool]
params.sys.k = 4; %[1/s^2]
params.sys.k_unknown = false;
params.sys.c = 1.5; %[1/s]
params.sys.c_unknown = false;
params.sys.l = 1.3; %[m]
params.sys.g = 9.81; %[m/s^2]
params.sys.dt = 0.1; %[s]
params.sys.p = 0; % Number of estimated parameters

% state constraints
params.sys.A_x = [1,0; -1,0; 0,1; 0,-1];
params.sys.b_x = [pi/6; pi/6; pi/4; pi/4]; %[rad]
% input constraints
params.sys.A_u = [1;-1];
params.sys.b_u = [5;5]; %[1/s^2]
% noise description
params.sys.A_w = [1,0; -1,0; 0,1; 0,-1];
params.sys.b_w = [0.4*pi/180; 0.4*pi/180; 0.5*pi/180; 0.5*pi/180]; %[rad, rad, rad/s, rad/s]
% parameter description
params.sys.A_theta = [];
params.sys.b_theta = [];
params.sys.paramVariance = 0;

params.sys.noiseProcess = false;
params.sys.noiseVariance = 1;

%% Simulation Parameters
params.sim.nrSteps = 30;
params.sim.nrTraj = 50;
params.sim.x_0 = [15*pi/180;0]; %[rad]

%% Plot Parameters
params.plot.height = 350;
params.plot.width = 900;
params.plot.alpha = 0.2;
params.plot.color = [0.7216, 0.1490, 0.0039];
params.plot.lw = 1;
end