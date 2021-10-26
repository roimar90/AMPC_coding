%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2021, ETH Zurich, {adidier, jsieber}@ethz.ch
% 
% This code is only made available for students taking the advanced MPC class
% in the fall semester of 2021 (151-0371-00L) and is NOT to be distributed.
%
% Authors: Alexandre Didier, Jérôme Sieber
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% reset everything and configure path
setup
% fix random number generator
rng(4);

%% get parameters and define system and controller
params = get_params('rec_04');
sys = System(params.sys);
ctrl = Nonlinear_RMPC(sys, params.ctrl);

%% Exercise 1a/1b
% 1a) Implement compute_tightening in the Nonlinear_RMPC class.
% 1b) Try various values for rho and observe the produced plots. Choose a
%     value for rho making sure the resulting RPI set is contained in the
%     state constraints
rho = 0.85; %
[x_tight, u_tight, P, K, delta] = ctrl.compute_tightening(rho);
X_tight = Polyhedron(sys.X.A, sys.X.b - x_tight);
U_tight = Polyhedron(sys.U.A, sys.U.b - u_tight);

figIdx = [1,2,3,4,5,6,7,8,9];
plot_rec04(figIdx(1), params.sim.nrSteps, X_tight, U_tight, ...
           sys.X, sys.U, P, delta, params.plot)

%% Exercise 1c
% Implement the nonlinear robust MPC problem in Nonlinear_RMPC.
% Then we simulate the closed-loop system for nrSteps time steps and nrTraj
% trajectories starting in x_0.

ctrl = Nonlinear_RMPC(sys, params.ctrl, rho);
nrSteps = params.sim.nrSteps;
nrTraj = params.sim.nrTraj;
x_0 = params.sim.x_0;

% allocate state and input trajectories
x = zeros(nrSteps+1,size(x_0,1),nrTraj); 
u = zeros(nrSteps,nrTraj);
x(1,:,:) = repmat(x_0,[1,1,nrTraj]);

% control-loop
for i=1:nrTraj
    for j=1:nrSteps
        [sol, ~, errmsg] = ctrl.solve(x(j,:,i)', {}, 0);
        if errmsg
            error(errmsg);
        end
        u(j,i) = sol(1);
        x(j+1,:,i) = sys.step(x(j,:,i)', u(j,i));
    end
end

%% plot results
plot_x('state-time', figIdx(2), x, sys.X, params.plot);
plot_x('state-state', figIdx(3), x, sys.X, params.plot);
plot_u(figIdx(4), u, sys.U, params.plot);

%% Exercise 2a
% Implement compute_min_tightening in the Nonlinear_RMPC class.
% Then run this cell and observe how the RPI set and the tightenings change
% compared to the implementation of compute_tightening.
[x_tight, u_tight, P, K, delta] = ctrl.compute_min_tightening(rho);
X_tight = Polyhedron(sys.X.A, sys.X.b - x_tight);
U_tight = Polyhedron(sys.U.A, sys.U.b - u_tight);

plot_rec04(figIdx(5), params.sim.nrSteps, X_tight, U_tight, ...
           sys.X, sys.U, P, delta, params.plot)

%% Let's increase the disturbance set
params.sys.b_w = [1*pi/180; 1*pi/180; 2*pi/180; 2*pi/180]; %[rad, rad, rad/s, rad/s]
params.ctrl.tightening = 'minimize';
sys = System(params.sys);
ctrl = Nonlinear_RMPC(sys, params.ctrl);

%% Exercise 2b - Part 1
% Try again various values for rho and observe the produced plots. Choose a
% value for rho again making sure the resulting RPI set is contained in the
% state constraints
rho = 0.75;
[x_tight, u_tight, P, K, delta] = ctrl.compute_min_tightening(rho);
X_tight = Polyhedron(sys.X.A, sys.X.b - x_tight);
U_tight = Polyhedron(sys.U.A, sys.U.b - u_tight);

plot_rec04(figIdx(6), params.sim.nrSteps, X_tight, U_tight, ...
           sys.X, sys.U, P, delta, params.plot)

%% Exercise 2b - Part 2
% simulate the closed-loop system again for the nonlinear robust MPC with
% improved tightenings.
ctrl = Nonlinear_RMPC(sys, params.ctrl, rho);
nrSteps = params.sim.nrSteps;
nrTraj = params.sim.nrTraj;
x_0 = params.sim.x_0;

% allocate state and input trajectories
x = zeros(nrSteps+1,size(x_0,1),nrTraj); 
u = zeros(nrSteps,nrTraj);
x(1,:,:) = repmat(x_0,[1,1,nrTraj]);

% control-loop
for i=1:nrTraj
    for j=1:nrSteps
        [sol, ~, errmsg] = ctrl.solve(x(j,:,i)', {}, 0);
        if errmsg
            error(errmsg);
        end
        u(j,i) = sol(1);
        x(j+1,:,i) = sys.step(x(j,:,i)', u(j,i));
    end
end

%% plot results
plot_x('state-time', figIdx(7), x, sys.X, params.plot);
plot_x('state-state', figIdx(8), x, sys.X, params.plot);
plot_u(figIdx(9), u, sys.U, params.plot);
