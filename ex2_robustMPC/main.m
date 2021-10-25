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
rng(2);

%% get parameters and define system and controller
params = get_params('rec_02');
sys = System(params.sys);
ctrl = Nonlinear_MPC(sys, params.ctrl); % TODO: Implement Nonlinear_MPC in the separate file

%% simulate the closed-loop system for nrSteps time steps and nrTraj
%  trajectories starting in x_0
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
figIdx = [1,2,3,4,5,6];
plot_x('state-time', figIdx(1), x, sys.X, params.plot);
plot_x('state-state', figIdx(2), x, sys.X, params.plot);
plot_u(figIdx(3), u, sys.U, params.plot);


%% Exercise 3a
% simulate the disturbed system with multiple initial states & plot results

% add disturbance to the nonlinear system
params.sys.b_w = [0.3; 0.3; 1e-6; 1e-6];
params.plot.alpha = 0.2;
unc_sys = System(params.sys);

nrSteps = params.sim.nrSteps;
% grid state constraints
x_0 = sys.X.grid(15)';

% allocate state trajectories
unc_x = nan(nrSteps+1,size(x_0,1),size(x_0,2));
% control-loop
for i=1:size(x_0,2)
    unc_x(1,:,i) = x_0(:,i);
    for j=1:nrSteps
        try
        [sol, ~, ~] = ctrl.solve(unc_x(j,:,i)', {}, 0);
        catch
            break;
        end
        unc_x(j+1,:,i) = unc_sys.step(unc_x(j,:,i)', sol(1));
    end
end

% identify unstable trajectories
infeas_traj = zeros(size(x_0,2),1);
for i=1:size(x_0,2)
    copy = unc_x(:,:,i);
    infeas_traj(i) = size(copy(~isnan(copy)),1)/2;
end

% plot results of disturbed system
p1 = plot_x('state-time', figIdx(4), unc_x(:,:,infeas_traj == 31), sys.X, params.plot);

% highlight unstable trajectories
params.plot.color = [0.392, 0.584, 0.929];
params.plot.alpha = 0.7;

figure(figIdx(4)); hold on;
p2 = plot_x('state-time', figIdx(4), unc_x(:,:,infeas_traj < 31), sys.X, params.plot);
legend([p1(1); p2(1)], {'feasible trajectories', 'infeasible trajectories'}, 'Location', 'southeast')

% plot feasible/infeasible initial states
infeas_traj(infeas_traj < 31) = 1;
plot_rec02(figIdx(5), x_0, infeas_traj, unc_sys.X, params.plot)

%% Exercise 3b (robustness of nominal nonlinear MPC)
% Test different initial states and disturbance sizes to get a feeling for
% the robustness of nominal nonlinear MPC.

% TODO: Change the values for initial state & disturbance size and observe
% how these two parameters affect the closed-loop trajectories and the cost
% decrease. Also try to change the number of trajectories. This will
% compute multiple trajectories with different noise realizations starting
% from the same initial state.

% Note: If you run this cell multiple times, please make sure to close the
% figure (Figure 6) before running the cell again!

x_0 = [0.3;0.15]; % initial state
params.sys.b_w = [0.08; 0.08; 1e-6; 1e-6]; % disturbance size
nrTraj = 40; % number of simulation trajectories
%%%

nrSteps = params.sim.nrSteps;
unc_sys = System(params.sys);
% allocate state trajectories
unc_x = nan(nrSteps+1,size(x_0,1),nrTraj);
unc_x(1,:,:) = repmat(x_0,[1,1,nrTraj]);
opt_cost = nan(nrSteps+1,nrTraj);
% control-loop
for i=1:nrTraj
    for j=1:nrSteps
        try
        [sol, ~, ~] = ctrl.solve(unc_x(j,:,i)', {}, 0);
        catch
            break;
        end
        unc_x(j+1,:,i) = unc_sys.step(unc_x(j,:,i)', sol(1));
        % compute optimal cost
        opt_cost(j,i) = ctrl.prob.value(ctrl.objective);
    end
end

% identify unstable trajectories
infeas_traj = zeros(nrTraj,1);
for i=1:nrTraj
    copy = unc_x(:,:,i);
    infeas_traj(i) = size(copy(~isnan(copy)),1)/2;
end

% plot results
params.plot.color = [0.7216, 0.1490, 0.0039]; params.plot.alpha = 0.2;
params.plot.width = 1400;
figure(figIdx(6))
subplot(1,2,1)
p3 = plot(x_0(1)*180/pi, x_0(2)*180/pi, 'ro');
hold on;
p4 = plot_x('state-state', figIdx(6), unc_x(:,:,infeas_traj == 31), sys.X, params.plot);
hold on;
params.plot.color = [0.392, 0.584, 0.929]; params.plot.alpha = 0.7;
p5 = plot_x('state-state', figIdx(6), unc_x(:,:,infeas_traj < 31), sys.X, params.plot);
if isempty(p5)
    legend([p3(1); p4(1)], {'initial state','feasible trajectories'}, 'Location', 'northwest')
else
    legend([p3(1); p4(1); p5(1)], {'initial state', 'feasible trajectories', 'infeasible trajectories'}, 'Location', 'northwest')
end
subplot(1,2,2)
params.plot.color = [0.7216, 0.1490, 0.0039]; params.plot.alpha = 0.2;
plot(opt_cost(:,infeas_traj == 31),'color',[params.plot.color,params.plot.alpha],'linewidth',params.plot.lw)
xlabel('time'); ylabel('optimal cost'); grid();
