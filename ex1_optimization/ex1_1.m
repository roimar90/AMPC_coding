%% Quadratic Programming using YALMIP

Q = diag([0.9,0.7]);
R = 1;
A = [0.2,0.9;0.7,0];
B = [0.6;0.8];
A_x = [1,0; -1,0; 0,1; 0,-1; 1,1];
b_x = [1; 1; 2; 2; 2];
%% Main
% Define optimization variables
x = sdpvar(2,6,'full');
u = sdpvar(1,5,'full');
x_k = sdpvar(2,1,'full');
% Define objective
objective = 0;
for i=1:5
    objective = objective + x(:,i)'*Q*x(:,i) + u(:,i)'*R*u(:,i);
end
% Define constraints
constraints = [x(:,1) ==  x_k];
for i=1:5
    constraints = [constraints, x(:,i+1) == A*x(:,i) + B*u(:,i)];
    constraints = [constraints, A_x*x(:,i) <= b_x];
end
constraints = [constraints, x(:,6) == zeros(2,1)];

% Define a YALMIP Optimizer
opt_MPC = optimizer(constraints, objective, [], {x_k}, ...
    {u(:,1)});

%% Simulation 30 steps
x_k = [0; 1.8];
for k=1:30
    u_k = opt_MPC(x_k(:,k));
    x_k(:,k+1) = A*x_k(:,k) + B*u_k;
end