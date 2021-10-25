% Define A_x, b_x
A_x = [1,0; -1,0; 0,1; 0,-1; 1,1];
b_x = [1; 1; 2; 2; 2];
% Define a YALMIP optimisation variable
E = sdpvar(2,2);

% Define objective
% Hint: Use the logdet() function
objective = -logdet(E);

% Define constraints
% Hint: The operator >= denotes matrix inequality
constraints = [E>=0];
for i=1:size(A_x,1)
    constraints = [constraints, A_x(i,:)*E*A_x(i,:)'<=b_x(i)^2];
end

% Solve the YALMIP Problem
optimize(constraints, objective)

% Display the result
disp(inv(value(E)))


%% Plot result
figure(1)
hold on
Polyhedron(A_x,b_x).plot('Color','white')
x=sdpvar(2,1);
YSet(x,[x'*inv(value(E))*x<=1]).plot('Color','b','alpha',0.5)
hold off