%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2021, ETH Zurich, {adidier, jsieber}@ethz.ch
% 
% This code is only made available for students taking the advanced MPC class
% in the fall semester of 2021 (151-0371-00L) and is NOT to be distributed.
%
% Authors: Alexandre Didier, Jérôme Sieber
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Nonlinear_RMPC < Controller
    %NLMPC Nonlinear MPC Controller Class
    %   Construct and solve nominal nonlinear MPC Problem
    
    properties
        % we need to define the following properties, for Ipopt to work
        % properly as a class object.
        x %state (optimization variable)
        x_0 %initial state (parameter)
        u %input (optimization variable)
        solver %NLP solver name
        
        % since this controller uses feedback, we also store K
        K % feedback gain
        
        % since we need to compute the RPI set, we additionally store the
        % system information
        sys % System object
    end
    
    methods
        function obj = Nonlinear_RMPC(sys, params, rho, solver)
            %NLMPC Construct an instance of this class
            %   Construct non-linear MPC Class and initialize solver
            
            %%% Parse inputs %%%
            switch nargin
                case 2
                    rho = 0.9;
                    solver = 'ipopt';
                    
                case 3
                    solver = 'ipopt';
                    
                case 4
                    
                otherwise 
                    error('Wrong number of inputs!')
            end
            
            %%%
            
            % call super class constructor
            obj@Controller(sys, params);
            
            % store system information
            obj.sys = sys;
            
            % initialize NLP solver name
            obj.solver = solver;
            
            %%% Initialize the non-linear MPC Controller
            % --- compute tightening ---
            switch params.tightening
                case 'simple'
                    [x_tight, u_tight, P, K, delta] = obj.compute_tightening(rho);
                    
                case 'minimize'
                    [x_tight, u_tight, P, K, delta] = obj.compute_min_tightening(rho);
                    
                otherwise
                    error("Not a valid tightening objective!");
            end
            obj.K = K;
            
            % --- initialize Opti stack ---
            obj.prob = casadi.Opti();
            
            % --- define optimization variables ---
            obj.x = obj.prob.variable(sys.n,obj.params.N+1);
            obj.x_0 = obj.prob.parameter(sys.n,1);
            obj.u = obj.prob.variable(sys.m,obj.params.N);
            
            % --- define objective ---
            objective = 0;
            for i = 1:obj.params.N
                objective = objective + obj.x(:,i)'*obj.params.Q*obj.x(:,i) + obj.u(:,i)'*obj.params.R*obj.u(:,i);
            end
            
            % --- minimize objective ---
            obj.prob.minimize(objective)
            
            % --- define constraints ---
            obj.prob.subject_to((obj.x_0 - obj.x(:,1))'*P*(obj.x_0 - obj.x(:,1)) <= delta^2);
            for i = 1:obj.params.N
                obj.prob.subject_to(obj.x(:,i+1) == sys.f(obj.x(:,i),obj.u(:,i)));
                obj.prob.subject_to(sys.X.A*obj.x(:,i)<=sys.X.b - x_tight);
                obj.prob.subject_to(sys.U.A*obj.u(:,i)<=sys.U.b - u_tight);
            end
            obj.prob.subject_to(obj.x(:,obj.params.N+1)==[0;0]);
            
            % --- setup NLP ---
            % initialize non-verbose
            obj.prob.solver(solver, struct('print_time', 0), struct('print_level', 0));          
        end
        
        function [x_tight, u_tight, P, K, delta] = compute_tightening(obj, rho)
            %COMPUTE TIGHTENING Computes an RPI set and the corresponding
            %   tightening.
            
            %%% Parse inputs %%%
            switch nargin
                case 1
                    rho = 0.9;
                    
                case 2
                    
                otherwise 
                    error('Wrong number of inputs!')
            end
            %%%
            
            % look up system parameters
            dt = obj.sys.params.dt; %sampling time
            k = obj.sys.params.k; %spring constant
            c = obj.sys.params.c; %damping constant
            l = obj.sys.params.l; %length
            g = obj.sys.params.g; %gravitational constant
            n = obj.sys.n; %state dimension
            m = obj.sys.m; %input dimension
            
            %% TODO Starts
            % define differential system matrices
            B = [0; dt];
            A_1 = [1,dt; -dt*k + dt*g/l*cos(0), 1 - dt*c];
            A_2 = [1, dt; -dt*k + dt*g/l*cos(pi/6), 1 - dt*c];
            
            %--- setup & solve offline optimization problem ---
            % define optimization variables
            E = sdpvar(n,n);
            Y = sdpvar(m,n);

            % define constraints
            constraints = [E >= eye(n)];
            constraints = [constraints, [rho^2*E, (A_1*E + B*Y)'; A_1*E + B*Y, E] >= 0];
            constraints = [constraints, [rho^2*E, (A_2*E + B*Y)'; A_2*E + B*Y, E] >= 0];
            
            % solve problem
            options = sdpsettings('solver','sedumi', 'verbose', 0);
            % only uncomment if you have Mosek installed
            %options = sdpsettings('solver','mosek', 'verbose', 0);
            optimize(constraints, trace(E), options);
            
            % recover Lyapunov function & controller
            P = inv(value(E));
            K = value(Y)*inv(value(E));
            
            %--- compute tightening ---
            % look up state, input, and disturbance sets
            X = obj.sys.X;
            U = obj.sys.U;
            W = obj.sys.W;
                        
            % compute delta
            delta = zeros(size(W.V,1),1);
            for i=1:size(W.V,1)
                delta(i) = norm(P^(1/2)*W.V(i,:)',2);
            end
            delta = max(delta)/(1 - rho);

            % compute tightening of state constraints
            x_tight = zeros(size(X.A,1),1);
            for i=1:size(X.A,1)
                x_tight(i) = norm(value(E)^(1/2)*X.A(i,:)',2);
            end
            x_tight = delta*x_tight;
            
            % compute tightening of input constraints
            u_tight = zeros(size(U.A,1),1);
            for i=1:size(U.A,1)
                u_tight(i) = norm(value(E)^(1/2)*K'*U.A(i,:)',2);
            end
            u_tight = delta*u_tight;
        end
        
        function [x_tight, u_tight, P, K, delta] = compute_min_tightening(obj, rho)
            %COMPUTE MIN TIGHTENING Computes an RPI set and the corresponding
            %   tightening, which minimizes the constraint tightening.
            
            %%% Parse inputs %%%
            switch nargin
                case 1
                    rho = 0.9;
                    
                case 2
                    
                otherwise 
                    error('Wrong number of inputs!')
            end
            %%%
            
            % look up system parameters
            dt = obj.sys.params.dt; %sampling time
            k = obj.sys.params.k; %spring constant
            c = obj.sys.params.c; %damping constant
            l = obj.sys.params.l; %length
            g = obj.sys.params.g; %gravitational constant
            n = obj.sys.n; %state dimension
            m = obj.sys.m; %input dimension
            
            % look up state, input, and disturbance sets
            X = obj.sys.X;
            U = obj.sys.U;
            W = obj.sys.W;
            
            % define differential system matrices
            B = [0; dt];
            A_1 = [1,dt; -dt*k + dt*g/l*cos(0), 1 - dt*c];
            A_2 = [1, dt; -dt*k + dt*g/l*cos(pi/6), 1 - dt*c];
            
            %--- setup & solve offline optimization problem ---
            % define optimization variables
            E = sdpvar(n,n);
            Y = sdpvar(m,n);
            gamma_x = sdpvar(size(X.A,1),1);
            gamma_u = sdpvar(size(U.A,1),1);
            gamma_w = sdpvar(1,1);
            
            % define constraints
            constraints = [E >= eye(n)];
            constraints = [constraints, [rho^2*E, (A_1*E + B*Y)'; A_1*E + B*Y, E] >= 0];
            constraints = [constraints, [rho^2*E, (A_2*E + B*Y)'; A_2*E + B*Y, E] >= 0];
            
            for i=1:size(X.A,1)
                constraints = [constraints, [gamma_x(i), X.A(i,:)*E; E'*X.A(i,:)', E] >= 0];
            end
            for i=1:size(U.A,1)
                constraints = [constraints, [gamma_u(i), U.A(i,:)*Y; Y'*U.A(i,:)', E] >= 0];
            end
            for i=1:size(W.V,1)
                constraints = [constraints, [gamma_w, W.V(i,:); W.V(i,:)', E] >= 0];
            end
            
            % solve problem
            options = sdpsettings('solver','sedumi', 'verbose', 0);
            % only uncomment if you have Mosek installed
            %options = sdpsettings('solver','mosek', 'verbose', 0);
            
            % Please note that we included here a weighting on the state
            % tightening, i.e., 50*sum(gamma_x). We did this since for this
            % specific example, the cost favours the input tightening and
            % including this weighting puts more emphasis on the state
            % tightening, therefore ensuring more balance between the two
            % terms.
            optimize(constraints, (50*sum(gamma_x) + sum(gamma_u) + (size(X.A,1) + size(U.A,1))*gamma_w)/(2*(1-rho)), options);
            
            % recover Lyapunov function & controller
            P = inv(value(E));
            K = value(Y)*inv(value(E));
            
            %--- compute tightening ---                       
            % compute delta
            delta = sqrt(value(gamma_w))/(1 - rho);

            % compute tightening of state constraints
            x_tight = delta*sqrt(value(gamma_x));
            
            % compute tightening of input constraints
            u_tight = delta*sqrt(value(gamma_u));
        end
        
        function [u, out, info] = solve(obj, x, vars, verbose)
            % call super class constructor
            [v, z, info] = solve@Controller(obj, x, vars, verbose);
            
            % apply feedback controller
            u = obj.K*(x - z(:,1)) + v(1);
            out = z;
        end
    end
end

