%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2021, ETH Zurich, {adidier, jsieber}@ethz.ch
% 
% This code is only made available for students taking the advanced MPC class
% in the fall semester of 2021 (151-0371-00L) and is NOT to be distributed.
%
% Authors: Alexandre Didier, Jérôme Sieber
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef System
    %SYSTEM System class implementing a segway system
    %   The System class defines the interfaces, which will be used by the
    %   controller and parameter classes implemented in the recitations.
    properties
        params %parameter struct
        X %MPT Polytope object describing the state constraints
        U %MPT Polytope object describing the input constraints
        W %MPT Polytope object describing the disturbance set
        noiseProcess %stochastic process according which disturbances are drawn
        linear %bool flag indicating if system is linear or non-linear
        A %dynamic system matrix (only linear system)
        B %input matrix (only linear system)
        f %system dynamics
        A_affine
        Omega %uncertain parameter vector
        n %state dimension
        m %input dimension
    end
    methods
        function obj = System(params)
            %SYSTEM System Class Constructor
            %   constructs system object according to the provided
            %   parameters
            
            %%% Parse inputs %%%
            switch nargin
                case 1
                    
                otherwise 
                    error('Wrong number of inputs!')
            end
            %%%%%%%%%%%%%%%%%%%
            
            obj.params = params;
            
            % define system constraints as Polytopes
            obj.X = Polyhedron(obj.params.A_x, obj.params.b_x);
            obj.U = Polyhedron(obj.params.A_u, obj.params.b_u);
            if ~isempty(obj.params.A_w) &&  ~isempty(obj.params.b_w)
                obj.W = Polyhedron(obj.params.A_w, obj.params.b_w);
            end
            if ~isempty(obj.params.A_theta) && ~isempty(obj.params.b_theta)
                obj.Omega = Polyhedron(obj.params.A_theta, obj.params.b_theta);
            end
            % get noise distribution & linearity of the system
            obj.noiseProcess = obj.params.noiseProcess;
            obj.linear = obj.params.linear;
            
            % generate update law
            c = obj.params.c;
            k = obj.params.k;
            l = obj.params.l;
            g = obj.params.g;
            dt = obj.params.dt;
            if obj.linear
                obj.A = eye(2) + dt*[0, 1; -k+g/l, -c];
                obj.B = dt*[0;1];
                obj.n = 2;
                obj.m = 1;
                
                % define A=A_affine{1}+sum_i A_affine{i}*[theta]_i
                if obj.params.k_unknown && obj.params.c_unknown
                    obj.A_affine{1} = eye(2)+dt*[0, 1; 0, 0];
                    obj.A_affine{2} = [0, 0; -1, 0];
                    obj.A_affine{3} = [0, 0; 0, -1];
                else if obj.params.k_unknown
                        obj.A_affine{1} = eye(2)+dt*[0, 1; 0, -c];
                        obj.A_affine{2} = [0, 0; -1, 0];
                    else if obj.params.c_unknown
                            obj.A_affine{1} = eye(2)+dt*[0, 1; -k+g/l, 0];
                            obj.A_affine{2} = [0, 0; 0, -1];
                        end
                    end
                end
                
                % define update law depending on disturbance type
                if obj.noiseProcess
                    obj.f = @(x,u) obj.A*x + obj.B*u + obj.noiseProcess();
                elseif ~isempty(obj.W)
                    obj.f = @(x,u) obj.A*x + obj.B*u + obj.W.randomPoint();
                else
                    obj.f = @(x,u) obj.A*x + obj.B*u;
                end
            else
                obj.A = false;
                obj.B = false;
                obj.n = 2;
                obj.m = 1;
                
                % define update law depending on disturbance type
                if obj.noiseProcess
                   obj.f = @(x,u) [x(1) + dt*x(2); x(2) + dt*(-k*x(1) - c*x(2) + sin(x(1))*g/l + u)] + obj.noiseProcess();
                elseif ~isempty(obj.W)
                   obj.f = @(x,u) [x(1) + dt*x(2); x(2) + dt*(-k*x(1) - c*x(2) + sin(x(1))*g/l + u)] + obj.W.randomPoint();
                else
                    obj.f = @(x,u) [x(1) + dt*x(2); x(2) + dt*(-k*x(1) - c*x(2) + sin(x(1))*g/l + u)];
                end
            end
        end
        function x1 = step(obj,x,u)
            %STEP Advance system from state x with input u
            
            %%% Parse inputs %%%
            switch nargin
                case 3
                    
                otherwise 
                    error('Wrong number of inputs!')
            end
            %%%%%%%%%%%%%%%%%%%
            
            x1 = obj.f(x,u);
        end
        function update_params(obj, params)
            %UPDATE_PARAMS
            
            %%% Parse inputs %%%
            switch nargin
                case 2
                    
                otherwise 
                    error('Wrong number of inputs!')
            end
            %%%%%%%%%%%%%%%%%%%
            
            error("Not implemented");
        end
    end
end