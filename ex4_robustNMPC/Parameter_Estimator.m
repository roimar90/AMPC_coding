classdef Parameter_Estimator
    %PARAMETER_ESTIMATOR Class for parameter estimation
    %   Class to determine parameters given inputs
    
    properties
        params % struct of parameters
    end
    
    methods
        function obj = Parameter_Estimator(sys)
            %PARAMETER_ESTIMATOR Construct an instance of this class
            %   Construct parameter estimator
            
            %%% Parse inputs %%%
            switch nargin
                case 1
                    
                otherwise 
                    error('Wrong number of inputs!')
            end
            %%%%%%%%%%%%%%%%%%%
            
            % initialize parameters
            obj.params = sys.params;
        end       
    end
end

