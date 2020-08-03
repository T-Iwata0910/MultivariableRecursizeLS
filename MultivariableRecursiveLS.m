classdef MultivariableRecursiveLS < handle
    properties
        ForgettingFactor (1, 1) double
        
    end
    properties(Access=private)
        Covariance
        Parameter
    end
    methods
        function obj = MultivariableRecursiveLS(inputSize, outputSize, options)
            arguments
                inputSize (1, 1) double {mustBePositive}
                outputSize (1, 1) double {mustBePositive}
                options.initialParameter double = double.empty
                options.initialCovariance double = double.empty
                options.ForgettingFactor (1, 1) double {mustBeNonnegative, mustBeLessThanOrEqual(options.ForgettingFactor,1)} = 1
            end
            
            if isempty(options.initialParameter)
                options.initialParameter = -0.05 + 0.1 * rand(inputSize, outputSize);
            end
            if isempty(options.initialCovariance)
                options.initialCovariance = eye(inputSize) * 1e4;
            end
            
            validateattributes(options.initialParameter, {'double'}, {'size', [inputSize, outputSize]}, 'initialParameter');
            validateattributes(options.initialCovariance, {'double'}, {'size', [inputSize, inputSize]}, 'initialCovariance');
            
            obj.Parameter = options.initialParameter;
            obj.Covariance = options.initialCovariance;
            obj.ForgettingFactor = options.ForgettingFactor;
        end
    end
    methods
        function [thetaHat, yHat] = step(obj, y, u)
            % Varidation
            % TODO: Imprement
            
            % 
            P = obj.Covariance;
            theta = obj.Parameter;
            lambda = obj.ForgettingFactor;
            
            % Calucrate estimate error
            % e[i] = y[i] - u'[i]theta[i-1]
            yHat =  theta' * u;
            err = y - yHat;
            
            % Update Parameter
            % 
            thetaHat = theta + (P * u) / (lambda + u'*P*u) * err';
            obj.Covariance = 1/lambda * (P - (P*u*u'*P)/(lambda + u'*P*u));
            
            obj.Parameter = thetaHat;
        end
    end
    
end