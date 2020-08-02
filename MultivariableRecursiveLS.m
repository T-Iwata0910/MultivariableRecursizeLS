classdef MultivariableRecursiveLS < handle
    properties
        ForgettingFactor (1, 1) double
        
    end
%     properties(Access=private)
    properties
        Covariance
        Parameter
    end
    methods
        function obj = MultivariableRecursiveLS(initialParameter, initialCovariance, forgettingFactor)
            arguments
                initialParameter double
                initialCovariance double
                forgettingFactor (1, 1) double = 1
            end
            
            obj.Parameter = initialParameter;
            obj.Covariance = initialCovariance;
            obj.ForgettingFactor = forgettingFactor;
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