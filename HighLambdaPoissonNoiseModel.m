classdef HighLambdaPoissonNoiseModel < NoiseModel
    properties
    end
    
    properties (GetAccess = public, SetAccess = immutable)
        Index = ParameterIndex({});
        LowerBound = [];
        UpperBound = [];
    end
     
    methods
        function img = AddNoise(~, img)
            img = normrnd(img, sqrt(img));
        end
        
        function [l, grad] = LogLikelihood(~, s, mu)
            l = - log(sqrt(2 * pi) * sqrt(mu)) - (s - mu).^2 ./ (2 * (mu));
            
            if nargout > 1
                grad.External = (s.^2 - mu .* (1 + mu)) ./ (2 * mu.^2);
                grad.Internal = [];
            end
        end
        
        function n = InitialEstimate(~, ~)
            n = [];
        end
    end
end