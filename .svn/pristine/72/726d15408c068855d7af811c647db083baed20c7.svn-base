classdef PoissonNoiseModel < NoiseModel
    properties
        PhotonsPerLevel = 1;
    end
    
    properties (GetAccess = public, SetAccess = immutable)
        Index = ParameterIndex({'PhotonsPerLevel', 1});
        LowerBound = 1;
        UpperBound = 5;
    end
    
    methods
        function this = set.PhotonsPerLevel(this, value)
            if value < 1 || value > 5
                error('PhotonsPerLevel value out of range');
            end
            
            this.PhotonsPerLevel = value;
        end
        
        function img = AddNoise(this, img)
            img = poissrnd(img) / this.PhotonsPerLevel;
        end
        
        function [l, grad] = LogLikelihood(this, s, lambda_s)
            k_s = this.PhotonsPerLevel;
            l = log(k_s) + k_s * s .* log(lambda_s) - lambda_s - gammaln(k_s * s + 1);
            
            if nargout > 1
                grad.Internal = 1 / k_s + s .* (log(lambda_s) - psi(k_s * s + 1));
                grad.External = k_s * s ./ lambda_s - 1;
            end
        end
        
        function n = InitialEstimate(this, ~)
            n(this.Index.PhotonsPerLevel) = 1;
        end
        
        function bg = EstimateBackground(this, img)
            bg = mean(img(:)) * this.PhotonsPerLevel;
        end
    end
end