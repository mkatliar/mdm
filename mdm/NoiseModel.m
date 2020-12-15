classdef NoiseModel
    properties (Abstract, GetAccess = public, SetAccess = immutable)
        Index;
        LowerBound;
        UpperBound;
    end
    
    properties (Dependent = true, GetAccess = public, SetAccess = public)
        ParameterVector;
    end
    
    methods
        function x = get.ParameterVector(this)
            x = this.Index.Collect(this);
        end
        
        function this = set.ParameterVector(this, x)
            this = this.Index.Distribute(x, this);
        end
        
        function bg = EstimateBackground(~, img)
            bg = mean(img(:));
        end
    end
    
    methods (Abstract)
        img = AddNoise(this, img);
        % Adds noise to vector img, according to NoiseModel.
        
        [L, dL_e, dL_theta] = LogLikelihood(this, x, expected);
        % Returns LogLikelihood and it's derivative. 
        % dL_theta = dL/d ParameterVector, and
        % dL_e = dL/d\expected.
        
        n = InitialEstimate(this, img);
        % Performs rough initial estimate of noise model based on sample
        % given.
    end
end