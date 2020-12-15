classdef PowerNoiseModel < NoiseModel
    properties
        Sigma = 10;
    end
    
    properties (GetAccess = public, SetAccess = immutable)
        Index = ParameterIndex({'Sigma', 1});
        LowerBound = 1;
        UpperBound = 50;
    end
     
    methods
        function this = set.Sigma(this, value)
            if value <= 0
                error('Sigma value out of range');
            end
            
            this.Sigma = value;
        end
        
        function y = AddNoise(this, x)
            delta = normrnd(0, this.Sigma, numel(x), 2);
            y = reshape((sqrt(x(:)) + delta(:, 1)).^2 + delta(:, 2).^2, size(x));
        end
        
        function [L, dL_e, dL_theta] = LogLikelihood(this, w, w0)
            %
            if ~all(w(:) > 0 | isnan(w(:))) || ~all(w0(:) > 0 | isnan(w0(:)))
                error('Invalid argument');
            end
            
            sigma = this.Sigma;
            ww = sqrt(w .* w0) / sigma^2;
            
            L = -(w + w0) / (2 * sigma^2) - log(2 * sigma^2) + logbesseli0(ww);
            
            if nargout > 1
                bessel_ratio = besseli(1, ww) ./ besseli(0, ww);
                dL_e = -1 / (2 * sigma^2) + sqrt(w ./ w0) .* bessel_ratio / (2 * sigma^2);
                dL_theta = (w + w0) / sigma^3 - 2 / sigma - 2 / sigma * ww .* bessel_ratio;
            end
        end
        
        function n = InitialEstimate(this, img)
            n(this.Index.Sigma) = std(img(:));
        end
    end
end