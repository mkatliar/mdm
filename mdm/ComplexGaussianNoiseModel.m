classdef ComplexGaussianNoiseModel < NoiseModel
    properties
        Sigma = 1e2;
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
        
        function img = AddNoise(this, img)
            img = img + complex(normrnd(0, this.Sigma, size(img)), normrnd(0, this.Sigma, size(img)));
        end
        
        function [L, dL_e, dL_theta] = LogLikelihood(this, s, lambda_s)
            delta = s - lambda_s;
            D = 2;
            
            L = -D * log(sqrt(2 * pi) * this.Sigma) - abs(delta).^2 ./ (2 * this.Sigma^2);
            
            if nargout > 1
                dL_e = delta ./ this.Sigma^2;
                dL_theta = -D / this.Sigma + abs(delta).^2 / this.Sigma^3;
            end
        end
        
        function [l, grad] = LogLikelihoodOptimal(~, s, h)
            %
            % Valid data range:
            valid = isfinite(s) & isfinite(h);
            s = s(valid);
            h = h(valid);
            
            % Optimal scale value:
            scale = s.' * h / sum(h.^2);      
            
            % Optimal sigma value:
            delta2 = (s - scale * h).^2;
            sigma = sqrt(mean(delta2(isfinite(delta2))));
            
            % Likelihood:
            l = - log(sqrt(2 * pi) * sigma) - delta2 ./ (2 * sigma^2);
            
            % Derivatives:
            if nargout > 1
                % В расчёте производной пренебрегаем зависимостью scale и
                % sigma от h. (!) Не совсем точно, но чем больше размер
                % векторов тем меньше вызванная этим погрешность (?)
                grad.External = scale * (s - scale * h) ./ sigma^2;
                grad.Internal = [];
            end
        end
        
        function n = InitialEstimate(this, img)
            n(this.Index.Sigma) = std(img(:));
        end
    end
end