classdef EllipticObject
    properties
        RadialProfile; % Radial profile function
        ShapeParameters = [1, 1, 0];  % Shape parameters
        ProfileParameters;
    end
    
    properties (Dependent = true)
        AMatrix;  % A = U^T U, where U = transform matrix
    end
    
    methods
        function this = EllipticObject(radial_profile, shape_parameters, profile_parameters)
            this.RadialProfile = radial_profile;
            this.ShapeParameters = shape_parameters;
            this.ProfileParameters = profile_parameters;
        end
        
        function this = set.RadialProfile(this, value)
            if ~isa(value, 'function_handle')
                error('RadialProfile must be a function handle');
            end
            
            this.RadialProfile = value;
        end
        
        function this = set.ShapeParameters(this, value)
            if ~(numel(value) == 3 && value(1) >= 0 && value(2) >= 0 && value(2) <= 1) %|| ~(value(3) >= 0 && value(3) < pi) % ~abs(value(3) <= pi/2)
                error('Invalid shape parameters value');
            end
            
            this.ShapeParameters = value;                
        end
        
        function A = get.AMatrix(this)
            a = this.ShapeParameters(1);
            r = this.ShapeParameters(2);
            theta = this.ShapeParameters(3);
            
            b = a * r;            
            A = [[(a^2 + b^2) + (a^2 - b^2) * cos(2*theta), (a^2 - b^2) * sin(2*theta)];
                [(a^2 - b^2) * sin(2*theta), (a^2 + b^2) - (a^2 - b^2) * cos(2*theta)]] / 2;
        end
        
        function r = R(this, nu)
            r = sqrt(sum((nu * this.AMatrix) .* nu, 2));
        end
        
        function [Y, dY_shape, dY_profile] = ComputeAt(this, nu)
            r = this.R(nu);
            
            if nargout > 1
                [Y, dY_r, dY_profile] = this.RadialProfile(r, this.ProfileParameters);
                
                %dA_shape = EllipticObject.ShapeParamD(dA_xi, nu, this.ShapeParameters);
                N = size(nu, 1);
                dY_A = zeros(N, 2, 2);

                for i = 1 : 2
                    for j = 1 : 2
                        dY_A(:, i, j) = dY_r .* (nu(:, i) .* nu(:, j)) ./ (2 * r);
                    end
                end

                % dAg/d(a,b,\theta)
                a = this.ShapeParameters(1);
                r = this.ShapeParameters(2);
                theta = this.ShapeParameters(3);

                dY_shape = reshape(dY_A, N, 4) * ...
                    [[a*(1+r^2-(-1+r^2)*cos(2*theta)), 2*a^2*r*sin(theta)^2, a^2*(-1+r^2)*sin(2*theta)];
                    [-a*(-1+r^2)*sin(2*theta), -a^2*r*sin(2*theta), -a^2*(-1+r^2)*cos(2*theta)];
                    [-a*(-1+r^2)*sin(2*theta), -a^2*r*sin(2*theta), -a^2*(-1+r^2)*cos(2*theta)];
                    [a*(1+r^2+(-1+r^2)*cos(2*theta)), 2*a^2*r*cos(theta)^2, -a^2*(-1+r^2)*sin(2*theta)]];
            else
                Y = this.RadialProfile(r, this.ProfileParameters);
            end
        end
    end
end