classdef ModelF
    properties
        PSF;
        Galaxy;
        OffsetS = [0, 0];
        OffsetG = [0, 0];
        
%         SersicIndex = 1;
%         MoffatBeta = 3;
        
        NoiseS = GaussianNoiseModel(); % Noise for star image
        NoiseG = GaussianNoiseModel(); % Noise for galaxy image
        
        ScaleS = 1e6; % Scale for star
        ScaleG = 1e7; % Scale for galaxy image
    end
    
    properties (Constant)
        Index = struct(...
            'PSFShape', 1 : 3, ...
            'PSFProfile', 4, ...
            'GalaxyShape', 5 : 7, ...
            'GalaxyProfile', [], ...
            'NoiseS', 8, ...
            'NoiseG', 9, ...
            'ScaleS', 10, ...
            'ScaleG', 11, ...
            'OffsetS', 12 : 13, ...
            'OffsetG', 14 : 15 ...
            );
        ParameterVectorLength = 15;
                
        StarRadiusMin = 0.5;
        StarRadiusMax = 10;
        StarRadius0 = 4;

        GalaxyRadiusMin = 0.5;
        GalaxyRadiusMax = 10;
        GalaxyRadius0 = 3;
        
        NoiseLevelMin = 10;
        NoiseLevelMax = 1000;
        ScaleMin = 0;
        ScaleMax = Inf;
    end
    
    properties (Dependent = true, SetAccess = public)
        ParameterVector;
    end
    
    methods
        function obj = ModelF()
            obj.PSF = EllipticObject(@MoffatF, [ModelF.StarRadius0, 0.9, 0], 3);    % MoffatBeta = 3;
            obj.Galaxy = EllipticObject(@SersicF, [ModelF.GalaxyRadius0, 0.9, 0], []);    % SersicIndex = const = 1;
            
            % Initial approximation.
            obj.NoiseS.Sigma = 200;
            obj.NoiseG.Sigma = 300;
            
            obj.ScaleS = 4000; % Scale for star
            obj.ScaleG = 15000; % Scale for galaxy image
        end
        
        function theta = get.ParameterVector(obj)
            theta = [...
                obj.PSF.ShapeParameters, ...
                obj.PSF.ProfileParameters, ...
                obj.Galaxy.ShapeParameters, ...
                obj.Galaxy.ProfileParameters, ...
                obj.NoiseS.Sigma, ...
                obj.NoiseG.Sigma, ...
                obj.ScaleS, ...
                obj.ScaleG, ...
                obj.OffsetS, ...
                obj.OffsetG
                ];
        end
        
        function this = set.NoiseS(this, value)
            if ~isa(value, 'NoiseModel')
                error('Invalid value for NoiseS');
            end
            
            this.NoiseS = value;
        end
        
        function this = set.NoiseG(this, value)
            if ~isa(value, 'NoiseModel')
                error('Invalid value for NoiseG');
            end
            
            this.NoiseG = value;
        end
        
        function this = set.ScaleS(this, value)
            if value <= 0
                error('Invalid value for ScaleS');
            end
            
            this.ScaleS = value;
        end
        
        function this = set.ScaleG(this, value)
            if value <= 0
                error('Invalid value for ScaleG');
            end
            
            this.ScaleG = value;
        end
        
        function this = set.PSF(this, value)
            if ~isa(value, 'EllipticObject')
                error('PSF must be an EllipticObject');
            end
            
            this.PSF = value;
        end
        
        function this = set.Galaxy(this, value)
            if ~isa(value, 'EllipticObject')
                error('Galaxy must be an EllipticObject');
            end
            
            this.Galaxy = value;
        end
        
        function obj = set.ParameterVector(obj, theta)
            obj.PSF.ShapeParameters = theta(ModelF.Index.PSFShape);
            obj.PSF.ProfileParameters = theta(ModelF.Index.PSFProfile);
            obj.Galaxy.ShapeParameters = theta(ModelF.Index.GalaxyShape);
            obj.Galaxy.ProfileParameters = theta(ModelF.Index.GalaxyProfile);
            obj.NoiseS.Sigma = theta(ModelF.Index.NoiseS);
            obj.NoiseG.Sigma = theta(ModelF.Index.NoiseG);
            obj.ScaleS = theta(ModelF.Index.ScaleS);
            obj.ScaleG = theta(ModelF.Index.ScaleG);
            obj.OffsetS = theta(ModelF.Index.OffsetS);
            obj.OffsetG = theta(ModelF.Index.OffsetG);
        end
        
        function img = GenerateImage(obj, noise, sz)
            %
            if nargin < 3
                sz = [48 48];
            end
            
            if nargin < 2
                noise = false;
            end

            xi = Image.SpatialFreq(sz);
            psf = obj.PSF.ComputeAt(xi);
            g = obj.Galaxy.ComputeAt(xi);
            
%             pix_s = prod(sinc(xi), 2);
%             pix_g = prod(sinc(xi), 2);
            
            S = reshape(obj.ScaleS * psf .* exp(-2 * pi * 1i * xi * obj.OffsetS.'), sz);
            G = reshape(obj.ScaleG * psf .* g .* exp(-2 * pi * 1i * xi * obj.OffsetG.'), sz);

            if noise
                S = complex(obj.NoiseS.AddNoise(real(S)), obj.NoiseS.AddNoise(imag(S)));
                G = complex(obj.NoiseG.AddNoise(real(G)), obj.NoiseG.AddNoise(imag(G)));
            end
            
            img = ImagePair(Image(S, Domain.FreqA), Image(G, Domain.FreqA));
        end
        
        function [L, grad] = LogLikelihood(this, img)
            %
%             sinc_s = prod(sinc(img.S.ValidNu), 2);
%             sinc_g = prod(sinc(img.G.ValidNu), 2);

            S = flat(img.S.ValidF);
            G = flat(img.G.ValidF);

%             pix_s = prod(sinc(img.S.ValidNu), 2);
%             pix_g = prod(sinc(img.G.ValidNu), 2);
            
            offset_s = flat(exp(-2i * pi * img.S.ValidNu * this.OffsetS.'));
            offset_g = flat(exp(-2i * pi * img.G.ValidNu * this.OffsetG.'));
            
            if nargout < 2
                psf = repmat(this.PSF.ComputeAt(img.S.ValidNu), 2, 1);
                g = repmat(this.Galaxy.ComputeAt(img.G.ValidNu), 2, 1);
                
%                 Ws_e = sinc_s.^2 .* Ws_e;
%                 Wg_e = sinc_g.^2 .* Wg_e;
                
                Ls = this.NoiseS.LogLikelihood(S, this.ScaleS * psf .* offset_s);  % star term
                Lg = this.NoiseG.LogLikelihood(G, this.ScaleG * psf .* g .* offset_g);   % galaxy term
                
            else
                [psf, ds_shape, ds_prof] = this.PSF.ComputeAt(img.S.ValidNu);
                psf = repmat(psf, 2, 1);
                ds_shape = repmat(ds_shape, 2, 1);
                ds_prof = repmat(ds_prof, 2, 1);
                
                [g, dg_shape, dg_prof] = this.Galaxy.ComputeAt(img.G.ValidNu);
                g = repmat(g, 2, 1);
                dg_shape = repmat(dg_shape, 2, 1);
                dg_prof = repmat(dg_prof, 2, 1);
                
                ds = [ds_shape, ds_prof];
                dg = [dg_shape, dg_prof];
                
                [Ls, dLs_s, dLs_noise] = this.NoiseS.LogLikelihood(S, this.ScaleS * psf .* offset_s);  % star term
                [Lg, dLg_g, dLg_noise] = this.NoiseG.LogLikelihood(G, this.ScaleG * psf .* g .* offset_g);   % galaxy term
                
                grad([ModelF.Index.PSFShape, ModelF.Index.PSFProfile]) = ...
                    dLs_s.' * (this.ScaleS * ds .* repmat(offset_s, 1, size(ds, 2))) ...
                    + dLg_g.' * (this.ScaleG * ds .* repmat(g .* offset_g, 1, size(ds, 2)));
                
                grad([ModelF.Index.GalaxyShape, ModelF.Index.GalaxyProfile]) = ...
                    + dLg_g.' * (this.ScaleG * dg .* repmat(psf .* offset_g, 1, size(dg, 2)));
                
                grad([ModelF.Index.NoiseS, ModelF.Index.NoiseG]) = [sum(dLs_noise), sum(dLg_noise)];
                
                grad(ModelF.Index.ScaleS) = dLs_s.' * (psf .* offset_s);
                grad(ModelF.Index.ScaleG) = dLg_g.' * (psf .* g .* offset_g);
                
                d_offset_s = flat(-2i * pi * img.S.ValidNu .* repmat(exp(-2i * pi * img.S.ValidNu * this.OffsetS.'), 1, 2));
                d_offset_g = flat(-2i * pi * img.G.ValidNu .* repmat(exp(-2i * pi * img.G.ValidNu * this.OffsetG.'), 1, 2));
            
                grad(ModelF.Index.OffsetS) = dLs_s.' * (this.ScaleS * repmat(psf, 1, 2) .* d_offset_s);
                grad(ModelF.Index.OffsetG) = dLg_g.' * (this.ScaleG * repmat(psf .* g, 1, 2) .* d_offset_g);
                
                if any(isinf(grad)) || any(isnan(grad))
                    error('5');
                end
            end
            
            L = sum(Ls) + sum(Lg);

            if isinf(L) || isnan(L)
                error('4');
            end
        end
        
        function e = Ellipticity(this)
            % Compute ellipticity from galaxy params.
            r = this.Galaxy.ShapeParameters(2);
            theta = this.Galaxy.ShapeParameters(3);
            
            e = (1 - r) / (1 + r) * [cos(2 * theta), -sin(2 * theta)];
        end
    end % methods
   
    methods (Static)
        function x0 = InitialValue()
            model = ModelF();
            x0 = model.ParameterVector;
        end
        
        function [lb, ub] = Bounds(finite)
            if nargin < 1
                finite = 0;
            end
            
            if finite
                angle_min = 0;
                angle_max = pi;
            else
                angle_min = -Inf;
                angle_max = Inf;
            end
            
            % Set problem lower bounds.
            lb(ModelF.Index.PSFShape) = [ModelF.StarRadiusMin, 0, angle_min];
            lb(ModelF.Index.GalaxyShape) = [ModelF.GalaxyRadiusMin, 0, angle_min];
            lb(ModelF.Index.PSFProfile) = 1.01;
            lb(ModelF.Index.GalaxyProfile) = 0.4;
            lb([ModelF.Index.NoiseS, ModelF.Index.NoiseG]) = 0.1;
            lb([ModelF.Index.ScaleS, ModelF.Index.ScaleG]) = 0.1;
            lb([ModelF.Index.OffsetS, ModelF.Index.OffsetG]) = -5;

            % Set problem upper bounds.
            ub(ModelF.Index.PSFShape) = [ModelF.StarRadiusMax, 1, angle_max];
            ub(ModelF.Index.GalaxyShape) = [ModelF.GalaxyRadiusMax, 1, angle_max];
            ub(ModelF.Index.PSFProfile) = 10.0;
            ub(ModelF.Index.GalaxyProfile) = 10;
            ub([ModelF.Index.NoiseS, ModelF.Index.NoiseG]) = Inf;
            ub([ModelF.Index.ScaleS, ModelF.Index.ScaleG]) = Inf;
            ub([ModelF.Index.OffsetS, ModelF.Index.OffsetG]) = 5;
        end
    end
end % classdef