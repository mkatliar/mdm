classdef Image
    properties (SetAccess = immutable)
        I;
        F;  % Amplitude spectrum
        ValidF;  % Part of amplitude spectrum included in analysis
        W;  % Power spectrum
        ValidW;  % Part of power spectrum included in analysis
        IV;  % Index of valid pixels in Fourier domain.
        Nu; % Spatial frequencies
        ValidNu; % Spatial frequencies having valid values of W
        Nyx;
        SourceDomain;
    end
    
    properties
        NoiseFreeI;
    end
    
    properties (Dependent = true)
        Nxy;
        Nx;
        Ny;
        Center;
    end
    
    methods 
        function c = get.Center(this)
            c = ceil((this.Nyx + 1) / 2);
        end
        
        function rv = get.Nxy(this)
            rv = this.Nyx([2 1]);
        end
        
        function rv = get.Nx(this)
            rv = this.Nyx(2);
        end
        
        function rv = get.Ny(this)
            rv = this.Nyx(1);
        end
        
        function Show(this, show_contour, spi, t)
            if nargin < 2
                show_contour = false;
            end
            
            if nargin < 3
                spi = [[2 1 1]; [2 1 2]];
            end
            
            if nargin < 4
                t = '%s';
            end
            
            n_contour = 4;
            line_spec = ':';
            
            if ~isempty(this.I)
                subplot(spi(1, 1), spi(1, 2), spi(1, 3));
                hold off;
                
                a = this.I;
                if this.SourceDomain ~= Domain.Space
                    a = real(a);
                    a(this.Center(1), this.Center(2)) = NaN;
                end
                    
                imagesc(a);
                colormap gray;
                
                axis equal ij;
                title(sprintf(t, 'I'));
                
                if show_contour && ~isempty(this.NoiseFreeI)
                    hold on;
                    contour(this.NoiseFreeI, n_contour, line_spec);
                end
            end
            
            if ~isempty(this.F)
                subplot(spi(2, 1), spi(2, 2), spi(2, 3));
                hold off;
                
                a = this.F;
                if this.SourceDomain == Domain.Space
                    a(this.Center(1), this.Center(2)) = NaN;
                end
                
                imagesc(abs(a));
                axis equal ij;
                title(sprintf(t, 'F'));
            end
        end
        
        function this = Image(X, type)
            %
            this.SourceDomain = type;
            
            % Set size            
            sz = size(X);
            
            if ~(ndims(X) == 2 && sz(1) == sz(2))
                error('Invalid argument(s)!');
            end
                        
            this.Nyx = sz;
            this.IV = Image.ValidPart(size(X));
            
            if type == Domain.Space
                this.I = X;
                this.F = fftshift(fft2(ifftshift(X)));
                this.W = this.F .* conj(this.F);
            else                    
                if type == Domain.FreqA
                    this.F = X;
                    this.W = this.F .* conj(this.F);
                    
                    a = X;
                elseif type == Domain.FreqW
                    if ~isreal(X)
                        error('W must be real!');
                    end
                    
                    this.W = X;
                    this.F = sqrt(X);
                    
                    a = this.F;
                else
                    error('Valid ''type'' valued are ''space'', ''freqA'', ''freqW''');
                end
                
                a(isnan(a)) = 0;

                if mod(sz(1), 2) == 0
                    a(:, 1) = 0;
                end

                if mod(sz(2), 2) == 0
                    a(1, :) = 0;
                end

                this.I = ifftshift(ifft2(fftshift(a)));
            end
            
            this.ValidW = this.W(this.IV);
            this.ValidF = this.F(this.IV);
            
            % Compute spatial frequencies
            this.Nu = Image.SpatialFreq(sz);
            this.ValidNu = this.Nu(this.IV, :);
        end
    end
    
    methods (Static)
        function xi = SpatialFreq(sz)
            Nxy = sz([2 1]);
            Nx = sz(2);
            Ny = sz(1);
            
            % Compute spatial frequencies
            xc = ceil((Nxy + 1) / 2);
            [xi1, xi2] = meshgrid((-xc(1) + 1 : -xc(1) + Nx) /Nx, ...
                (-xc(2) + 1 : -xc(2) + Ny) / Ny);
            
            xi = [xi1(:), xi2(:)];
        end
        
        function valid = ValidPart(sz)
            %
            % Check size
            
            if ~(numel(sz) == 2 && sz(1) == sz(2))
                error('Invalid argument(s)!');
            end
            
            % Set valid index.
            valid = logical(triu(ones(sz), 1));
                        
            ofs = mod(sz + 1, 2);
            center = ceil((sz + 1) / 2);
            
            valid(1 : center(1), 1 : center(2)) = ...
                valid(1 : center(1), 1 : center(2)) | eye(center);
            valid(1 : ofs(1), :) = 0;
            valid(:, 1 : ofs(2)) = 0;
            valid(center(1), center(2)) = 0;
        end
        
        function img = Load(file_name)
            s = double(imread(file_name));
            
            bp = s(:) == 0;
            nbp = sum(bp);
            if nbp > 0
                s(bp) = 256;
                fprintf('%d black pixels in image %s set to 256\n', nbp, file_name);
            end
            
            img = Image(s, Domain.Space);
        end
    end
end