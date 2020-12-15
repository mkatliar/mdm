classdef ImagePair
    properties
        S;  % STAR image
        G;  % GALAXY image
    end
    
    properties (Constant)
        FileNamePattern = fullfile('..', 'mdm_images', '%s_postage', '%s', 'mdm_%s_%s_%d.png');
    end
    
    methods 
        function Show(this, show_contour)
            if nargin < 2
                show_contour = false;
            end
            
            this.S.Show(show_contour, [[2 2 1]; [2 2 3]], 'Star %s');
            this.G.Show(show_contour, [[2 2 2]; [2 2 4]], 'Galaxy %s');
        end
        
        function this = set.S(this, value)
            if ~isa(value, 'Image')
                error('Invalid argument!');
            end
            
            this.S = value;            
        end
        
        function this = set.G(this, value)
            if ~isa(value, 'Image')
                error('Invalid argument!');
            end
            
            this.G = value;            
        end
        
        function this = ImagePair(s, g)
            this.S = s;
            this.G = g;
        end
    end
    
    methods (Static)
        function imp = Load(mode, n)
            %
            ptn = strrep(ImagePair.FileNamePattern, '\', '\\');
            
            file_name = sprintf(ptn, 'star', mode, 'star', mode, n); 
            s = Image.Load(file_name);
            
            file_name = sprintf(ptn, 'galaxy', mode, 'galaxy', mode, n); 
            g = Image.Load(file_name);
            
            imp = ImagePair(s, g);
        end
    end
end