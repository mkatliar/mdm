classdef ParameterIndex < dynamicprops
    properties (SetAccess = immutable)
        VectorLength;
    end
    
    properties (SetAccess = immutable, GetAccess = private)
        Index;
    end
    
    methods
        function this = ParameterIndex(c)
            %
            this.Index = struct();
            
            i = 1;
            for k = 1 : size(c, 1)
                propname = c{k, 1};
                n = c{k, 2};
                
                P = addprop(this, propname);
                P.SetAccess = 'immutable';
                P.Dependent = true;
                
                if isa(n, 'ParameterIndex')
                    this.Index.(propname).Range = i : i + n.VectorLength - 1;
                    this.Index.(propname).Index = n;
                    i = i + n.VectorLength;
                else
                    this.Index.(propname).Range = i : i + n - 1;
                    i = i + n;
                end
                
                P.GetMethod = @(obj) obj.Index.(propname).Range;
            end
            
            this.VectorLength = i - 1;
        end
        
        function par = Distribute(this, x, par)
            %            
            if length(x) ~= this.VectorLength
                error('Invalid parameter vector length');
            end
            
            f = fields(this.Index);
            for k = 1 : length(f)
                propname = f{k};
                ind = this.Index.(propname);
                range = ind.Range;
                
                if isfield(ind, 'Index')
                    par.(propname) = ind.Index.Distribute(x(range), par.(propname));
                else
                    par.(propname) = x(range);
                end
            end
        end
        
        function x = Collect(this, par)
            %            
            x = zeros(1, this.VectorLength);            
            
            f = fields(this.Index);
            for k = 1 : length(f)
                propname = f{k};
                ind = this.Index.(propname);
                range = ind.Range;
                
                if isfield(ind, 'Index')
                    x(range) = ind.Index.Collect(par.(propname));
                else
                    x(range) = par.(propname);
                end
            end
        end
    end
end