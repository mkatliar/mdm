classdef Estimate
    properties (SetAccess = immutable)
        Value;
        FunctionValue;
        Gradient;
        Hessian;
        ExitFlag;
    end
    
    properties (Dependent = true)
        Sigma;
        Delta;
    end
    
    methods 
        function this = Estimate(img, options, global_opt, varargin)
            % Estimate parameters from star+galaxy image pair.

            % Initialize problem
            model = ModelF();
            
            % Set initial approximation.
            problem.x0 = model.ParameterVector.';
            problem.objective = @(x) Estimate.Objective(img, model, x);
%             problem.x0(ModelF.Index.NoiseS) = 500;
%             problem.x0(ModelF.Index.NoiseG) = 500;
            
            [lb, ub] = ModelF.Bounds();
            
            % Set problem bounds.
            problem.lb = lb.';
            problem.ub = ub.';

            % Check parameters and bounds.
            if ~all(problem.lb < problem.x0)
                error('lb < x0 condition must be satisfied!');
            end

            if ~all(problem.x0 <= problem.ub)
                error('x0 < ub condition must be satisfied!');
            end
            

%             hess_pattern = zeros(ModelF.ParameterVectorLength);
%             hess_pattern(1 : ModelF.ParameterVectorLength / 2, 1 : ModelF.ParameterVectorLength / 2) = 1;
%             hess_pattern(ModelF.ParameterVectorLength / 2 + 1 : end, ModelF.ParameterVectorLength / 2 + 1 : end) = 1;

            problem.options = optimset('Algorithm', 'interior-point', 'Display', 'notify', 'TolFun', 1e-8, ...
                'FunValCheck', 'on', 'GradObj', 'on', 'Hessian', 'off');
            
            if nargin > 1
                problem.options = optimset(problem.options, options);
            end

            problem.solver = 'fmincon';
            
            % Start optimization
            if nargin > 2
                [x, fval, exitflag, output, solutions] = run(global_opt, problem, varargin{:});
                [~, grad] = problem.objective(x);
                hessian = [];
            else
                [x, fval, exitflag, output, lambda, grad, hessian] = fmincon(problem);
            end
            
            ind = [ModelF.Index.PSFShape(3), ModelF.Index.GalaxyShape(3)];
            x(ind) = mod(x(ind), pi);   % Limit rotation angle by [0, pi)
            
            this.Value = x.';
            this.FunctionValue = fval;
            this.Gradient = grad.';
            this.Hessian = hessian;
            this.ExitFlag = exitflag;
        end
        
        function s = get.Sigma(this)
            s = inv(this.Hessian);
        end
        
        function d = get.Delta(this)
            d = sqrt(diag(this.Sigma)).';
        end
        
        function model = Model(this)
            model = ModelF();
            model.ParameterVector = this.Value;
        end
    end
    
    methods (Static, Access = private)
        function [l, grad, Hess] = Objective(img, model, x)
            %
            model.ParameterVector = x.';

            if nargout > 2
                [l, grad, Hess] = model.LogLikelihood(img);
                l = -l;
                grad = -grad.';
                Hess = -Hess;
            elseif nargout > 1
                [l, grad] = model.LogLikelihood(img);
                l = -l;
                grad = -grad.';
            else
                l = -model.LogLikelihood(img);
            end

        end
    end
end