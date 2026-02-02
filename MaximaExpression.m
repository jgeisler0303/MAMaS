classdef MaximaExpression < handle
    % MaximaExpression
    % MATLAB class for symbolic Maxima expressions

    properties (SetAccess = private)
        identifier
    end

    properties (Access = private)
        cachedOutput
        maximaInstance  % Cached MaximaInterface instance
    end

    methods
        function obj = MaximaExpression(name, internal, maxima)
            % Constructor for new symbols
            arguments
                name
                internal (1,1) logical = false
                maxima = []
            end

            if internal
                obj.identifier = char(name);
                obj.cachedOutput = '';
                obj.maximaInstance = maxima;
                return;
            end

            if ~(ischar(name) || (isstring(name) && isscalar(name)))
                error('Symbol name must be a string.');
            end
            name = char(name);

            if ~MaximaExpression.isValidSymbolName(name)
                error('Invalid symbol name: "%s".', name);
            end

            obj.identifier = name;
            obj.cachedOutput = '';
            obj.maximaInstance = MaximaInterface.getInstance();
        end

        function s = char(obj)
            s = obj.getDisplayString();
        end

        function s = string(obj)
            s = string(obj.getDisplayString());
        end

        function disp(obj)
            disp(obj.getDisplayString());
        end

        function display(obj)
            disp(obj.getDisplayString());
        end

        % Operator Overloads
        function y = plus(a, b)
            y = MaximaExpression.binaryOp(a, b, '+');
        end

        function y = minus(a, b)
            y = MaximaExpression.binaryOp(a, b, '-');
        end

        function y = times(a, b)
            y = MaximaExpression.binaryOp(a, b, '*');
        end

        function y = mtimes(a, b)
            y = MaximaExpression.binaryOp(a, b, '*');
        end

        function y = rdivide(a, b)
            y = MaximaExpression.binaryOp(a, b, '/');
        end

        function y = mrdivide(a, b)
            y = MaximaExpression.binaryOp(a, b, '/');
        end

        function y = power(a, b)
            y = MaximaExpression.binaryOp(a, b, '^');
        end

        function y = mpower(a, b)
            y = MaximaExpression.binaryOp(a, b, '^');
        end

        function y = uminus(a)
            y = MaximaExpression.unaryOp(a, '-');
        end

        function y = uplus(a)
            y = a;
        end

        % Common math functions
        function y = sin(x)
            y = MaximaExpression.funcOp('sin', x);
        end

        function y = cos(x)
            y = MaximaExpression.funcOp('cos', x);
        end

        function y = tan(x)
            y = MaximaExpression.funcOp('tan', x);
        end

        function y = asin(x)
            y = MaximaExpression.funcOp('asin', x);
        end

        function y = acos(x)
            y = MaximaExpression.funcOp('acos', x);
        end

        function y = atan(x)
            y = MaximaExpression.funcOp('atan', x);
        end

        function y = exp(x)
            y = MaximaExpression.funcOp('exp', x);
        end

        function y = log(x)
            y = MaximaExpression.funcOp('log', x);
        end

        function y = sqrt(x)
            y = MaximaExpression.funcOp('sqrt', x);
        end

        function y = abs(x)
            y = MaximaExpression.funcOp('abs', x);
        end

        function y = ratsimp(x)
            y = MaximaExpression.funcOp('ratsimp', x);
        end

        function y = trigsimp(x)
            y = MaximaExpression.funcOp('trigsimp', x);
        end

        function y = simplify(x)
            % Calls trigsimp(ratsimp(x)) in Maxima
            exprX = MaximaExpression.toMaxima(x);
            
            % Get and validate instance from operand
            if isa(x, 'MaximaExpression')
                x.validateMaximaInstance();
                maxima = x.maximaInstance;
            else
                maxima = MaximaInterface.getInstance();
            end
            
            cmd = ['trigsimp(ratsimp(', exprX, '))']; 
            id = maxima.sendNoWait(cmd);
            y = MaximaExpression.fromId(id, maxima);
        end

        function validateMaximaInstance(obj)
            % Check if the cached Maxima instance is still valid
            if isempty(obj.maximaInstance) || ~isvalid(obj.maximaInstance)
                error('MaximaExpression: The Maxima interface instance has been deleted. Expression is invalid.');
            end
        end
    end

    methods (Static)
        function y = func(fname, varargin)
            % Apply a Maxima function to one or more arguments
            % func(fname) - function name only
            % func(fname, arg1) - single argument
            % func(fname, arg1, arg2, ...) - multiple arguments
            % func(fname, [arg1, arg2, ...]) - vector of arguments
            % func(fname, {arg1, arg2, ...}) - cell array of arguments
            
            if nargin < 1
                error('func: at least function name required');
            end
            
            maxima = [];
            if nargin == 1
                % No arguments provided
                exprArgs = '';
            else
                if nargin == 2
                    arg = varargin{1};
                    
                    if isvector(arg)
                        if isa(arg, 'MaximaExpression')
                            arg = {arg(:)}; 
                        elseif isstring(arg)
                            arg =  convertStringsToChars(arg);
                        else
                            error('Unsupported argument type for vector input.');
                        end
                    elseif ~iscell(arg)
                        arg= {arg};
                    end
                else
                    arg = varargin;
                end
                argStrs = cellfun(@MaximaExpression.toMaxima, arg, 'UniformOutput', false);
                exprArgs = strjoin(argStrs, ', ');

                for i = 1:length(arg)
                    if isa(arg{i}, 'MaximaExpression')
                        if isempty(maxima)
                            % if all instances are the same, vaidating the first is sufficient
                            arg{i}.validateMaximaInstance();
                            maxima = arg{i}.maximaInstance;
                        else
                            if maxima ~= arg{i}.maximaInstance
                                error('All MaximaExpression arguments must belong to the same MaximaInterface instance.');
                            end
                        end
                    end
                end
            end
            if isempty(maxima)
                maxima = MaximaInterface.getInstance();
            end
            
            cmd = [fname, '(', exprArgs, ')'];
            id = maxima.sendNoWait(cmd);
            y = MaximaExpression.fromId(id, maxima);
        end

        function obj = const(name)
            % Factory for common Maxima constants (with aliases)
            if ~(ischar(name) || (isstring(name) && isscalar(name)))
                error('Constant name must be a string.');
            end

            name = char(name);
            key = lower(strtrim(name));
            if startsWith(key, '%')
                key = key(2:end);
            end

            switch key
                case 'pi'
                    id = '%pi';
                case 'e'
                    id = '%e';
                case 'gamma'
                    id = '%gamma';
                case 'i'
                    id = '%i';
                case 'inf'
                    id = '%inf';
                case 'minf'
                    id = '%minf';
                case 'phi'
                    id = '%phi';
                otherwise
                    error('Unknown constant name: "%s".', name);
            end

            obj = MaximaExpression(id, true);
        end

        function obj = pi()
            obj = MaximaExpression.const('pi');
        end

        function obj = e()
            obj = MaximaExpression.const('e');
        end

        function obj = gamma()
            obj = MaximaExpression.const('gamma');
        end

        function obj = i()
            obj = MaximaExpression.const('i');
        end

        function obj = inf()
            obj = MaximaExpression.const('inf');
        end

        function obj = minf()
            obj = MaximaExpression.const('minf');
        end

        function obj = phi()
            obj = MaximaExpression.const('phi');
        end
    end

    methods (Access = private)
        function s = getDisplayString(obj)
            if ~isempty(obj.cachedOutput)
                s = obj.cachedOutput;
                return;
            end

            obj.validateMaximaInstance();
            result = obj.maximaInstance.sendAndParse(['string(', obj.identifier, ');']);
            if iscell(result)
                s = strjoin(result, newline);
            else
                s = result;
            end
            obj.cachedOutput = s;
        end
    end

    methods (Static, Access = private)
        function tf = isValidSymbolName(name)
            tf = ~isempty(regexp(name, '^[A-Za-z][A-Za-z0-9_]*$', 'once'));
        end

        function y = binaryOp(a, b, op)
            exprA = MaximaExpression.toMaxima(a);
            exprB = MaximaExpression.toMaxima(b);
            
            % Get and validate instance from first MaximaExpression operand
            maxima = MaximaExpression.getMaximaFromOperands(a, b);
            
            % no need for parentheses because operands are always symbols or maxima identifiers
            cmd = [exprA, op, exprB];
            id = maxima.sendNoWait(cmd);
            y = MaximaExpression.fromId(id, maxima);
        end

        function y = unaryOp(a, op)
            exprA = MaximaExpression.toMaxima(a);
            
            % Get and validate instance from operand
            if isa(a, 'MaximaExpression')
                a.validateMaximaInstance();
                maxima = a.maximaInstance;
            else
                maxima = MaximaInterface.getInstance();
            end
            
            % no need for parentheses because operands are always symbols or maxima identifiers
            cmd = [op, exprA];
            id = maxima.sendNoWait(cmd);
            y = MaximaExpression.fromId(id, maxima);
        end
        
        function y = funcOp(fname, arg)
            % Apply a function operation to a single argument
            % fname: function name (string)
            % arg: argument (can be MaximaExpression, number, or string)
            
            exprArg = MaximaExpression.toMaxima(arg);
            
            % Get and validate instance from operand
            if isa(arg, 'MaximaExpression')
                arg.validateMaximaInstance();
                maxima = arg.maximaInstance;
            else
                maxima = MaximaInterface.getInstance();
            end
            
            cmd = [fname, '(', exprArg, ')'];
            id = maxima.sendNoWait(cmd);
            y = MaximaExpression.fromId(id, maxima);
        end

        function maxima = getMaximaFromOperands(a, b)
            % Get MaximaInterface instance from operands and validate
            if isa(a, 'MaximaExpression') && isa(b, 'MaximaExpression')
                if a.maximaInstance ~= b.maximaInstance
                    error('Both MaximaExpression operands must belong to the same MaximaInterface instance.');
                end
                % if both instances are the same, vaidating one is sufficient
                a.validateMaximaInstance();
                maxima = a.maximaInstance;
            elseif isa(a, 'MaximaExpression')
                a.validateMaximaInstance();
                maxima = a.maximaInstance;
            elseif isa(b, 'MaximaExpression')
                b.validateMaximaInstance();
                maxima = b.maximaInstance;
            else
                maxima = MaximaInterface.getInstance();
            end
        end

        function s = toMaxima(x)
            if isa(x, 'MaximaExpression')
                s = x.identifier;
            elseif isnumeric(x) && isscalar(x)
                s = num2str(x, 16);
            elseif ischar(x) || (isstring(x) && isscalar(x))
                x = char(x);
                if ~MaximaExpression.isValidSymbolName(x)
                    error('String operands must be valid symbol names.');
                end
                s = x;
            else
                error('Unsupported operand type.');
            end
        end

        function obj = fromId(id, maxima)
            obj = MaximaExpression(id, true, maxima);
        end
    end
end
