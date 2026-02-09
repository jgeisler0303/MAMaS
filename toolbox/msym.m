classdef msym < handle
    % sym
    % MATLAB class for symbolic Maxima expressions

    properties (SetAccess = private)
        identifier = 'nan'                   % Maxima expression identifier (e.g., %o1, %o2, etc.), or symbol name for simple symbols or numeric expressions
        isNumber (1,1) logical = false       % True if this is a numeric value
        isConst  (1,1) logical = false       % True if this is a constant
        isSymbol (1,1) logical = false       % True if this is a symbol
        isExpression (1,1) logical = false   % True if this is a compound expression
        isMatrix (1,1) logical = false       % True if this is a Maxima matrix
        matrixRows = 0                       % Number of rows (0 if not a matrix)
        matrixCols = 0                       % Number of columns (0 if not a matrix)
    end

    properties (Access = private)
        cachedOutput = ''    % Cached string representation
        maximaInstance = []  % Cached MaximaInterface instance
    end

    methods
        function obj = msym(x, varargin)
            % Constructor for new symbols

            if nargin == 0
                error('At least one argument expected for constructing an msym.')
            end

            if isempty(x)
                % delegate to empty
                obj = msym.empty(0,1);
                return
            end

            % Parse varargin: [dim], [sym_type], [maxima]
            dim = [];
            sym_type = [];
            maxima = [];
            
            nextArgIdx = 1;
            
            % First vararg: dimension (if numeric)
            if nextArgIdx <= length(varargin) && isnumeric(varargin{nextArgIdx})
                dim = varargin{nextArgIdx};
                nextArgIdx = nextArgIdx + 1;
            end
            
            % Next vararg: sym_type (if char/string/cell)
            if nextArgIdx <= length(varargin)
                arg = varargin{nextArgIdx};
                if ischar(arg) || isstring(arg) || (iscell(arg) && all(cellfun(@(x) ischar(x) || isstring(x), arg)))
                    sym_type = string(arg);
                    nextArgIdx = nextArgIdx + 1;
                end
            end
            
            % Final vararg: maxima instance
            if nextArgIdx <= length(varargin)
                maxima = varargin{nextArgIdx};
                if ~isa(maxima, 'MaximaInterface')
                    error('The maxima argument must be a MaximaInterface instance.');
                end
                nextArgIdx = nextArgIdx + 1;
            end
            
            % Check for extra arguments
            if nextArgIdx <= length(varargin)
                error('Too many input arguments.');
            end
            
            % Set default maxima if not provided
            if isempty(maxima)
                maxima = MaximaInterface.getInstance();
            end
            
            % Process dimension
            if ~isempty(dim)
                if all(dim == 0) % Special case: zero dimension means scalar
                    dim = [];
                else
                    if ~isvector(dim) || any(dim < 1) || any(mod(dim, 1) ~= 0)
                        error('Dimension must be a positive integer or vector of positive integers.');
                    end
                    if isscalar(dim)
                        dim = [dim, dim];
                    elseif numel(dim) > 2
                        error('Dimension vector must have at most 2 elements.');
                    end
                end
            end

            if ~isempty(sym_type) && (ismember("expr", sym_type) || ismember("const", sym_type))
                obj.identifier = char(x);
                obj.cachedOutput = '';
                obj.maximaInstance = maxima;
                obj.isMatrix = prod(dim) > 1;
                if ~isempty(dim)
                    obj.matrixRows = dim(1);
                    obj.matrixCols = dim(2);
                end
                % For internal expressions, assume they are expressions (result of operations)
                if ismember("const", sym_type)
                    obj.isConst = true;
                else
                    obj.isExpression = true;
                end
                return;
            end

            % Check if input is matrix data (numeric array, string array, or cell array)
            if ~ischar(x) && ~isscalar_matlab(x)
                % TODO: error if dim is also provided?
                % TODO: validate sym_type?

                % Delegate to matrix creation
                obj = msym.matrix(x, maxima);
                return;
            end

            % Check if the string represents a number
            numValue = str2double(x);
            if ~isnan(numValue)
                obj.identifier = x;
                obj.isNumber = true;
                return
            end

            if isnumeric(x)
                % For numeric input, dim is not allowed
                if ~all(dim == 1)
                    error('Dimension argument cannot be used when name is numeric.');
                end
                
                % Validate and process sym_type for numeric input (limited options)
                if isempty(sym_type)
                    sym_type = 'r'; % Default to rational
                end
                if numel(sym_type) > 1
                    error('For numeric input, sym_type must be a single character: "r", "d", "e", or "f".');
                end
                % Handle different numeric formats
                switch char(sym_type)
                    case 'r'
                        % Rational representation
                        obj.identifier = ['rationalize(', num2str(x, 16), ')'];
                    case 'd'
                        % Decimal/float representation
                        obj.identifier = num2str(x, 16);
                    case 'e'
                        % Exact representation (keep as rational if possible)
                        obj.identifier = num2str(x, 16);
                    case 'f'
                        % Force floating point
                        obj.identifier = ['float(', num2str(x, 16), ')'];
                    otherwise
                        error('For numeric input, sym_type must be "r" (rational), "d" (decimal), "e" (exact), or "f" (float).');
                end
                
                obj.maximaInstance = maxima;
                obj.isNumber = true;
                return;
            end

            if ~ischar(x) && ~isstring(x)
                error('Input must be a string or numeric.');
            end

            x = char(x);
            if ~msym.isValidSymbolName(x)
                error('Invalid symbol name: "%s".', x);
            end

            % Handle dimension: create matrix or scalar symbol
            if ~isempty(dim)
                % Create symbolic matrix
                matObj = msym.matrixSymbolic(x, dim(1), dim(2), maxima);
                obj.identifier = matObj.identifier;
                obj.cachedOutput = matObj.cachedOutput;
                obj.maximaInstance = matObj.maximaInstance;
                obj.isMatrix = matObj.isMatrix;
                obj.matrixRows = matObj.matrixRows;
                obj.matrixCols = matObj.matrixCols;
                obj.isExpression = matObj.isExpression;
                obj.isSymbol = false; % Matrix is not a simple symbol
            else
                % Create scalar symbol
                obj.identifier = x;
                obj.maximaInstance = maxima;
                obj.isSymbol = true;
            end
            
            % Validate and process sym_type assumptions
            if ~isempty(sym_type)
                % Check for "clear"
                if any(strcmp(sym_type, "clear"))
                    if length(sym_type) > 1
                        error('sym_type "clear" cannot be combined with other assumptions.');
                    end
                    % Clear assumptions (for future implementation)
                    % Currently no-op as variables are real by default
                else
                    % Validate sym_type values
                    validTypes =["real", "positive", "integer", "rational", "expr"];
                    if any(~(ismember(sym_type, validTypes)))
                        error('Invalid sym_type input. Must be a string or cell array of strings with values: %s, or "clear".', strjoin(validTypes, ', '));
                    end
                    % Apply assumptions (for future implementation)
                    % Currently no-op as variables are real by default
                    % Future: call Maxima's assume() function with appropriate properties
                end
            end
        end

        function s = char(obj)
            if isempty(obj)
                s = '';
            elseif isscalar_matlab(obj)
                s = obj.getDisplayString();
            else
                % For arrays, create cell array of strings
                s = arrayfun(@(x) x.getDisplayString(), obj, 'UniformOutput', false);
            end
        end

        function s = string(obj)
            if isempty(obj)
                s = "";
            elseif isscalar_matlab(obj)
                s = string(obj.getDisplayString());
            else
                % For arrays, create string array
                s = arrayfun(@(x) string(x.getDisplayString()), obj);
            end
        end

        function disp(obj)
            if isempty(obj)
                fprintf('empty %dx%d symbolic array\n\n', size(obj, 1), size(obj, 2));
            elseif isscalar_matlab(obj)
                disp(obj.getDisplayString());
            else
                % For arrays, display size and elements
                fprintf('%dx%d symbolic array:\n\n', size(obj, 1), size(obj, 2));
                for i = 1:size(obj, 1)
                    for j = 1:size(obj, 2)
                        fprintf('  (%d,%d): %s\n', i, j, obj(i,j).getDisplayString());
                    end
                end
            end
        end

        function s = ccode(obj)
            if ~isscalar_matlab(obj)
                s = arrayfun(@ccode, obj, UniformOutput=false);
                return
            end

            cmd = sprintf('gentran(eval(scanmap(''float,%s)))$', obj.identifier);
            obj.validateMaximaInstance;
            [~, extraLines] = obj.maximaInstance.sendAndParse(cmd);
            s = strjoin(extraLines);
        end            

        % Operator Overloads
        function y = plus(a, b)
            % Element-wise addition
            if isempty(a)
                y = b;
                return
            end
            if isempty(b)
                y = a;
                return
            end

            if ~isscalar_matlab(a)
                a = msym.matrix(a);
            end
            if ~isscalar_matlab(b)
                b = msym.matrix(b);
            end

            if isa(a, 'msym') && a.isMatrix && isa(b, 'msym') && b.isMatrix
                if a.matrixRows ~= b.matrixRows || a.matrixCols ~= b.matrixCols
                    error('Matrix dimensions mismatch for addition: [%d x %d] + [%d x %d].', a.matrixRows, a.matrixCols, b.matrixRows, b.matrixCols);
                end
            end
            y = msym.binaryOp(a, b, '+');
        end

        function y = minus(a, b)
            if isempty(a)
                y= b;
                return
            end
            if isempty(b)
                y= a;
                return
            end

            if ~isscalar_matlab(a)
                a = msym.matrix(a);
            end
            if ~isscalar_matlab(b)
                b = msym.matrix(b);
            end
            if isa(a, 'msym') && a.isMatrix && isa(b, 'msym') && b.isMatrix
                if a.matrixRows ~= b.matrixRows || a.matrixCols ~= b.matrixCols
                    error('Matrix dimensions mismatch for subtraction: [%d x %d] - [%d x %d].', a.matrixRows, a.matrixCols, b.matrixRows, b.matrixCols);
                end
            end
            y = msym.binaryOp(a, b, '-');
        end

        function y = times(a, b)
            % Element-wise multiplication .*
            if isempty(a)
                if isscalar(b) || isempty(b)
                    y = a;
                    return
                else
                    error('Trying to multiply non scalar with empty.')
                end
            end
            if isempty(b) 
                if isscalar(a) || isempty(a)
                    y = b;
                    return
                else
                    error('Trying to multiply non scalar with empty.')
                end
            end

            if ~isscalar_matlab(a)
                a = msym.matrix(a);
            end
            if ~isscalar_matlab(b)
                b = msym.matrix(b);
            end

            if isa(a, 'msym') && a.isMatrix && isa(b, 'msym') && b.isMatrix
                if a.matrixRows ~= b.matrixRows || a.matrixCols ~= b.matrixCols
                    error('Matrix dimensions mismatch for element-wise multiplication: [%d x %d] .* [%d x %d].', a.matrixRows, a.matrixCols, b.matrixRows, b.matrixCols);
                end
            end
            y = msym.binaryOp(a, b, '*');
        end

        function y = mtimes(a, b)
            % Matrix multiplication uses dot operator in Maxima (a . b)
            % Check if both operands are matrices
            if isempty(a)
                if isscalar(b) || isempty(b)
                    y = a;
                    return
                else
                    error('Trying to multiply non scalar with empty.')
                end
            end
            if isempty(b) 
                if isscalar(a) || isempty(a)
                    y = b;
                    return
                else
                    error('Trying to multiply non scalar with empty.')
                end
            end

            if ~isscalar_matlab(a)
                a = msym.matrix(a);
            end
            if ~isscalar_matlab(b)
                b = msym.matrix(b);
            end

            aIsMatrix = isa(a, 'msym') && a.isMatrix;
            bIsMatrix = isa(b, 'msym') && b.isMatrix;
            
            if aIsMatrix && bIsMatrix
                % Matrix-matrix multiplication uses "." in Maxima
                if a.matrixCols ~= b.matrixRows
                    error('Matrix dimensions mismatch for multiplication: [%d x %d] * [%d x %d].', a.matrixRows, a.matrixCols, b.matrixRows, b.matrixCols);
                end
                y = msym.binaryOp(a, b, '.');
                % Result is a matrix with dimensions from outer dimensions
                y.isMatrix = true;
                y.matrixRows = a.matrixRows;
                y.matrixCols = b.matrixCols;

                % automatic back conversion to maxima scalar
                if y.matrixRows==1 && y.matrixCols==1
                    y = y(1, 1);
                end
            else
                % Scalar multiplication uses "*"
                y = msym.binaryOp(a, b, '*');
            end
        end

        function y = rdivide(a, b)
            % Element-wise right division ./
            if isempty(a)
                if isscalar(b) || isempty(b)
                    y = a;
                    return
                else
                    error('Trying to divide empty by non scalar.')
                end
            end
            if isempty(b) 
                error('Trying to divide by empty.')
            end

            if ~isscalar_matlab(a)
                a = msym.matrix(a);
            end
            if ~isscalar_matlab(b)
                b = msym.matrix(b);
            end
            
            if isa(a, 'msym') && a.isMatrix && isa(b, 'msym') && b.isMatrix
                if a.matrixRows ~= b.matrixRows || a.matrixCols ~= b.matrixCols
                    error('Matrix dimensions mismatch for element-wise multiplication: [%d x %d] .* [%d x %d].', a.matrixRows, a.matrixCols, b.matrixRows, b.matrixCols);
                end
            end
            y = msym.binaryOp(a, b, '/');
        end

        function y = mrdivide(a, b)
            % Matrix right division a / b
            if isempty(a)
                if isscalar(b) || isempty(b)
                    y = a;
                    return
                else
                    error('Trying to divide empty by non scalar.')
                end
            end
            if isempty(b) 
                error('Trying to divide by empty.')
            end

            if ~isscalar_matlab(a)
                a = msym.matrix(a);
            end
            if ~isscalar_matlab(b)
                b = msym.matrix(b);
            end

            if isa(b, 'msym') && b.isMatrix
                error('Matrix right division is not supported.');
            end
            y = msym.binaryOp(a, b, '/');
        end

        function y = power(a, b)
            % Element-wise power .^
            if isempty(a)
                if isscalar(b)
                    y = a;
                    return
                else
                    error('Trying to raise empty symbol to non scalar or empty power.')
                end
            end
            if isempty(b)
                if isscalar(a)
                    y = b;
                    return
                else
                    error('Trying to raise non scalar or empty symbol to empty power.')
                end
            end

            if ~isscalar_matlab(a)
                a = msym.matrix(a);
            end
            if ~isscalar_matlab(b)
                b = msym.matrix(b);
            end

            if isa(a, 'msym') && a.isMatrix && isa(b, 'msym') && b.isMatrix
                error('Base and exponent cannot both be a matrix.');
            end
            y = msym.binaryOp(a, b, '^');
        end

        function y = mpower(a, b)
            % Matrix power a ^ b
            if isempty(a) || isempty(b)
                error('Trying to raise empty symbol to power or to empty power.')
            end

            if ~isscalar_matlab(a)
                a = msym.matrix(a);
            end
            
            if (isa(b, 'msym') && b.isMatrix) || ~isscalar_matlab(b)
                error('Exponent must not be a matrix.');
            end
            if isa(a, 'msym') && a.isMatrix
                if a.matrixRows ~= a.matrixCols
                    error('Matrix exponentiation requires a square matrix.');
                end
                y = msym.binaryOp(a, b, '^^');
            else
                y = msym.binaryOp(a, b, '^');
            end
        end

        function y = uminus(a)
            if isempty(a)
                y = a;
                return
            end
            if ~isscalar_matlab(a)
                a = msym.matrix(a);
            end

            y = msym.unaryOp(a, '-');
        end

        function y = uplus(a)
            y = a;
        end

        function y = transpose(obj)
            % Array transpose operator .'
            % For symbolic matrix, calls matrixTranspose
            if ~isscalar_matlab(obj)
                y = builtin('transpose', obj);
                return
            end

            if obj.isMatrix
                % Transpose the Maxima matrix
                y = obj.matrixTranspose();
            else
                % Scalar transpose is identity
                y = obj;
            end
        end

        function y = ctranspose(obj)
            % Conjugate transpose operator '
            % For sym (real symbolic expressions), same as transpose
            % For sym matrix, calls matrixTranspose
            if ~isscalar_matlab(obj)
                y = builtin('transpose', msym.funcOp('conjugate', obj));
                return
            end
            
            y = msym.funcOp('conjugate', obj);
            if obj.isMatrix
                % Transpose the Maxima matrix
                y = y.matrixTranspose();
            end
        end

        % Common math functions
        function y = sin(x)
            if ~isscalar_matlab(x)
                y = arrayfun(@sin, x);
                return
            end
            y = msym.funcOp('sin', x);
        end

        function y = cos(x)
            if ~isscalar_matlab(x)
                y = arrayfun(@cos, x);
                return
            end
            y = msym.funcOp('cos', x);
        end

        function y = tan(x)
            if ~isscalar_matlab(x)
                y = arrayfun(@tan, x);
                return
            end
            y = msym.funcOp('tan', x);
        end

        function y = asin(x)
            if ~isscalar_matlab(x)
                y = arrayfun(@asin, x);
                return
            end
            y = msym.funcOp('asin', x);
        end

        function y = acos(x)
            if ~isscalar_matlab(x)
                y = arrayfun(@acos, x);
                return
            end
            y = msym.funcOp('acos', x);
        end

        function y = atan(x)
            if ~isscalar_matlab(x)
                y = arrayfun(@atan, x);
                return
            end
            y = msym.funcOp('atan', x);
        end

        function y = exp(x)
            if ~isscalar_matlab(x)
                y = arrayfun(@exp, x);
                return
            end
            y = msym.funcOp('exp', x);
        end

        function y = expm(x)
            % Matrix exponential or scalar exponential
            % If x is a matrix, computes matrix exponential using matrixexp
            % Otherwise, computes scalar exponential
            if ~isscalar_matlab(x)
                x = msym.matrix(x);
            end

            if isa(x, 'msym') && x.isMatrix
                x.validateMaximaInstance();
                if x.matrixRows ~= x.matrixCols
                    error('Matrix exponential requires a square matrix.');
                end
                cmd = ['matrixexp(', x.identifier, ')'];
                id = x.maximaInstance.sendNoWait(cmd);
                y = msym.fromId(id, x.maximaInstance, x.getDimensions);
            else
                y = msym.funcOp('exp', x);
            end
        end

        function y = log(x)
            if ~isscalar_matlab(x)
                y = arrayfun(@log, x);
                return
            end        
            y = msym.funcOp('log', x);
        end

        function y = sqrt(x)
            if ~isscalar_matlab(x)
                y = arrayfun(@sqrt, x);
                return
            end
            y = msym.funcOp('sqrt', x);
        end

        function y = abs(x)
            if ~isscalar_matlab(x)
                y = arrayfun(@abs, x);
                return
            end
            y = msym.funcOp('abs', x);
        end

        function y = ratsimp(x)
            if ~isscalar_matlab(x)
                y = arrayfun(@ratsimp, x);
                return
            end
            y = msym.funcOp('ratsimp', x);
        end

        function y = trigsimp(x)
            if ~isscalar_matlab(x)
                y = arrayfun(@trigsimp, x);
                return
            end
            y = msym.funcOp('trigsimp', x);
        end

        function y = simplify(x)
            % Calls trigsimp(ratsimp(x)) in Maxima
            if ~isscalar_matlab(x)
                y = arrayfun(@simplify, x);
                return;
            end
            exprX = msym.toMaximaStr(x);
            
            % Get and validate instance from operand
            if isa(x, 'msym')
                x.validateMaximaInstance();
                maxima = x.maximaInstance;
            else
                maxima = MaximaInterface.getInstance();
            end
            
            cmd = ['trigsimp(ratsimp(', exprX, '))']; 
            id = maxima.sendNoWait(cmd);
            y = msym.fromId(id, maxima, x.getDimensions());
        %     s0: trigsimp(fullratsubst([eps_rot^3=0], x)),
        %     s1: ratsimp(trigreduce(s0)),
        %     s2: ratsimp(s0),
        %     if exp_size(s1)>exp_size(s2) then s3: s2 else s3:s1,
        % /*    s4: ratsimp(ev(s3, float)),*/
        %     s4: fullratsimp(s3),
        %     if exp_size(s3)>exp_size(s4) then s4 else s3
        end

        function y = diff(x, var, n)
            arguments
                x msym
                var msym
                n (1,1) double {mustBeInteger, mustBePositive} = 1
            end
            % Differentiate expression with respect to a variable
            % y = diff(expr, var)      - first derivative
            % y = diff(expr, var, n)   - nth derivative
            if isempty(x)
                y = x;
                return
            end
            if ~isscalar_matlab(x)
                y = arrayfun(@(xi)diff(xi,var,n), x);
                return
            end
            x.validateMaximaInstance();
            maxima = x.maximaInstance;

            if isempty(var)
                error('Var cannot be empty.');
            end
            if ~isscalar_matlab(var) || var.isMatrix
                error('Can only differentiate with respect to a single variable.');
            end
            if var.maximaInstance ~= maxima
                error('All msym arguments must belong to the same Maxima instance.');
            end
            
            res = sendAndParse(maxima, ['every(lambda([e], (atom(e) and not numberp(e)) or (not atom(e) and not ratp(e))),[', var.identifier, '])']);
            if isempty(res) || ~strcmpi(res{1}, 'true')
                error('All var must be scalar symbol or function like diff.');
            end

            if n==1
                cmd = sprintf('diff(%s,%s)', x.identifier, var.identifier);
            else
                cmd = sprintf('diff(%s,%s,%d)', x.identifier, var.identifier, n);
            end
            id = maxima.sendNoWait(cmd);
            y = msym.fromId(id, maxima, x.getDimensions);
        end

        function y = jacobian(exprs, vars)
            arguments
                exprs msym
                vars msym
            end
            % Compute Jacobian matrix (matrix of partial derivatives)
            % y = jacobian(exprs, vars)   - Jacobian of expressions w.r.t. variables
            % exprs: symbolic expressions
            % vars:  symbols
            if isempty(exprs)
                y = exprs;
                return
            end
            if isscalar_matlab(exprs)
                if exprs.isMatrix
                    exprsStr = ['flatten(args(', exprs.identifier, '))'];
                    numRows = exprs.matrixRows*exprs.matrixCols;
                else
                    exprsStr = ['[', exprs.identifier, ']'];
                    numRows = 1;
                end
                exprs.validateMaximaInstance();
                maxima = exprs.maximaInstance;
            else
                if any(arrayfun(@(e) e.isMatrix, exprs))
                    error('All elements in exprs must be scalar syms.');
                end
                exprs(1).validateMaximaInstance();
                maxima = exprs(1).maximaInstance;
                for i = 2:numel(exprs)
                    if exprs(i).maximaInstance ~= maxima
                        error('All sym arguments must belong to the same Maxima instance.');
                    end
                end
                exprsStr = ['[', strjoin(arrayfun(@(e) e.identifier, exprs(:), 'UniformOutput', false), ', '), ']'];
                numRows = numel(exprs);
            end

            if isempty(vars)
                error('vars cannot be empty.');
            end
            if isscalar_matlab(vars)
                if vars.maximaInstance ~= maxima
                    error('All sym arguments must belong to the same Maxima instance.');
                end
                if vars.isMatrix
                    cmd = sprintf('flatten(args(%s))', vars.identifier);
                    varsId = maxima.sendNoWait(cmd);
                    numCols = vars.matrixRows*vars.matrixCols;
                else
                    cmd = ['[', vars.identifier, ']'];
                    varsId = maxima.sendNoWait(cmd);
                    numCols = 1;
                end
            else
                for i = 1:numel(vars)
                    if vars(i).maximaInstance ~= maxima
                        error('All sym arguments must belong to the same Maxima instance.');
                    end
                end
                if any(arrayfun(@(v) v.isMatrix, vars))
                    error('All elements in vars must be scalar syms.');
                end
                
                cmd = ['[', strjoin(arrayfun(@(v) v.identifier, vars(:), 'UniformOutput', false), ', '), ']'];
                varsId = maxima.sendNoWait(cmd);
                numCols = numel(vars);
            end
            res = sendAndParse(maxima, ['every(lambda([e], (atom(e) and not numberp(e)) or (not atom(e) and not ratp(e))),', varsId, ')']);
            if isempty(res) || ~strcmpi(res{1}, 'true')
                error('All elements in vars must be scalar symbols or funcitions like diff.');
            end
            
            if numRows == 1 && numCols == 1
                if vars.isMatrix
                    cmd = ['diff(', exprs.identifier, ', ', vars.identifier, '[1,1])'];
                else
                    cmd = ['diff(', exprs.identifier, ', ', vars.identifier, ')'];
                end
                numRows = 0; numCols = 0; % convention for scalars
            else
                cmd = ['jacobian(', exprsStr, ', ', varsId, ')'];
            end
            id = maxima.sendNoWait(cmd);

            y = msym.fromId(id, maxima, [numRows, numCols]);
        end

        function y = subs(x, old, new, fullrat)
            arguments
                x msym
                old
                new
                fullrat (1,1) logical = false
            end
            % Substitute symbols/expressions in x with new values
            % MATLAB API: subs(expr, old, new)

            if isempty(x) || isempty(old) || isempty(new)
                y = x;
                return;
            end

            if ~isscalar_matlab(x)
                y = arrayfun(@(xi) subs(xi, old, new), x);
                return;
            end

            x.validateMaximaInstance();
            maxima = x.maximaInstance;

            oldSubsStr= msym.toMaximaSubstArg(old, maxima);
            newSubsStr= msym.toMaximaSubstArg(new, maxima);

            if isscalar(newSubsStr)
                subsList = cellfun(@(o) [o, '=', newSubsStr{1}], oldSubsStr, 'UniformOutput', false);
            elseif length(oldSubsStr) == length(newSubsStr)
                subsList = cellfun(@(o,n) [o, '=', n], oldSubsStr, newSubsStr, 'UniformOutput', false);
            else
                error('subs: old and new must have the same number of elements.');
            end
            if fullrat
                % Use fullratsimp after substitution
                cmd = ['fullratsubst([', strjoin(subsList, ', '), '], ', x.identifier, ')'];
            else
                cmd = ['subst([', strjoin(subsList, ', '), '], ', x.identifier, ')'];
            end
            id = maxima.sendNoWait(cmd);

            y = msym.fromId(id, maxima, x.getDimensions());
        end

        function vars = symvar(x)
            arguments
                x msym
            end
            % Find symbolic variables in expression
            % MATLAB API: symvar(expr) returns symbolic variables in alphabetical order
            % This implementation returns a string array of variable names
            if isempty(x)
                vars = string.empty(0, 1);
                return
            end
            if ~isscalar_matlab(x)
                for i = 2:numel(x)
                    if x(i).maximaInstance ~= x(1).maximaInstance
                        error('All sym arguments must belong to the same Maxima instance.');
                    end
                end
            end
            
            % Use Maxima's listofvars() to get list of variables
            idStrs = arrayfun(@(e) e.identifier, x(:), 'UniformOutput', false);
            idList = strjoin(idStrs, ', ');
            cmd = ['listofvars([', idList, '])'];
            x(1).validateMaximaInstance();
            result = x(1).maximaInstance.sendAndParse([cmd, ';']);
            
            if isempty(result)
                vars = string.empty(0, 1);
                return;
            end
            
            % Should always be a single line result
            resultStr = result{1};

            % Parse the result - Maxima returns a list like [x, y, z]
            % Remove leading/trailing brackets and whitespace
            resultStr = strtrim(resultStr);
            if startsWith(resultStr, '[')
                resultStr = resultStr(2:end);
            end
            if endsWith(resultStr, ']')
                resultStr = resultStr(1:end-1);
            end
            
            % Split by comma and clean up
            if isempty(strtrim(resultStr))
                vars = string.empty(0, 1);
            else
                varList = strsplit(resultStr, ',');
                vars = string(strtrim(varList));
                % Sort alphabetically to match MATLAB behavior
                vars = sort(vars);
            end
        end

        function tf = contains(x, varargin)
            % Check whether expression contains one or more subexpressions/symbols
            % Uses Maxima's freeof() and returns the negation

            if isempty(varargin)
                error('contains: at least one search expression must be supplied.');
            end

            if isempty(x)
                tf = false;
                return
            end

            if ~isscalar_matlab(x)
                tf = any(arrayfun(@(xi) contains(xi, varargin{:}), x));
                return
            end

            x.validateMaximaInstance();
            maxima = x.maximaInstance;

            % Build list of search terms for freeof
            itemsStrs = {};
            for k = 1:numel(varargin)
                s = msym.toMaximaStr(varargin{k}, false); % don't check for valid symbol strings because operands are also allowed
                if iscell(s)
                    itemsStrs = [itemsStrs; s(:)]; %#ok<AGROW>
                else
                    itemsStrs = [itemsStrs; {s}]; %#ok<AGROW>
                end
            end

            if isempty(itemsStrs)
                error('contains: at least one search expression must be supplied.');
            end

            cmd = ['freeof(', strjoin(itemsStrs(:)', ', '), ', ', x.identifier, ')'];
            res = maxima.sendAndParse([cmd, ';']);
            if isempty(res)
                error('contains: failed to parse freeof result.');
            end
            if strcmpi(res{1}, 'true')
                tf = false;
            elseif strcmpi(res{1}, 'false')
                tf = true;
            else
                error('contains: unexpected freeof result "%s".', res{1});
            end
        end

        function assume(~, varargin) %#ok<INUSD>
            % Declare assumptions about variables (STUB)
            % assume(var, property) - e.g., assume(x, 'real')
            % 
            % Note: Maxima assumes variables are real by default, which is the 
            % default behavior in CADynM. This method is a stub for future enhancement
            % and compatibility with MATLAB Symbolic Toolbox API.
            
            % For now, this is a no-op since all variables are real by default in Maxima
            % Future enhancement can delegate to Maxima's assume() function if needed
        end

        function assumeAlso(~, varargin) %#ok<INUSD>
            % Add additional assumptions about variables (STUB)
            % assumeAlso(var, property) - add property without clearing previous assumptions
            %
            % Note: This is a stub for future enhancement and compatibility with 
            % MATLAB Symbolic Toolbox API.
            
            % For now, this is a no-op since all variables are real by default in Maxima
        end

        function [opt_obj, tnames, tvalues] = optimize(obj)
            is_list = false;
            if ~isscalar_matlab(obj)
                objStrs = arrayfun(@(x) x.identifier, obj, 'UniformOutput', false);
                objStr = ['[', strjoin(objStrs, ', '), ']'];
                maxima = obj(1).maximaInstance;
                obj(1).validateMaximaInstance();
                for i = 2:numel(obj)
                    if obj(i).maximaInstance ~= maxima
                        error('All sym arguments must belong to the same Maxima instance.');
                    end
                end
                is_list = true;
            else
                objStr = obj.identifier;
                maxima = obj.maximaInstance;
                obj.validateMaximaInstance();
            end

            cmd = ['optimizeWithIntermediates(', objStr, ', temp)'];
            id = maxima.sendNoWait(cmd);

            cmd = [id '[1]'];
            res = maxima.sendAndParse(cmd);
            if isempty(res)
                error('Failed to parse response 1 from optimizeWithIntermediates.');
            end
            n_temps = str2double(res{1});

            if n_temps == 0
                tnames = string.empty(0, 1);
                tvalues = msym.empty(0, 1);
            else
                cmd = [id, '[2]'];
                idTempNames = maxima.sendNoWait(cmd);
                tnames = msym.unpackList(idTempNames, maxima);

                cmd = [id, '[3]'];
                idTempExprs = maxima.sendNoWait(cmd);
                tvalues = msym.unpackList(idTempExprs, maxima);
            end

            cmd = [id, '[4]'];
            idOpts = maxima.sendNoWait(cmd);
            if ~is_list
                opt_obj = msym.fromId(idOpts, maxima, obj.getDimensions());
            else
                opt_obj = msym.unpackList(idOpts, maxima);
            end
        end

        function validateMaximaInstance(obj)
            % Check if the cached Maxima instance is still valid
            % This should never be called obj being an array
            if isempty(obj.maximaInstance) || ~isvalid(obj.maximaInstance)
                error('sym: The Maxima interface instance has been deleted. Expression is invalid.');
            end
        end

        function newObj = copy(obj)
            % Create a copy of the sym
            % For matrices, creates a new matrix in Maxima using copymatrix
            % For other expressions, creates a shallow copy with the same identifier
            if ~isscalar_matlab(obj)
                newObj = arrayfun(@msym.copy, obj);
                return
            end
            
            obj.validateMaximaInstance();
            
            if obj.isMatrix
                % For matrices, use Maxima's copymatrix to create a true copy
                cmd = ['copymatrix(', obj.identifier, ')'];
                id = obj.maximaInstance.sendNoWait(cmd);
                newObj = msym.fromId(id, obj.maximaInstance, obj.getDimensions());
            else
                % For non-matrices, shallow copy is sufficient
                newObj = msym(obj.identifier, "expr", obj.maximaInstance);
                % Copy type flags
                newObj.isNumber = obj.isNumber;
                newObj.isSymbol = obj.isSymbol;
                newObj.isConst = obj.isConst;
                newObj.isExpression = obj.isExpression;
            end
            
            % Note: cachedOutput is intentionally not copied to allow fresh evaluation
        end

        function y = matrixElement(obj, row, col)
            % Get a matrix element
            % row: row index (1-based)
            % col: column index (1-based)
            arguments
                obj msym
                row (1,1) double {mustBeInteger, mustBePositive}
                col (1,1) double {mustBeInteger, mustBePositive}
            end
            
            if ~isscalar_matlab(obj) || ~obj.isMatrix
                error('matrixElement can only be called on matrix sym.');
            end
            
            if row < 1 || row > obj.matrixRows || col < 1 || col > obj.matrixCols
                error('Matrix index (%d, %d) out of bounds [1:%d, 1:%d].', row, col, obj.matrixRows, obj.matrixCols);
            end
            
            obj.validateMaximaInstance();
            cmd = [obj.identifier, '[', num2str(row), ',', num2str(col), ']'];
            id = obj.maximaInstance.sendNoWait(cmd);
            y = msym.fromId(id, obj.maximaInstance, [0, 0]);
        end

        function y = matrixTranspose(obj)
            arguments
                obj msym
            end 
            % Transpose a matrix
            
            if ~isscalar_matlab(obj) || ~obj.isMatrix
                error('matrixTranspose can only be called on matrix sym.');
            end
            
            obj.validateMaximaInstance();
            cmd = ['transpose(', obj.identifier, ')'];
            id = obj.maximaInstance.sendNoWait(cmd);
            
            y = msym.fromId(id, obj.maximaInstance, [obj.matrixCols, obj.matrixRows]);
        end

        function y = matrixDeterminant(obj)
            arguments
                obj msym
            end
            % Calculate matrix determinant
            
            if ~isscalar_matlab(obj) || ~obj.isMatrix
                error('matrixDeterminant can only be called on matrix sym.');
            end
            
            if obj.matrixRows ~= obj.matrixCols
                error('Determinant requires a square matrix.');
            end
            
            obj.validateMaximaInstance();
            cmd = ['determinant(', obj.identifier, ')'];
            id = obj.maximaInstance.sendNoWait(cmd);
            y = msym.fromId(id, obj.maximaInstance, [0, 0]);
        end

        function y = matrixInverse(obj)
            arguments
                obj msym
            end
            
            if ~isscalar_matlab(obj) || ~obj.isMatrix
                error('matrixInverse can only be called on matrix sym.');
            end
            
            if obj.matrixRows ~= obj.matrixCols
                error('Inverse requires a square matrix.');
            end
            
            obj.validateMaximaInstance();
            cmd = ['invert(', obj.identifier, ')'];
            id = obj.maximaInstance.sendNoWait(cmd);
            
            y = msym.fromId(id, obj.maximaInstance, obj.getDimensions());
        end

        function y = diag(obj)
            % Extract diagonal from matrix or create diagonal matrix from vector
            %
            % For a matrix: y = diag(M) extracts the main diagonal as a column vector
            % For a vector: y = diag(v) creates a diagonal matrix with v on the diagonal
            
            if isempty(obj)
                y = obj;
                return
            end

            if isscalar_matlab(obj)
                obj.validateMaximaInstance();
                if obj.isMatrix
                    if obj.matrixRows == 1  || obj.matrixCols == 1
                        if obj.matrixRows == 1
                            vecLen = obj.matrixCols;
                            cmd = sprintf('apply(''matrix, makelist(makelist(if i=j then %s[1,i] else 0, j,1,%d), i,1,%d))', obj.identifier, vecLen, vecLen);
                        else
                            vecLen = obj.matrixRows;
                            cmd = sprintf('apply(''matrix, makelist(makelist(if i=j then %s[i,1] else 0, j,1,%d), i,1,%d))', obj.identifier, vecLen, vecLen);
                        end
                        id = obj.maximaInstance.sendNoWait(cmd);
                        y = msym.fromId(id, obj.maximaInstance, [vecLen, vecLen]);
                    else
                        % Input is a matrix - extract diagonal
                        minDim = min(obj.matrixRows, obj.matrixCols);
                        cmd = sprintf('transpose(matrix(makelist(%s[i,i], i, 1, %d)))', obj.identifier, minDim);
                        id = obj.maximaInstance.sendNoWait(cmd);
                        y = msym.fromId(id, obj.maximaInstance, [minDim, 1]);
                    end
                else
                    y = obj; % Scalar is its own diagonal
                end
            else
                if any(arrayfun(@(e) e.isMatrix, obj))
                    error('All elements in the array must be non-matrix syms.');
                end
                obj(1).validateMaximaInstance();
                maxima = obj(1).maximaInstance;
                for i = 2:numel(obj)
                    if obj(i).maximaInstance ~= maxima
                        error('All sym arguments must belong to the same Maxima instance.');
                    end
                end
                
                if isvector_matlab(obj)
                    listStr = strjoin(arrayfun(@(e) e.identifier, obj(:), 'UniformOutput', false), ', ');
                    vecLen = numel(obj);
                    cmd = sprintf('[%s]', listStr);
                    id = maxima.sendNoWait(cmd);
                    cmd = sprintf('apply(''matrix,makelist(makelist(if i=j then %s[i] else 0, j,1,%d), i,1,%d))', id, vecLen, vecLen);
                    id = maxima.sendNoWait(cmd);
                    y = msym.fromId(id, maxima, [vecLen, vecLen]);
                else
                    minDim = min(size(obj, 1), size(obj, 2));
                    listStr = cell(1, minDim);
                    for k = 1:minDim
                        listStr{k} = obj(k, k).identifier;
                    end
                    listStr = strjoin(listStr, ', ');
                    cmd = sprintf('transpose(matrix([%s]))', listStr);
                    id = maxima.sendNoWait(cmd);
                    y = msym.fromId(id, maxima, [minDim, 1]);
                end
            end
        end

        function depends(var, varargin)
            % Declare dependencies between variables
            % var: sym that is a symbol (dependent variable)
            % varargin: list, array or cell array of symbols (independent variables)
            
            if isempty(varargin)
                error('At least one independent variable must be supplied.')
            end

            if ~isscalar_matlab(var)
                arrayfun(@(v) v.depends(varargin{:}), var);
                return;
            end

            % Validate var is a symbol
            if ~var.isSymbol
                error('First argument to depends must be a symbol.');
            end
            
            if isscalar(varargin)
                args = varargin{1};
            else
                args = varargin;
            end
            % Convert indepVars to cell array if needed
            if ~iscell(args)
                if isscalar_matlab(args)
                    args = {args};
                else
                    args = num2cell(args);
                end
            end
            
            % Validate all independent variables are symbols
            for i = 1:length(args)
                if ~isa(args{i}, 'msym') || ~args{i}.isSymbol
                    error('All elements in in the independent variables must be symbols, found "%s".', string(args{i}));
                end
            end
            
            % Get Maxima instance from var
            var.validateMaximaInstance();
            maxima = var.maximaInstance;
            
            % Build the list of independent variables
            indepStrs = cellfun(@(x) x.identifier, args, 'UniformOutput', false);
            indepList = strjoin(indepStrs, ', ');
            
            % Build and execute the command
            cmd = ['depends(', var.identifier, ', [', indepList, '])'];
            maxima.sendNoWait(cmd);
        end
    
        function gradef(var, indep, dsym)
            % Declare dependencie between variables and the symbol to use
            % as the derivative
            % var: sym that is a symbol (dependent variable)
            % indep: symbol (independent variable)
            % dsym: symbol to use as the derivative
            if ~isscalar_matlab(var)
                arrayfun(@(v) v.gradef(indep, dsym), var);
                return;
            end

            % TODO: allow list of dep and dsym

            % Validate var is a symbol
            if ~var.isSymbol
                error('First argument to gradef must be a symbol.');
            end
            
            % check arguments
            if isa(indep, 'msym')
                if ~isscalar_matlab(indep) || indep.isMatrix || ~indep.isSymbol
                    error('Independent variable must be a scalar symbol, found "%s".', string(indep))
                end
                indepStr = indep.identifier;
            elseif ischar(indep) || (isstring(indep) && isscalar(indep))
                indepStr = char(indep);
                if ~mysm.isValidSymbolName(indepStr)
                    error('"%s" is not a valid symbol name.', indepStr)
                end
            else
                error('Independent variable must be a scalar symbol, found "%s".', string(indep))
            end

            if isa(dsym, 'msym')
                if ~isscalar_matlab(dsym) || dsym.isMatrix || ~dsym.isSymbol
                    error('Derivative symbol must be a scalar symbol, found "%s".', string(dsym))
                end
                dsymStr = dsym.identifier;
            elseif ischar(dsym) || (isstring(dsym) && isscalar(indsymdep))
                dsymStr = char(dsym);
                if ~msym.isValidSymbolName(dsymStr)
                    error('"%s" is not a valid symbol name.', dsymStr)
                end                
            else
                error('Derivative symbol must be a scalar symbol, found "%s".', string(dsym))
            end
            
            % Get Maxima instance from var
            var.validateMaximaInstance();
            maxima = var.maximaInstance;
            
            % Build and execute the command
            cmd = ['gradef(', var.identifier, ', ', indepStr, ', ', dsymStr, ')'];
            maxima.sendNoWait(cmd);
        end
    
        function tf = isscalar(obj)
            % Check if sym is a true scalar (not a matrix and not an array)
            % Returns false if:
            %   - obj is a MATLAB array (more than one element)
            %   - obj contains a Maxima matrix
            % Returns true only if obj is a single scalar sym (not a matrix)
            if ~isscalar_matlab(obj)
                % MATLAB array or empty
                tf = false;
            elseif obj.isMatrix
                % Single MATLAB object but it's a Maxima matrix
                tf = false;
            else
                % Single MATLAB object and not a Maxima matrix
                tf = true;
            end
        end

        function tf = isvector(obj)
            % Check if sym is a vector (1D array or Maxima vector) or a MATLAB vector
            if isscalar_matlab(obj)
                if obj.isMatrix
                    % Single Maxima matrix - check if it's a vector
                    tf = obj.matrixRows == 1 || obj.matrixCols == 1;
                else
                    tf = true; % Single scalar sym (not a matrix) is considered a vector
                end
            else
                tf = isvector_matlab(obj);
            end
        end

        function varargout = size(obj, varargin)
            % SIZE  Size of msym object
            %   S = SIZE(A) returns the matrix dimensions of A as a row vector
            %   [M,N] = SIZE(A) returns the number of rows and columns separately
            %   SIZE(A,DIM) returns the size in dimension DIM
            
            if isscalar_matlab(obj) && obj.isMatrix
                % Single matrix object
                dims = obj.getDimensions();
            elseif ~isscalar_matlab(obj)
                % Array of msym objects
                dims = builtin('size', obj);
            else
                % Scalar non-matrix object
                dims = [1, 1];
            end
            
            if nargin == 1
                % No dimension specified - return full size vector
                if nargout <= 1
                    varargout{1} = dims;
                else
                    % Multiple output arguments: [m, n] = size(A)
                    for i = 1:nargout
                        if i <= length(dims)
                            varargout{i} = dims(i); %#ok<AGROW>
                        else
                            varargout{i} = 1; %#OK<AGROW>
                        end
                    end
                end
            else
                % Dimension specified
                dim = varargin{1};
                if ~isscalar(dim) || dim < 1 || mod(dim, 1) ~= 0
                    error('Size argument must be a positive integer scalar.');
                end
                
                if dim <= length(dims)
                    varargout{1} = dims(dim);
                else
                    varargout{1} = 1;
                end
            end
        end

        function n = numel(obj, varargin)
            % NUMEL  Number of elements in msym object
            %   N = NUMEL(A) returns the total number of elements in A
            %   NUMEL(A, INDEX1, INDEX2, ...) returns the number of elements in the indexed region
            
            if isscalar_matlab(obj) && obj.isMatrix
                % Single matrix object
                dims = obj.getDimensions();
                n = prod(dims);
            elseif ~isscalar_matlab(obj)
                % Array of msym objects
                n = builtin('numel', obj, varargin{:});
            else
                % Scalar non-matrix object
                n = 1;
            end
        end

        function n = length(obj)
            % LENGTH  Length of msym object
            %   N = LENGTH(A) returns the largest dimension of A
            %   For a vector, returns its length
            %   For a scalar or matrix, returns the larger of the two dimensions
            
            if isscalar_matlab(obj) && obj.isMatrix
                % Single matrix object
                dims = obj.getDimensions();
                n = max(dims);
            elseif ~isscalar_matlab(obj)
                % Array of msym objects
                n = builtin('length', obj);
            else
                % Scalar non-matrix object
                n = 1;
            end
        end

        function varargout = subsref(obj, s)
            % Subscript reference - handle indexing like obj(i,j)
            if isempty(s)
                varargout = {obj};
                return;
            end
            if ~isscalar_matlab(obj)
                varargout = {builtin('subsref', obj, s)};
                return
            end
            
            if isscalar_matlab(obj) && obj.isMatrix && strcmp(s(1).type, '()')
                if length(s(1).subs) > 2
                    error('Cannot assign to matrix beyond 2 subscripts (row, col).');
                end
                [rows, cols, scalar_index] = checksubs(obj, s, false);

                if ~isscalar(rows) || ~isscalar(cols)
                    % compose a command string to create a new matrix from the selected elements
                    if scalar_index
                        dataStr = cell(numel(rows), 1);
                        for i = 1:numel(rows)
                            dataStr{i,1} = [obj.identifier, '[', num2str(rows(i)), ',', num2str(cols(i)), ']'];
                        end
                    else
                        dataStr = cell(numel(rows), numel(cols));
                        for i = 1:numel(rows)
                            for j = 1:numel(cols)
                                dataStr{i,j} = [obj.identifier, '[', num2str(rows(i)), ',', num2str(cols(j)), ']'];
                            end
                        end
                    end
                    % Build Maxima matrix syntax: matrix([row1], [row2], ...)
                    rowStrs = cell(size(dataStr,1), 1);
                    for i = 1:size(dataStr,1)
                        rowStr = strjoin(dataStr(i,:), ', ');
                        rowStrs{i} = ['[', rowStr, ']'];
                    end
                    matrixStr = strjoin(rowStrs, ', ');
                    
                    % Send matrix creation command
                    cmd = ['matrix(', matrixStr, ')'];
                    obj.validateMaximaInstance();
                    id = obj.maximaInstance.sendNoWait(cmd);
                    
                    % Create sym object and set matrix properties
                    if scalar_index
                        y = msym.fromId(id, obj.maximaInstance, [numel(rows), 1]);
                    else
                        y = msym.fromId(id, obj.maximaInstance, [numel(rows), numel(cols)]);
                    end
                else
                    % Get one element
                    y = obj.matrixElement(rows, cols);
                end

                % Chained indexing doesn't make sense for symbols: throw error
                if length(s) > 1
                    error('Chained indexing not supported for symbols.');
                else
                    varargout = {y};
                end
            else
                % Handle dot notation (property access)
                try
                    % try calling expecting outputs
                    nret = max(1,nargout);
                    varargout{1:nret} = builtin('subsref', obj, s);
                catch ME
                    if strcmp(ME.identifier, 'MATLAB:TooManyOutputs')
                        % With this error, the function was never executed,
                        % so we need to do that here
                        builtin('subsref', obj, s);
                    else
                        rethrow(ME)
                    end
                end
            end
        end

        function obj = subsasgn(obj, s, val)
            % Subscript assignment - handle assignment like obj(i,j) = expr
            if isempty(s)
                return;
            end
            if isscalar_matlab(obj) && obj.isMatrix && strcmp(s(1).type, '()')
                if length(s(1).subs) > 2
                    error('Cannot assign to matrix beyond 2 subscripts (row, col).');
                end
                [rows, cols, scalar_index] = checksubs(obj, s);
                
                % prepare data
                if scalar_index
                    dataStr = cell(numel(rows), 1);
                else
                    dataStr = cell(numel(rows), numel(cols));
                end
                if isa(val, 'msym') && isscalar_matlab(val) && val.isMatrix
                    if ~isscalar(val)
                        if scalar_index
                            % TODO: possibly allow to assign row vector
                            if numel(rows) ~= val.matrixRows || 1 ~= val.matrixCols
                                error('Assigned matrix size does not match target submatrix size.');
                            end
                            for i = 1:numel(rows)
                                dataStr{i,1} = [val.identifier, '[', num2str(i), ',', num2str(i), ']'];
                            end
                        else
                            if numel(rows) ~= val.matrixRows || numel(cols) ~= val.matrixCols
                                error('Assigned matrix size does not match target submatrix size.');
                            end
                            for i = 1:numel(rows)
                                for j = 1:numel(cols)
                                    dataStr{i,j} = [val.identifier, '[', num2str(i), ',', num2str(j), ']'];
                                end
                            end
                        end
                    else
                        dataStr = [val.identifier, '[1,1]'];
                    end
                else
                    if ischar(val)
                        val = cellstr(val);
                    end
                    if ~iscell(val)
                        if isscalar_matlab(val)
                            val = {val};
                        else
                            val = num2cell(val);
                        end
                    end
                    if ~isscalar(val)
                        if scalar_index
                            % TODO: possibly allow to assign row vector
                            if size(val, 1) ~= numel(rows) || 1 ~= numel(cols)
                                error('Assigned value size does not match target submatrix size.');
                            end
                            for i = 1:numel(rows)
                                dataStr{i,1} = msym.toMaximaStr(val{i});
                            end
                        else
                            if size(val, 1) ~= numel(rows) || size(val, 2) ~= numel(cols)
                                error('Assigned value size does not match target submatrix size.');
                            end
                            for i = 1:numel(rows)
                                for j = 1:numel(cols)
                                    dataStr{i,j} = msym.toMaximaStr(val{i,j});
                                end
                            end
                        end
                    else
                        dataStr = msym.toMaximaStr(val{1,1});
                    end
                end

                % Create assignment command in Maxima
                obj.validateMaximaInstance();
                if scalar_index
                    for i = 1:numel(rows)
                        if iscell(dataStr)
                            cmd = [obj.identifier, '[', num2str(rows(i)), ',', num2str(cols(i)), '] : ', dataStr{i,1}];
                        else
                            cmd = [obj.identifier, '[', num2str(rows(i)), ',', num2str(cols(i)), '] : ', dataStr];
                        end
                        obj.maximaInstance.sendNoWait(cmd);
                    end
                else
                    for i = 1:numel(rows)
                        for j = 1:numel(cols)
                            if iscell(dataStr)
                                cmd = [obj.identifier, '[', num2str(rows(i)), ',', num2str(cols(j)), '] : ', dataStr{i,j}];
                            else
                                cmd = [obj.identifier, '[', num2str(rows(i)), ',', num2str(cols(j)), '] : ', dataStr];
                            end
                            obj.maximaInstance.sendNoWait(cmd);
                        end
                    end
                end                    
        
                % Clear cached output since matrix changed
                obj.cachedOutput = '';
            else
                obj = builtin('subsasgn', obj, s, val);
            end
        end

        function idx = subsindex(obj)
            % Convert sym to an index for use in subscripting
            % For numeric expressions, convert to a MATLAB number
            if ~isscalar_matlab(obj)
                error('Can only use scalars or maxima matrices as index.')
            end
            if obj.isNumber
                % Try to extract numeric value
                obj.validateMaximaInstance();
                result = obj.maximaInstance.sendAndParse(['float(', obj.identifier, ');']);
                if ischar(result)
                    idx = str2double(result);
                else
                    idx = result;
                end
                % MATLAB indices must be positive integers
                if idx ~= floor(idx) || idx < 1
                    error('Array indices must be positive integers.');
                end
            else
                error('Only numeric syms can be used as array indices.');
            end
        end

        function idx = end(obj, k, n)
            % Support the use of 'end' keyword in subscripting
            % k: dimension being indexed (1 for rows, 2 for columns)
            % n: total number of subscripts
            
            if isscalar_matlab(obj) && obj.isMatrix
                if n ~= 2
                    error('Matrix subscripting currently only supports row and column indices.');
                end
                switch k
                    case 1
                        % Row dimension
                        idx = obj.matrixRows;
                    case 2
                        % Column dimension
                        idx = obj.matrixCols;
                    otherwise
                        error('Matrix subscripting only supports 2 dimensions.');
                end
            else
                % Array of sym objects
                sz = size(obj);
                if k <= length(sz)
                    idx = sz(k);
                else
                    idx = 1;
                end
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
                        if isa(arg, 'msym')
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
                argStrs = cellfun(@msym.toMaximaStr, arg, 'UniformOutput', false);
                exprArgs = strjoin(argStrs, ', ');

                for i = 1:length(arg)
                    if isa(arg{i}, 'msym')
                        if isempty(maxima)
                            % if all instances are the same, vaidating the first is sufficient
                            arg{i}.validateMaximaInstance();
                            maxima = arg{i}.maximaInstance;
                        else
                            if maxima ~= arg{i}.maximaInstance
                                error('All sym arguments must belong to the same MaximaInterface instance.');
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
            dims = maximaMatrixDimensions(id, maxima);

            y = msym.fromId(id, maxima, dims);
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

            obj = msym(id, "const");
        end

        function obj = pi()
            obj = msym.const('pi');
        end

        function obj = e()
            obj = msym.const('e');
        end

        function obj = gamma()
            obj = msym.const('gamma');
        end

        function obj = i()
            obj = msym.const('i');
        end

        function obj = inf()
            obj = msym.const('inf');
        end

        function obj = minf()
            obj = msym.const('minf');
        end

        function obj = phi()
            obj = msym.const('phi');
        end

        function y = matrix(data, maxima)
            % Create a Maxima matrix from MATLAB data
            % data: MATLAB matrix of numbers, strings, or syms
            % maxima: optional MaximaInterface instance
            
            arguments
                data
                maxima = []
            end
            
            if isscalar_matlab(data)
                error('Will not produce a 1x1 matrix. Use sym() constructor instead.');
            end

            if ~isnumeric(data) && ~isstring(data) && ~iscell(data) && ~isa(data, 'msym')
                error('Matrix data must be sym, numeric, string array, cell array, or contain syms.');
            end
            
            % Convert to cell array for uniform handling
            if ~iscell(data)
                data = num2cell(data);
            end
            
            [rows, cols] = size(data);
            
            
            % Convert each element to Maxima representation
            dataStr = cell(rows, cols);
            for i = 1:rows
                for j = 1:cols
                    if isa(data{i,j}, 'msym')
                        if isempty(maxima)
                            data{i,j}.validateMaximaInstance();
                            maxima = data{i,j}.maximaInstance;
                        else
                            if maxima ~= data{i,j}.maximaInstance
                                error('All sym elements must belong to the same MaximaInterface instance.');
                            end
                        end
                    end
                    dataStr{i,j} = msym.toMaximaStr(data{i,j});
                end
            end
            if isempty(maxima)
                maxima = MaximaInterface.getInstance();
            end
            
            % Build Maxima matrix syntax: matrix([row1], [row2], ...)
            rowStrs = cell(rows, 1);
            for i = 1:rows
                rowStr = strjoin(dataStr(i,:), ', ');
                rowStrs{i} = ['[', rowStr, ']'];
            end
            matrixStr = strjoin(rowStrs, ', ');
            
            % Send matrix creation command
            cmd = ['matrix(', matrixStr, ')'];
            id = maxima.sendNoWait(cmd);
            
            % Create sym object and set matrix properties
            y = msym.fromId(id, maxima, size(data));
        end

        function y = matrixSymbolic(name, rows, cols, maxima)
            % Create a symbolic matrix with unspecified entries
            % name: base name for matrix symbols
            % rows: number of rows
            % cols: number of columns
            % maxima: optional MaximaInterface instance
            
            arguments
                name
                rows (1,1) double {mustBeInteger, mustBePositive}
                cols (1,1) double {mustBeInteger, mustBePositive}
                maxima = []
            end
            
            if rows == 1 && cols == 1
                error('Will not produce a 1x1 matrix. Use sym() constructor instead.');
            end

            if isempty(maxima)
                maxima = MaximaInterface.getInstance();
            end
            
            name = char(name);
            
            % Create symbolic matrix using genmatrix
            cmd = ['genmatrix(lambda([i,j], ', name, '[i,j]), ', num2str(rows), ', ', num2str(cols), ')'];
            id = maxima.sendNoWait(cmd);
            
            % Create sym object and set matrix properties
            y = msym.fromId(id, maxima, [rows, cols]);
        end

        function y = zeros(rows, cols, maxima)
            % Create a zero matrix
            % rows: number of rows
            % cols: number of columns
            % maxima: optional MaximaInterface instance
            
            arguments
                rows (1,1) double {mustBeInteger, mustBePositive}
                cols (1,1) double {mustBeInteger, mustBePositive}
                maxima = []
            end
            
            if rows == 1 && cols == 1
                error('Will not produce a 1x1 matrix. Use sym() constructor instead.');
            end

            if isempty(maxima)
                maxima = MaximaInterface.getInstance();
            end
            
            cmd = ['zeromatrix(', num2str(rows), ', ', num2str(cols), ')'];
            id = maxima.sendNoWait(cmd);
            
            y = msym.fromId(id, maxima, [rows, cols]);
        end

        function y = ones(rows, cols, maxima)
            % Create a matrix of ones
            % rows: number of rows
            % cols: number of columns
            % maxima: optional MaximaInterface instance
            
            arguments
                rows (1,1) double {mustBeInteger, mustBePositive}
                cols (1,1) double {mustBeInteger, mustBePositive}
                maxima = []
            end
            
            if rows == 1 && cols == 1
                error('Will not produce a 1x1 matrix. Use sym() constructor instead.');
            end

            if isempty(maxima)
                maxima = MaximaInterface.getInstance();
            end
            
            cmd = ['genmatrix(1, ', num2str(rows), ', ', num2str(cols), ')'];
            id = maxima.sendNoWait(cmd);
            
            y = msym.fromId(id, maxima, [rows, cols]);
        end

        function y = eye(n, maxima)
            % Create an identity matrix
            % n: matrix dimension (n x n)
            % maxima: optional MaximaInterface instance
            
            arguments
                n (1,1) double {mustBeInteger, mustBePositive}
                maxima = []
            end
            
            if n==1
                error('Will not produce a 1x1 matrix. Use sym() constructor instead.');
            end

            if isempty(maxima)
                maxima = MaximaInterface.getInstance();
            end
            
            cmd = ['ident(', num2str(n), ')'];
            id = maxima.sendNoWait(cmd);
            
            y = msym.fromId(id, maxima, [n, n]);
        end

        function dims = maximaMatrixDimensions(id, maxima)
            % Get dimensions of a Maxima matrix given its identifier
            % id: Maxima identifier for the matrix
            % maxima: MaximaInterface instance
            
            result = maxima.sendAndParse(['if matrixp(', id, ') then [length(', id, '), length(first(', id, '))] else [0, 0]']);
            dims = str2num(result{1}); %#ok<ST2NM>
        end

        function y = unpackList(idList, maxima)
            % Helper function to unpack a Maxima list into a MATLAB array of msym objects
            % id: Maxima identifier for the list
            % maxima: MaximaInterface instance
            
            result = maxima.sendAndParse(['length(', idList, ')']);
            n = str2double(result{1});
            
            if n == 0
                y = [];
                return;
            end
            
            y = msym.empty;
            for i = 1:n
                idElem = maxima.sendNoWait([idList '[', num2str(i), ']']);
                dims = msym.maximaMatrixDimensions(idElem, maxima);

                y(i) = msym.fromId(idElem, maxima, dims); %#ok<AGROW>
            end
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

        function dims = getDimensions(obj)
            % just a small helper function
            % but beware: obj.matrixRow==0, obj.matrixCols==0 means scalar
            % in the future obj.matrixRow==-1, obj.matrixCols==-1 may mean
            % empty
            dims = [obj.matrixRows, obj.matrixCols];
        end

        function [rows, cols, scalar_index] = checksubs(obj, s, assigning)
            scalar_index = true;
            if isscalar(s(1).subs)
                if islogical(s(1).subs{1})
                    if ~isequal(size(s(1).subs{1}), size(obj))
                        error('Logical index size does not match matrix size.')
                    end
                    [rows, cols] = ind2sub([obj.matrixRows, obj.matrixCols], find(s(1).subs{1}));
                elseif ischar(s(1).subs{1})
                    if s(1).subs{1}==':'
                        [rows, cols] = meshgrid(1:obj.matrixRows, 1:obj.matrixCols);
                        rows = rows(:);
                        cols = cols(:);
                    else
                        error('Unknown index "%s".', s(1).subs{1})
                    end
                else
                    [rows, cols] = ind2sub([obj.matrixRows, obj.matrixCols], s(1).subs{1});
                end
            else
                scalar_index = false;
                rows = s(1).subs{1};
                cols = s(1).subs{2};

                if ischar(rows)
                    if rows==':'
                        rows = 1:obj.matrixRows;
                    else
                        error('Unknown index "%s".', rows)
                    end
                end
                if ischar(cols) 
                    if cols==':'
                        cols = 1:obj.matrixCols;
                    else
                        error('Unknown index "%s".', cols)
                    end
                end
            end

            % Validate indices
            if any(rows < 1) || any(rows > obj.matrixRows) || any(cols < 1) || any(cols > obj.matrixCols)
                if assigning
                    error('One ore more matrix indeces out of bounds [1:%d, 1:%d]. Automatic expansion currently not supported.', obj.matrixRows, obj.matrixCols);
                else
                    error('One ore more matrix indeces out of bounds [1:%d, 1:%d].', obj.matrixRows, obj.matrixCols);
                end
            end
        end

    end

    methods (Static, Access = private)
        function tf = isValidSymbolName(name)
            tf = ~isempty(regexp(name, '^[A-Za-z][A-Za-z0-9_]*$', 'once'));
        end

        function subsStr = toMaximaSubstArg(val, maxima)
            % Convert substitution arguments to Maxima strings
            % Returns:
            %  subsStr: cell array of Maxima strings

            if isa(val, 'msym')
                if isscalar_matlab(val)
                    if val.maximaInstance ~= maxima
                        error('All msym arguments must belong to the same Maxima instance.');
                    end
                    if val.isMatrix
                        [R,C] = meshgrid(1:val.matrixRows, 1:val.matrixCols);
                        subsStr = arrayfun(@(r,c)sprintf('%s[%d,%d]', val.identifier, r, c), R, C, UniformOutput=false);
                        subsStr = subsStr(:);
                    else
                        subsStr = {val.identifier};
                    end
                    return;
                end

                if any(arrayfun(@(v) v.isMatrix, val))
                    error('Array inputs to subs cannot contain matrix msyms.');
                end
                for i = 1:numel(val)
                    if val(i).maximaInstance ~= maxima
                        error('All msym arguments must belong to the same Maxima instance.');
                    end
                end

                subsStr = arrayfun(@(v) v.identifier, val(:), 'UniformOutput', false);
                return;
            end

            if iscell(val)
                subsStr = cellfun(@(v) msym.toMaximaStr(v), val(:), 'UniformOutput', false);
                return;
            end

            if ischar(val) || isscalar_matlab(val)
                subsStr = {msym.toMaximaStr(val)};
                return;
            end

            if isnumeric(val) || isstring(val)
                subsStr = arrayfun(@(v) msym.toMaximaStr(v), val(:), 'UniformOutput', false);
                return
            end

            error('Unsupported substitution argument type.');
        end

        function y = binaryOp(a, b, op)
            exprA = msym.toMaximaStr(a);
            exprB = msym.toMaximaStr(b);
            
            % Get and validate instance from first sym operand
            maxima = msym.getMaximaFromOperands(a, b);
            
            % no need for parentheses because operands are always symbols or maxima identifiers
            cmd = [exprA, op, exprB];
            id = maxima.sendNoWait(cmd);


            % Handle element-wise matrix operations
            dims = [0, 0];
            if isa(a, 'msym') && a.isMatrix
                dims = a.getDimensions();
            elseif isa(b, 'msym') && b.isMatrix
                dims = b.getDimensions();
            end

            y = msym.fromId(id, maxima, dims);
        end

        function y = unaryOp(a, op)
            exprA = msym.toMaximaStr(a);
            
            % Get and validate instance from operand
            if isa(a, 'msym')
                a.validateMaximaInstance();
                maxima = a.maximaInstance;
            else
                maxima = MaximaInterface.getInstance();
            end
            
            % no need for parentheses because operands are always symbols or maxima identifiers
            cmd = [op, exprA];
            id = maxima.sendNoWait(cmd);

            dims = [0, 0];
            if isa(a, 'msym') && a.isMatrix
                dims = a.getDimensions();
            end
            y = msym.fromId(id, maxima, dims);
        end
        
        function y = funcOp(fname, arg)
            % Apply a function operation to a single argument
            % fname: function name (string)
            % arg: argument (can be sym, number, or string)
            
            if ~isscalar_matlab(arg)
                y = arrayfun(@(e) msym.funcOp(fname, e), arg);
                return;
            end

            exprArg = msym.toMaximaStr(arg);
            
            % Get and validate instance from operand
            if isa(arg, 'msym')
                arg.validateMaximaInstance();
                maxima = arg.maximaInstance;
            else
                maxima = MaximaInterface.getInstance();
            end
            
            cmd = [fname, '(', exprArg, ')'];
            id = maxima.sendNoWait(cmd);
            y = msym.fromId(id, maxima, arg.getDimensions());
        end

        function maxima = getMaximaFromOperands(a, b)
            % Get MaximaInterface instance from operands and validate
            if isa(a, 'msym') && isa(b, 'msym')
                if a.maximaInstance ~= b.maximaInstance
                    error('Both sym operands must belong to the same MaximaInterface instance.');
                end
                % if both instances are the same, vaidating one is sufficient
                a.validateMaximaInstance();
                maxima = a.maximaInstance;
            elseif isa(a, 'msym')
                a.validateMaximaInstance();
                maxima = a.maximaInstance;
            elseif isa(b, 'msym')
                b.validateMaximaInstance();
                maxima = b.maximaInstance;
            else
                maxima = MaximaInterface.getInstance();
            end
        end

        function s = toMaximaStr(x, check_symbol_names)
            arguments
                x
                check_symbol_names (1,1) logical = true
            end
            % Convert a MATLAB value to a Maxima string representation that can be directly used in Maxima commands or as msym identifiers
            %   If x represents a Maxima symbol, a char string is returned
            %   otherwise, a cell of char strings is returned

            if isa(x, 'msym')
                if isscalar_matlab(x)
                    s = x.identifier;
                else
                    s = arrayfun(@(v) v.identifier, x(:), 'UniformOutput', false);
                end
            elseif isnumeric(x)
                if isscalar_matlab(x)
                    s = num2str(x, 16);
                else
                    s = arrayfun(@(v) num2str(v, 16), x(:), 'UniformOutput', false);
                end
            elseif ischar(x) || (isstring(x) && isscalar(x))
                x = char(x);
                if check_symbol_names && ~msym.isValidSymbolName(x)
                    % TODO: also allow numerals
                    error('String operands must be valid symbol names.');
                end
                s = x;
            elseif isstring(x)
                s = arrayfun(@(v) msym.toMaximaStr(v, check_symbol_names), x, 'UniformOutput', false);
            elseif iscell(x)
                s = cellfun(@(v) msym.toMaximaStr(v, check_symbol_names), x, 'UniformOutput', false);
            else
                error('Unsupported operand type.');
            end
        end

        function obj = fromId(id, maxima, dims)
            obj = msym(id, dims, "expr", maxima);
        end
    end
end
