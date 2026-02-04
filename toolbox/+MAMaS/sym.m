classdef sym < handle
    % sym
    % MATLAB class for symbolic Maxima expressions

    properties (SetAccess = private)
        identifier
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
        function obj = sym(x, varargin)
            % Constructor for new symbols

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
                if ~isa(maxima, 'MAMaS.MaximaInterface')
                    error('The maxima argument must be a MAMaS.MaximaInterface instance.');
                end
                nextArgIdx = nextArgIdx + 1;
            end
            
            % Check for extra arguments
            if nextArgIdx <= length(varargin)
                error('Too many input arguments.');
            end
            
            % Set default maxima if not provided
            if isempty(maxima)
                maxima = MAMaS.MaximaInterface.getInstance();
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
            if ~ischar(x) && ~isscalar(x)
                % TODO: error if dim is also provided?
                % TODO: validate sym_type?
                % Delegate to matrix creation
                matObj = MAMaS.sym.matrix(x, maxima);
                % Copy properties from the created matrix object
                obj.identifier = matObj.identifier;
                obj.cachedOutput = matObj.cachedOutput;
                obj.maximaInstance = matObj.maximaInstance;
                obj.isMatrix = matObj.isMatrix;
                obj.matrixRows = matObj.matrixRows;
                obj.matrixCols = matObj.matrixCols;
                obj.isExpression = matObj.isExpression;
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
            if ~MAMaS.sym.isValidSymbolName(x)
                error('Invalid symbol name: "%s".', x);
            end

            % Handle dimension: create matrix or scalar symbol
            if ~isempty(dim)
                % Create symbolic matrix
                matObj = MAMaS.sym.matrixSymbolic(x, dim(1), dim(2), maxima);
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
            if isscalar(obj)
                s = obj.getDisplayString();
            else
                % For arrays, create cell array of strings
                s = arrayfun(@(x) x.getDisplayString(), obj, 'UniformOutput', false);
            end
        end

        function s = string(obj)
            if isscalar(obj)
                s = string(obj.getDisplayString());
            else
                % For arrays, create string array
                s = arrayfun(@(x) string(x.getDisplayString()), obj);
            end
        end

        function disp(obj)
            if isscalar(obj)
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

        % Operator Overloads
        function y = plus(a, b)
            % Element-wise addition
            if ~isscalar(a)
                a = MAMaS.sym.matrix(a);
            end
            if ~isscalar(b)
                b = MAMaS.sym.matrix(b);
            end

            if isa(a, 'MAMaS.sym') && a.isMatrix && isa(b, 'MAMaS.sym') && b.isMatrix
                if a.matrixRows ~= b.matrixRows || a.matrixCols ~= b.matrixCols
                    error('Matrix dimensions mismatch for addition: [%d x %d] + [%d x %d].', a.matrixRows, a.matrixCols, b.matrixRows, b.matrixCols);
                end
            end
            y = MAMaS.sym.binaryOp(a, b, '+');
        end

        function y = minus(a, b)
            if ~isscalar(a)
                a = MAMaS.sym.matrix(a);
            end
            if ~isscalar(b)
                b = MAMaS.sym.matrix(b);
            end
            if isa(a, 'MAMaS.sym') && a.isMatrix && isa(b, 'MAMaS.sym') && b.isMatrix
                if a.matrixRows ~= b.matrixRows || a.matrixCols ~= b.matrixCols
                    error('Matrix dimensions mismatch for subtraction: [%d x %d] - [%d x %d].', a.matrixRows, a.matrixCols, b.matrixRows, b.matrixCols);
                end
            end
            y = MAMaS.sym.binaryOp(a, b, '-');
        end

        function y = times(a, b)
            % Element-wise multiplication .*
            if ~isscalar(a)
                a = MAMaS.sym.matrix(a);
            end
            if ~isscalar(b)
                b = MAMaS.sym.matrix(b);
            end

            if isa(a, 'MAMaS.sym') && a.isMatrix && isa(b, 'MAMaS.sym') && b.isMatrix
                if a.matrixRows ~= b.matrixRows || a.matrixCols ~= b.matrixCols
                    error('Matrix dimensions mismatch for element-wise multiplication: [%d x %d] .* [%d x %d].', a.matrixRows, a.matrixCols, b.matrixRows, b.matrixCols);
                end
            end
            y = MAMaS.sym.binaryOp(a, b, '*');
        end

        function y = mtimes(a, b)
            % Matrix multiplication uses dot operator in Maxima (a . b)
            % Check if both operands are matrices
            if ~isscalar(a)
                a = MAMaS.sym.matrix(a);
            end
            if ~isscalar(b)
                b = MAMaS.sym.matrix(b);
            end

            aIsMatrix = isa(a, 'MAMaS.sym') && a.isMatrix;
            bIsMatrix = isa(b, 'MAMaS.sym') && b.isMatrix;
            
            if aIsMatrix && bIsMatrix
                % Matrix-matrix multiplication uses "." in Maxima
                if a.matrixCols ~= b.matrixRows
                    error('Matrix dimensions mismatch for multiplication: [%d x %d] * [%d x %d].', a.matrixRows, a.matrixCols, b.matrixRows, b.matrixCols);
                end
                y = MAMaS.sym.binaryOp(a, b, '.');
                % Result is a matrix with dimensions from outer dimensions
                y.isMatrix = true;
                y.matrixRows = a.matrixRows;
                y.matrixCols = b.matrixCols;
            else
                % Scalar multiplication uses "*"
                y = MAMaS.sym.binaryOp(a, b, '*');
            end
        end

        function y = rdivide(a, b)
            % Element-wise right division ./
            if ~isscalar(a)
                a = MAMaS.sym.matrix(a);
            end
            if ~isscalar(b)
                b = MAMaS.sym.matrix(b);
            end
            
            if isa(a, 'MAMaS.sym') && a.isMatrix && isa(b, 'MAMaS.sym') && b.isMatrix
                if a.matrixRows ~= b.matrixRows || a.matrixCols ~= b.matrixCols
                    error('Matrix dimensions mismatch for element-wise multiplication: [%d x %d] .* [%d x %d].', a.matrixRows, a.matrixCols, b.matrixRows, b.matrixCols);
                end
            end
            y = MAMaS.sym.binaryOp(a, b, '/');
        end

        function y = mrdivide(a, b)
            % Matrix right division a / b
            if ~isscalar(a)
                a = MAMaS.sym.matrix(a);
            end
            if ~isscalar(b)
                b = MAMaS.sym.matrix(b);
            end

            if isa(b, 'MAMaS.sym') && b.isMatrix
                error('Matrix right division is not supported.');
            end
            y = MAMaS.sym.binaryOp(a, b, '/');
        end

        function y = power(a, b)
            % Element-wise power .^
            if ~isscalar(a)
                a = MAMaS.sym.matrix(a);
            end
            if ~isscalar(b)
                b = MAMaS.sym.matrix(b);
            end

            if isa(a, 'MAMaS.sym') && a.isMatrix && isa(b, 'MAMaS.sym') && b.isMatrix
                error('Base and exponent cannot both be a matrix.');
            end
            y = MAMaS.sym.binaryOp(a, b, '^');
        end

        function y = mpower(a, b)
            % Matrix power a ^ b
            if ~isscalar(a)
                a = MAMaS.sym.matrix(a);
            end
            
            if (isa(b, 'MAMaS.sym') && b.isMatrix) || ~isscalar(b)
                error('Exponent must not be a matrix.');
            end
            if isa(a, 'MAMaS.sym') && a.isMatrix
                if a.matrixRows ~= a.matrixCols
                    error('Matrix exponentiation requires a square matrix.');
                end
                y = MAMaS.sym.binaryOp(a, b, '^^');
            else
                y = MAMaS.sym.binaryOp(a, b, '^');
            end
        end

        function y = uminus(a)
            if ~isscalar(a)
                a = MAMaS.sym.matrix(a);
            end

            y = MAMaS.sym.unaryOp(a, '-');
        end

        function y = uplus(a)
            y = a;
        end

        function y = transpose(obj)
            % Array transpose operator .'
            % For symbolic matrix, calls matrixTranspose
            if ~isscalar(obj)
                error('Operations with matrices of syms are not supported. Please consider using symbolic matrices instead.');
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
            if ~isscalar(obj)
                error('Operations with matrices of syms are not supported. Please consider using symbolic matrices instead.');
            end
            
            y = MAMaS.sym.funcOp('conjugate', obj);
            if obj.isMatrix
                % Transpose the Maxima matrix
                y = y.matrixTranspose();
            end
        end

        % Common math functions
        function y = sin(x)
            if ~isscalar(x)
                y = arrayfun(@sin, x);
                return
            end
            y = MAMaS.sym.funcOp('sin', x);
        end

        function y = cos(x)
            if ~isscalar(x)
                y = arrayfun(@cos, x);
                return
            end
            y = MAMaS.sym.funcOp('cos', x);
        end

        function y = tan(x)
            if ~isscalar(x)
                y = arrayfun(@tan, x);
                return
            end
            y = MAMaS.sym.funcOp('tan', x);
        end

        function y = asin(x)
            if ~isscalar(x)
                y = arrayfun(@asin, x);
                return
            end
            y = MAMaS.sym.funcOp('asin', x);
        end

        function y = acos(x)
            if ~isscalar(x)
                y = arrayfun(@acos, x);
                return
            end
            y = MAMaS.sym.funcOp('acos', x);
        end

        function y = atan(x)
            if ~isscalar(x)
                y = arrayfun(@atan, x);
                return
            end
            y = MAMaS.sym.funcOp('atan', x);
        end

        function y = exp(x)
            if ~isscalar(x)
                y = arrayfun(@exp, x);
                return
            end
            y = MAMaS.sym.funcOp('exp', x);
        end

        function y = expm(x)
            % Matrix exponential or scalar exponential
            % If x is a matrix, computes matrix exponential using matrixexp
            % Otherwise, computes scalar exponential
            if ~isscalar(x)
                x = MAMaS.sym.matrix(x);
            end

            if isa(x, 'MAMaS.sym') && x.isMatrix
                x.validateMaximaInstance();
                if x.matrixRows ~= x.matrixCols
                    error('Matrix exponential requires a square matrix.');
                end
                cmd = ['matrixexp(', x.identifier, ')'];
                id = x.maximaInstance.sendNoWait(cmd);
                y = MAMaS.sym.fromId(id, x.maximaInstance, x.getDimensions);
            else
                y = MAMaS.sym.funcOp('exp', x);
            end
        end

        function y = log(x)
            if ~isscalar(x)
                y = arrayfun(@log, x);
                return
            end        
            y = MAMaS.sym.funcOp('log', x);
        end

        function y = sqrt(x)
            if ~isscalar(x)
                y = arrayfun(@sqrt, x);
                return
            end
            y = MAMaS.sym.funcOp('sqrt', x);
        end

        function y = abs(x)
            if ~isscalar(x)
                y = arrayfun(@abs, x);
                return
            end
            y = MAMaS.sym.funcOp('abs', x);
        end

        function y = ratsimp(x)
            if ~isscalar(x)
                y = arrayfun(@ratsimp, x);
                return
            end
            y = MAMaS.sym.funcOp('ratsimp', x);
        end

        function y = trigsimp(x)
            if ~isscalar(x)
                y = arrayfun(@trigsimp, x);
                return
            end
            y = MAMaS.sym.funcOp('trigsimp', x);
        end

        function y = simplify(x)
            % Calls trigsimp(ratsimp(x)) in Maxima
            if ~isscalar(x)
                y = arrayfun(@simplify, x);
                return;
            end
            exprX = MAMaS.sym.toMaxima(x);
            
            % Get and validate instance from operand
            if isa(x, 'MAMaS.sym')
                x.validateMaximaInstance();
                maxima = x.maximaInstance;
            else
                maxima = MAMaS.MaximaInterface.getInstance();
            end
            
            cmd = ['trigsimp(ratsimp(', exprX, '))']; 
            id = maxima.sendNoWait(cmd);
            y = MAMaS.sym.fromId(id, maxima, x.getDimensions());
        end

        function y = diff(x, var, n)
            arguments
                x (1,1) MAMaS.sym
                var (1,1) MAMaS.sym
                n (1,1) double {mustBeInteger, mustBePositive} = 1
            end
            % Differentiate expression with respect to a variable
            % y = diff(expr, var)      - first derivative
            % y = diff(expr, var, n)   - nth derivative
            if ~isscalar(x)
                x = MAMaS.sym.matrix(x);
            end
            if ~isscalar(var)
                error('Can only differentiate with respect to a single variable.');
            end
            
            if ~var.isSymbol || var.isMatrix
                error('The variable argument must be a sym scalar symbol.');
            end
            if n==1
                y = MAMaS.sym.func('diff', x, var);
            else
                y = MAMaS.sym.func('diff', x, var, n);
            end
            if x.isMatrix
                y.isMatrix = true;
                y.matrixRows = x.matrixRows;
                y.matrixCols = x.matrixCols;
            end
        end

        function y = jacobian(exprs, vars)
            arguments
                exprs MAMaS.sym
                vars MAMaS.sym
            end
            % Compute Jacobian matrix (matrix of partial derivatives)
            % y = jacobian(exprs, vars)   - Jacobian of expressions w.r.t. variables
            % exprs: symbolic expressions
            % vars:  symbols
            
            if isscalar(exprs)
                if exprs.isMatrix
                    exprsStr = ['flatten(args(', exprs.identifier, '))'];
                    numRows = exprs.matrixRows*exprs.matrixCols;
                else
                    exprsStr = ['[', exprs.identifier, ']'];
                    numRows = 1;
                end
                maxima = exprs.maximaInstance;
            else
                if any(arrayfun(@(e) e.isMatrix, exprs))
                    error('All elements in exprs must be scalar syms.');
                end
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
            if isscalar(vars)
                if vars.maximaInstance ~= maxima
                    error('All sym arguments must belong to the same Maxima instance.');
                end
                if vars.isMatrix
                    error('vars cannot be a matrix symbol.');
                end
                if ~vars.isSymbol
                    error('vars must be a scalar symbol.');
                end

                varsStr = ['[', vars.identifier, ']'];
                numCols = 1;
            else
                for i = 1:numel(vars)
                    if vars(i).maximaInstance ~= maxima
                        error('All sym arguments must belong to the same Maxima instance.');
                    end
                end
                if any(arrayfun(@(v) v.isMatrix, vars))
                    error('All elements in vars must be scalar syms.');
                end
                if any(arrayfun(@(v) ~v.isSymbol, vars))
                    error('All elements in vars must be scalar symbols.');
                end
                
                varsStr = ['[', strjoin(arrayfun(@(v) v.identifier, vars(:), 'UniformOutput', false), ', '), ']'];
                numCols = numel(vars);
            end
            
            exprs(1).validateMaximaInstance();
            if numRows == 1 && numCols == 1
                cmd = ['diff(', exprs.identifier, ', ', vars.identifier, ')'];
            else
                cmd = ['jacobian(', exprsStr, ', ', varsStr, ')'];
            end
            id = maxima.sendNoWait(cmd);

            y = MAMaS.sym.fromId(id, maxima, [numRows, numCols]);
        end

        function y = subs(x, old, new, fullrat)
            arguments
                x MAMaS.sym
                old
                new
                fullrat (1,1) logical = false
            end
            % Substitute symbols/expressions in x with new values
            % MATLAB API: subs(expr, old, new)

            if isempty(old) || isempty(new)
                y = x;
                return;
            end

            if ~isscalar(x)
                y = arrayfun(@(xi) subs(xi, old, new), x);
                return;
            end

            x.validateMaximaInstance();
            maxima = x.maximaInstance;

            oldSubsStr= MAMaS.sym.toMaximaSubstArg(old, maxima);
            newSubsStr= MAMaS.sym.toMaximaSubstArg(new, maxima);

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

            y = MAMaS.sym.fromId(id, maxima, x.getDimensions());
        end

        function vars = symvar(x)
            arguments
                x MAMaS.sym
            end
            % Find symbolic variables in expression
            % MATLAB API: symvar(expr) returns symbolic variables in alphabetical order
            % This implementation returns a string array of variable names
            
            if ~isscalar(x)
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
            if ~isscalar(obj)
                error('Operations with matrices of syms are not supported. Please consider using symbolic matrices instead.');
            end
            
            obj.validateMaximaInstance();
            
            if obj.isMatrix
                % For matrices, use Maxima's copymatrix to create a true copy
                cmd = ['copymatrix(', obj.identifier, ')'];
                id = obj.maximaInstance.sendNoWait(cmd);
                newObj = MAMaS.sym.fromId(id, obj.maximaInstance, obj.getDimensions());
            else
                % For non-matrices, shallow copy is sufficient
                newObj = MAMaS.sym(obj.identifier, "expr", obj.maximaInstance);
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
                obj (1,1) MAMaS.sym
                row (1,1) double {mustBeInteger, mustBePositive}
                col (1,1) double {mustBeInteger, mustBePositive}
            end
            
            if ~isscalar(obj) || ~obj.isMatrix
                error('matrixElement can only be called on matrix sym.');
            end
            
            if row < 1 || row > obj.matrixRows || col < 1 || col > obj.matrixCols
                error('Matrix index (%d, %d) out of bounds [1:%d, 1:%d].', row, col, obj.matrixRows, obj.matrixCols);
            end
            
            obj.validateMaximaInstance();
            cmd = [obj.identifier, '[', num2str(row), ',', num2str(col), ']'];
            id = obj.maximaInstance.sendNoWait(cmd);
            y = MAMaS.sym.fromId(id, obj.maximaInstance, [1, 1]);
        end

        function y = matrixTranspose(obj)
            arguments
                obj (1,1) MAMaS.sym
            end 
            % Transpose a matrix
            
            if ~isscalar(obj) || ~obj.isMatrix
                error('matrixTranspose can only be called on matrix sym.');
            end
            
            obj.validateMaximaInstance();
            cmd = ['transpose(', obj.identifier, ')'];
            id = obj.maximaInstance.sendNoWait(cmd);
            
            y = MAMaS.sym.fromId(id, obj.maximaInstance, obj.getDimensions());
        end

        function y = matrixDeterminant(obj)
            arguments
                obj (1,1) MAMaS.sym
            end
            % Calculate matrix determinant
            
            if ~isscalar(obj) || ~obj.isMatrix
                error('matrixDeterminant can only be called on matrix sym.');
            end
            
            if obj.matrixRows ~= obj.matrixCols
                error('Determinant requires a square matrix.');
            end
            
            obj.validateMaximaInstance();
            cmd = ['determinant(', obj.identifier, ')'];
            id = obj.maximaInstance.sendNoWait(cmd);
            y = MAMaS.sym.fromId(id, obj.maximaInstance, [1, 1]);
        end

        function y = matrixInverse(obj)
            arguments
                obj (1,1) MAMaS.sym
            end
            
            if ~obj.isMatrix
                error('matrixInverse can only be called on matrix sym.');
            end
            
            if obj.matrixRows ~= obj.matrixCols
                error('Inverse requires a square matrix.');
            end
            
            obj.validateMaximaInstance();
            cmd = ['invert(', obj.identifier, ')'];
            id = obj.maximaInstance.sendNoWait(cmd);
            
            y = MAMaS.sym.fromId(id, obj.maximaInstance, obj.getDimensions());
        end

        function depends(var, varargin)
            arguments
                var MAMaS.sym
            end
            arguments (Repeating)
                varargin
            end
            % Declare dependencies between variables
            % var: sym that is a symbol (dependent variable)
            % varargin: list, array or cell array of symbols (independent variables)
            
            if isempty(varargin)
                error('At least one independent variable must be supplied.')
            end

            if ~isscalar(var)
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
                if isscalar(args)
                    args = {args};
                else
                    args = num2cell(args);
                end
            end
            
            % Validate all independent variables are symbols
            for i = 1:length(args)
                if ~isa(args{i}, 'MAMaS.sym') || ~args{i}.isSymbol
                    error('All elements in indepVars must be symbols.');
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

        function varargout = subsref(obj, s)
            % Subscript reference - handle indexing like obj(i,j)
            if isempty(s)
                varargout = {obj};
                return;
            end
            if ~isscalar(obj)
                varargout = {builtin('subsref', obj, s)};
                return
            end
            
            switch s(1).type
                case '()'
                    % special case, mostly for indexing element 1 of a
                    % scalar object
                    if isscalar(s(1).subs)
                        varargout = {builtin('subsref', obj, s)};
                        return
                    end

                    if length(s(1).subs) ~= 2
                        error('Matrix indexing requires exactly 2 subscripts (row, col).');
                    end

                    % Handle parentheses indexing
                    rows = s(1).subs{1};
                    cols = s(1).subs{2};
                    if ischar(rows) && rows==':'
                        rows = 1:obj.matrixRows;
                    end
                    if ischar(cols) && cols==':'
                        cols = 1:obj.matrixCols;
                    end

                    if ~obj.isMatrix
                        if isscalar(rows) && rows==1 && isscalar(cols) && cols==1
                            varargout = {obj};
                            return
                        else
                            error('Subscript indexing only supported for matrices.');
                        end
                    end
                    
                    if ~isscalar(rows) || ~isscalar(cols)
                        % compose a command string to create a new matrix from the selected elements
                        dataStr = cell(numel(rows), numel(cols));
                        for i = 1:numel(rows)
                            for j = 1:numel(cols)
                                dataStr{i,j} = [obj.identifier, '[', num2str(rows(i)), ',', num2str(cols(j)), ']'];
                            end
                        end

                        % Build Maxima matrix syntax: matrix([row1], [row2], ...)
                        rowStrs = cell(numel(rows), 1);
                        for i = 1:numel(rows)
                            rowStr = strjoin(dataStr(i,:), ', ');
                            rowStrs{i} = ['[', rowStr, ']'];
                        end
                        matrixStr = strjoin(rowStrs, ', ');
                        
                        % Send matrix creation command
                        cmd = ['matrix(', matrixStr, ')'];
                        obj.validateMaximaInstance();
                        id = obj.maximaInstance.sendNoWait(cmd);
                        
                        % Create sym object and set matrix properties
                        y = MAMaS.sym.fromId(id, obj.maximaInstance, [numel(rows), numel(cols)]);
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
                    
                case '.'
                    % Handle dot notation (property access)
                    varargout = {builtin('subsref', obj, s)};
                    
                case '{}'
                    % Cell array indexing not supported
                    error('Cell array indexing {} not supported for syms.');
                    
                otherwise
                    varargout = {builtin('subsref', obj, s)};
            end
        end

        function obj = subsasgn(obj, s, val)
            % Subscript assignment - handle assignment like obj(i,j) = expr
            if isempty(s)
                return;
            end
            if ~isscalar(obj)
                obj = builtin('subsasgn', obj, s, val);
                return
            end
            
            switch s(1).type
                case '()'
                    % Handle parentheses indexing assignment
                    if ~obj.isMatrix
                        error('Subscript assignment only supported for matrices.');
                    end
                    
                    if length(s(1).subs) ~= 2
                        error('Matrix assignment requires exactly 2 subscripts (row, col).');
                    end
                    
                    rows = s(1).subs{1};
                    cols = s(1).subs{2};
                    if ischar(rows) && rows==':'
                        rows = 1:obj.matrixRows;
                    end
                    if ischar(cols) && cols==':'
                        cols = 1:obj.matrixCols;
                    end
                    
                    % Validate indices
                    if any(rows < 1) || any(rows > obj.matrixRows) || any(cols < 1) || any(cols > obj.matrixCols)
                        error('One ore more matrix indeces out of bounds [1:%d, 1:%d]. Automatic expansion currently not supported.', obj.matrixRows, obj.matrixCols);
                    end
                    
                    dataStr = cell(numel(rows), numel(cols));
                    if isa(val, 'MAMaS.sym') && val.isMatrix
                        if numel(rows) ~= val.matrixRows || numel(cols) ~= val.matrixCols
                            error('Assigned matrix size does not match target submatrix size.');
                        end
                        for i = 1:numel(rows)
                            for j = 1:numel(cols)
                                dataStr{i,j} = [val.identifier, '[', num2str(rows(i)), ',', num2str(cols(j)), ']'];
                            end
                        end
                    else
                        if ischar(val)
                            val = cellstr(val);
                        end
                        if ~iscell(val)
                            if isscalar(val)
                                val = {val};
                            else
                                val = num2cell(val);
                            end
                        end
                        if size(val, 1) ~= numel(rows) || size(val, 2) ~= numel(cols)
                            error('Assigned value size does not match target submatrix size.');
                        end
                        for i = 1:numel(rows)
                            for j = 1:numel(cols)
                                dataStr{i,j} = MAMaS.sym.toMaxima(val{i,j});
                            end
                        end
                    end

                    % Create assignment command in Maxima
                    obj.validateMaximaInstance();
                    for i = 1:numel(rows)
                        for j = 1:numel(cols)
                            cmd = [obj.identifier, '[', num2str(rows(i)), ',', num2str(cols(j)), '] : ', dataStr{i,j}];
                            obj.maximaInstance.sendNoWait(cmd);
                        end
                    end
            
                    % Clear cached output since matrix changed
                    obj.cachedOutput = '';
                    
                case '.'
                    % Handle dot notation (property assignment)
                    obj = builtin('subsasgn', obj, s, val);
                    
                case '{}'
                    % Cell array indexing not supported
                    error('Cell array indexing {} not supported for syms.');
                    
                otherwise
                    obj = builtin('subsasgn', obj, s, val);
            end
        end

        function idx = subsindex(obj)
            % Convert sym to an index for use in subscripting
            % For numeric expressions, convert to a MATLAB number
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
            
            if isscalar(obj)
                % Single sym object
                if ~obj.isMatrix
                    error('end() indexing only supported for matrices.');
                end
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
                        if isa(arg, 'MAMaS.sym')
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
                argStrs = cellfun(@MAMaS.sym.toMaxima, arg, 'UniformOutput', false);
                exprArgs = strjoin(argStrs, ', ');

                for i = 1:length(arg)
                    if isa(arg{i}, 'MAMaS.sym')
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
                maxima = MAMaS.MaximaInterface.getInstance();
            end
            
            cmd = [fname, '(', exprArgs, ')'];
            id = maxima.sendNoWait(cmd);
            res = maxima.sendAndParse(sprintf('if matrixp(%s) then [length(%s), length(first(%s))] else [0, 0]', id, id, id));
            % TODO: handle 1x1 matrices properly
            % parse dimensions from result string "[r, c]"
            dims = str2num(res{1}); %#ok<ST2NM>
            y = MAMaS.sym.fromId(id, maxima, dims);
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

            obj = MAMaS.sym(id, "const");
        end

        function obj = pi()
            obj = MAMaS.sym.const('pi');
        end

        function obj = e()
            obj = MAMaS.sym.const('e');
        end

        function obj = gamma()
            obj = MAMaS.sym.const('gamma');
        end

        function obj = i()
            obj = MAMaS.sym.const('i');
        end

        function obj = inf()
            obj = MAMaS.sym.const('inf');
        end

        function obj = minf()
            obj = MAMaS.sym.const('minf');
        end

        function obj = phi()
            obj = MAMaS.sym.const('phi');
        end

        function y = matrix(data, maxima)
            % Create a Maxima matrix from MATLAB data
            % data: MATLAB matrix of numbers, strings, or syms
            % maxima: optional MaximaInterface instance
            
            arguments
                data
                maxima = []
            end
            
            if isscalar(data)
                error('Will not produce a 1x1 matrix. Use sym() constructor instead.');
            end

            if ~isnumeric(data) && ~isstring(data) && ~iscell(data) && ~isa(data, 'MAMaS.sym')
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
                    if isa(data{i,j}, 'MAMaS.sym')
                        if isempty(maxima)
                            data{i,j}.validateMaximaInstance();
                            maxima = data{i,j}.maximaInstance;
                        else
                            if maxima ~= data{i,j}.maximaInstance
                                error('All sym elements must belong to the same MaximaInterface instance.');
                            end
                        end
                    end
                    dataStr{i,j} = MAMaS.sym.toMaxima(data{i,j});
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
            y = MAMaS.sym.fromId(id, maxima, size(data));
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
            y = MAMaS.sym.fromId(id, maxima, [rows, cols]);
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
            
            y = MAMaS.sym.fromId(id, maxima, [rows, cols]);
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
            
            y = MAMaS.sym.fromId(id, maxima, [rows, cols]);
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
            
            y = MAMaS.sym.fromId(id, maxima, [n, n]);
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
            if obj.isMatrix
                dims = [obj.matrixRows, obj.matrixCols];
            else
                dims = [1, 1];
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

            if isa(val, 'MAMaS.sym')
                if isscalar(val)
                    if val.maximaInstance ~= maxima
                        error('All sym arguments must belong to the same Maxima instance.');
                    end
                    if val.isMatrix
                        error('Scalar sym expected for substitution argument, but matrix found.');
                    else
                        subsStr = {val.identifier};
                    end
                    return;
                end

                if any(arrayfun(@(v) v.isMatrix, val))
                    error('Array inputs to subs cannot contain matrix syms.');
                end
                for i = 1:numel(val)
                    if val(i).maximaInstance ~= maxima
                        error('All sym arguments must belong to the same Maxima instance.');
                    end
                end

                subsStr = arrayfun(@(v) v.identifier, val(:), 'UniformOutput', false);
                return;
            end

            if iscell(val)
                subsStr = cellfun(@(v) MAMaS.sym.toMaxima(v), val(:), 'UniformOutput', false);
                return;
            end

            if ischar(val) || isscalar(val)
                subsStr = {MAMaS.sym.toMaxima(val)};
                return;
            end

            if isnumeric(val) || isstring(val)
                subsStr = arrayfun(@(v) MAMaS.sym.toMaxima(v), val(:), 'UniformOutput', false);
                return
            end

            error('Unsupported substitution argument type.');
        end

        function y = binaryOp(a, b, op)
            exprA = MAMaS.sym.toMaxima(a);
            exprB = MAMaS.sym.toMaxima(b);
            
            % Get and validate instance from first sym operand
            maxima = MAMaS.sym.getMaximaFromOperands(a, b);
            
            % no need for parentheses because operands are always symbols or maxima identifiers
            cmd = [exprA, op, exprB];
            id = maxima.sendNoWait(cmd);


            % Handle element-wise matrix operations
            dims = [1, 1];
            if isa(a, 'MAMaS.sym') && a.isMatrix
                dims = a.getDimensions();
            elseif isa(b, 'MAMaS.sym') && b.isMatrix
                dims = b.getDimensions();
            end

            y = MAMaS.sym.fromId(id, maxima, dims);
        end

        function y = unaryOp(a, op)
            exprA = MAMaS.sym.toMaxima(a);
            
            % Get and validate instance from operand
            if isa(a, 'MAMaS.sym')
                a.validateMaximaInstance();
                maxima = a.maximaInstance;
            else
                maxima = MaximaInterface.getInstance();
            end
            
            % no need for parentheses because operands are always symbols or maxima identifiers
            cmd = [op, exprA];
            id = maxima.sendNoWait(cmd);

            dims = [1, 1];
            if isa(a, 'MAMaS.sym') && a.isMatrix
                dims = a.getDimensions();
            end
            y = MAMaS.sym.fromId(id, maxima, dims);
        end
        
        function y = funcOp(fname, arg)
            % Apply a function operation to a single argument
            % fname: function name (string)
            % arg: argument (can be sym, number, or string)
            
            exprArg = MAMaS.sym.toMaxima(arg);
            
            % Get and validate instance from operand
            if isa(arg, 'MAMaS.sym')
                arg.validateMaximaInstance();
                maxima = arg.maximaInstance;
            else
                maxima = MAMaS.MaximaInterface.getInstance();
            end
            
            cmd = [fname, '(', exprArg, ')'];
            id = maxima.sendNoWait(cmd);
            y = MAMaS.sym.fromId(id, maxima, arg.getDimensions());
        end

        function maxima = getMaximaFromOperands(a, b)
            % Get MaximaInterface instance from operands and validate
            if isa(a, 'MAMaS.sym') && isa(b, 'MAMaS.sym')
                if a.maximaInstance ~= b.maximaInstance
                    error('Both sym operands must belong to the same MaximaInterface instance.');
                end
                % if both instances are the same, vaidating one is sufficient
                a.validateMaximaInstance();
                maxima = a.maximaInstance;
            elseif isa(a, 'MAMaS.sym')
                a.validateMaximaInstance();
                maxima = a.maximaInstance;
            elseif isa(b, 'MAMaS.sym')
                b.validateMaximaInstance();
                maxima = b.maximaInstance;
            else
                maxima = MAMaS.MaximaInterface.getInstance();
            end
        end

        function s = toMaxima(x)
            if isa(x, 'MAMaS.sym')
                s = x.identifier;
            elseif isnumeric(x) && isscalar(x)
                s = num2str(x, 16);
            elseif ischar(x) || (isstring(x) && isscalar(x))
                x = char(x);
                if ~MAMaS.sym.isValidSymbolName(x)
                    error('String operands must be valid symbol names.');
                end
                s = x;
            else
                error('Unsupported operand type.');
            end
        end

        function obj = fromId(id, maxima, dims)
            obj = MAMaS.sym(id, dims, "expr", maxima);
        end
    end
end
