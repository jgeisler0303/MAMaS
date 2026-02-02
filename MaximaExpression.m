classdef MaximaExpression < handle
    % MaximaExpression
    % MATLAB class for symbolic Maxima expressions

    properties (SetAccess = private)
        identifier
        isNumber (1,1) logical = false       % True if this is a numeric value
        isSymbol (1,1) logical = false       % True if this is a symbol
        isExpression (1,1) logical = false   % True if this is a compound expression
        isMatrix (1,1) logical = false       % True if this is a Maxima matrix
        matrixRows = 0                       % Number of rows (0 if not a matrix)
        matrixCols = 0                       % Number of columns (0 if not a matrix)
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
                % For internal expressions, assume they are expressions (result of operations)
                obj.isExpression = true;
                return;
            end

            % Check if input is matrix data (numeric array, string array, or cell array)
            if ~ischar(name) && ~isscalar(name)
                % Delegate to matrix creation
                if isempty(maxima)
                    maxima = MaximaInterface.getInstance();
                end
                matObj = MaximaExpression.matrix(name, maxima);
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
            numValue = str2double(name);
            if ~isnan(numValue)
                name = numValue;
            end
            if isnumeric(name) && isscalar(name)
                % Convert numeric input to Maxima representation
                obj.identifier = num2str(name, 16);
                obj.cachedOutput = '';
                obj.maximaInstance = MaximaInterface.getInstance();
                obj.isNumber = true;
                return;
            end

            if ~(ischar(name) || (isstring(name) && isscalar(name)))
                error('Input must be a string or numeric scalar.');
            end
            name = char(name);

            if ~MaximaExpression.isValidSymbolName(name)
                error('Invalid symbol name: "%s".', name);
            end

            obj.identifier = name;
            obj.cachedOutput = '';
            obj.maximaInstance = MaximaInterface.getInstance();
            obj.isSymbol = true;
        end

        function s = char(obj)
            if numel(obj) == 1
                s = obj.getDisplayString();
            else
                % For arrays, create cell array of strings
                s = arrayfun(@(x) x.getDisplayString(), obj, 'UniformOutput', false);
            end
        end

        function s = string(obj)
            if numel(obj) == 1
                s = string(obj.getDisplayString());
            else
                % For arrays, create string array
                s = arrayfun(@(x) string(x.getDisplayString()), obj);
            end
        end

        function disp(obj)
            if numel(obj) == 1
                disp(obj.getDisplayString());
            else
                % For arrays, display size and elements
                fprintf('%dx%d MaximaExpression array:\n\n', size(obj, 1), size(obj, 2));
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
            if ~isscalar(a) || ~isscalar(b)
                error('Operations with matrices of MaximaExpressions are not supported. Please consider using symbolic matrices instead.');
            end
            if isa(a, 'MaximaExpression') && a.isMatrix && isa(b, 'MaximaExpression') && b.isMatrix
                if a.matrixRows ~= b.matrixRows || a.matrixCols ~= b.matrixCols
                    error('Matrix dimensions mismatch for addition: [%d x %d] + [%d x %d].', a.matrixRows, a.matrixCols, b.matrixRows, b.matrixCols);
                end
            end
            y = MaximaExpression.binaryOp(a, b, '+');
        end

        function y = minus(a, b)
            if ~isscalar(a) || ~isscalar(b)
                error('Operations with matrices of MaximaExpressions are not supported. Please consider using symbolic matrices instead.');
            end
            if isa(a, 'MaximaExpression') && a.isMatrix && isa(b, 'MaximaExpression') && b.isMatrix
                if a.matrixRows ~= b.matrixRows || a.matrixCols ~= b.matrixCols
                    error('Matrix dimensions mismatch for subtraction: [%d x %d] - [%d x %d].', a.matrixRows, a.matrixCols, b.matrixRows, b.matrixCols);
                end
            end
            y = MaximaExpression.binaryOp(a, b, '-');
        end

        function y = times(a, b)
            % Element-wise multiplication .*
            if ~isscalar(a) || ~isscalar(b)
                error('Operations with matrices of MaximaExpressions are not supported. Please consider using symbolic matrices instead.');
            end

            if isa(a, 'MaximaExpression') && a.isMatrix && isa(b, 'MaximaExpression') && b.isMatrix
                if a.matrixRows ~= b.matrixRows || a.matrixCols ~= b.matrixCols
                    error('Matrix dimensions mismatch for element-wise multiplication: [%d x %d] .* [%d x %d].', a.matrixRows, a.matrixCols, b.matrixRows, b.matrixCols);
                end
            end
            y = MaximaExpression.binaryOp(a, b, '*');
        end

        function y = mtimes(a, b)
            % Matrix multiplication uses dot operator in Maxima (a . b)
            % Check if both operands are matrices
            if ~isscalar(a) || ~isscalar(b)
                error('Operations with matrices of MaximaExpressions are not supported. Please consider using symbolic matrices instead.');
            end

            aIsMatrix = isa(a, 'MaximaExpression') && a.isMatrix;
            bIsMatrix = isa(b, 'MaximaExpression') && b.isMatrix;
            
            if aIsMatrix && bIsMatrix
                % Matrix-matrix multiplication uses "." in Maxima
                if a.matrixCols ~= b.matrixRows
                    error('Matrix dimensions mismatch for multiplication: [%d x %d] * [%d x %d].', a.matrixRows, a.matrixCols, b.matrixRows, b.matrixCols);
                end
                y = MaximaExpression.binaryOp(a, b, '.');
                % Result is a matrix with dimensions from outer dimensions
                y.isMatrix = true;
                y.matrixRows = a.matrixRows;
                y.matrixCols = b.matrixCols;
            else
                % Scalar multiplication uses "*"
                y = MaximaExpression.binaryOp(a, b, '*');
            end
        end

        function y = rdivide(a, b)
            % Element-wise right division ./
            if ~isscalar(a) || ~isscalar(b)
                error('Operations with matrices of MaximaExpressions are not supported. Please consider using symbolic matrices instead.');
            end

            if isa(a, 'MaximaExpression') && a.isMatrix && isa(b, 'MaximaExpression') && b.isMatrix
                if a.matrixRows ~= b.matrixRows || a.matrixCols ~= b.matrixCols
                    error('Matrix dimensions mismatch for element-wise multiplication: [%d x %d] .* [%d x %d].', a.matrixRows, a.matrixCols, b.matrixRows, b.matrixCols);
                end
            end
            y = MaximaExpression.binaryOp(a, b, '/');
        end

        function y = mrdivide(a, b)
            % Matrix right division a / b
            if ~isscalar(a) || ~isscalar(b)
                error('Operations with matrices of MaximaExpressions are not supported. Please consider using symbolic matrices instead.');
            end

            if isa(b, 'MaximaExpression') && b.isMatrix
                error('Matrix right division is not supported.');
            end
            y = MaximaExpression.binaryOp(a, b, '/');
        end

        function y = power(a, b)
            % Element-wise power .^
            if ~isscalar(a) || ~isscalar(b)
                error('Operations with matrices of MaximaExpressions are not supported. Please consider using symbolic matrices instead.');
            end

            if isa(a, 'MaximaExpression') && a.isMatrix && isa(b, 'MaximaExpression') && b.isMatrix
                error('Base and exponent cannot both be a matrix.');
            end
            y = MaximaExpression.binaryOp(a, b, '^');
        end

        function y = mpower(a, b)
            % Matrix power a ^ b
            if ~isscalar(a) || ~isscalar(b)
                error('Operations with matrices of MaximaExpressions are not supported. Please consider using symbolic matrices instead.');
            end
            
            if isa(b, 'MaximaExpression') && b.isMatrix
                error('Exponent must not be a matrix.');
            end
            if isa(a, 'MaximaExpression') && a.isMatrix
                if isa(b, 'MaximaExpression') && ~b.isMatrix
                    error('Matrix exponent must be scalar.');
                end
                if a.matrixRows ~= a.matrixCols
                    error('Matrix exponentiation requires a square matrix.');
                end
                y = MaximaExpression.binaryOp(a, b, '^^');
            else
                y = MaximaExpression.binaryOp(a, b, '^');
            end
        end

        function y = uminus(a)
            if ~isscalar(a)
                error('Operations with matrices of MaximaExpressions are not supported. Please consider using symbolic matrices instead.');
            end

            y = MaximaExpression.unaryOp(a, '-');
        end

        function y = uplus(a)
            y = a;
        end

        function y = transpose(obj)
            % Array transpose operator .'
            % For MaximaExpression matrix, calls matrixTranspose
            if ~isscalar(obj)
                error('Operations with matrices of MaximaExpressions are not supported. Please consider using symbolic matrices instead.');
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
            % For MaximaExpression (real symbolic expressions), same as transpose
            % For MaximaExpression matrix, calls matrixTranspose
            if ~isscalar(obj)
                error('Operations with matrices of MaximaExpressions are not supported. Please consider using symbolic matrices instead.');
            end
            
            y = MaximaExpression.funcOp('conjugate', obj);
            if obj.isMatrix
                % Transpose the Maxima matrix
                y = obj.matrixTranspose();
            else
                y = obj;
            end
        end

        % Common math functions
        function y = sin(x)
            if ~isscalar(x)
                error('Operations with matrices of MaximaExpressions are not supported. Please consider using symbolic matrices instead.');
            end
            y = MaximaExpression.funcOp('sin', x);
        end

        function y = cos(x)
            if ~isscalar(x)
                error('Operations with matrices of MaximaExpressions are not supported. Please consider using symbolic matrices instead.');
            end
            y = MaximaExpression.funcOp('cos', x);
        end

        function y = tan(x)
            if ~isscalar(x)
                error('Operations with matrices of MaximaExpressions are not supported. Please consider using symbolic matrices instead.');
            end
            y = MaximaExpression.funcOp('tan', x);
        end

        function y = asin(x)
            if ~isscalar(x)
                error('Operations with matrices of MaximaExpressions are not supported. Please consider using symbolic matrices instead.');
            end
            y = MaximaExpression.funcOp('asin', x);
        end

        function y = acos(x)
            if ~isscalar(x)
                error('Operations with matrices of MaximaExpressions are not supported. Please consider using symbolic matrices instead.');
            end
            y = MaximaExpression.funcOp('acos', x);
        end

        function y = atan(x)
            if ~isscalar(x)
                error('Operations with matrices of MaximaExpressions are not supported. Please consider using symbolic matrices instead.');
            end
            y = MaximaExpression.funcOp('atan', x);
        end

        function y = exp(x)
            if ~isscalar(x)
                error('Operations with matrices of MaximaExpressions are not supported. Please consider using symbolic matrices instead.');
            end
            y = MaximaExpression.funcOp('exp', x);
        end

        function y = expm(x)
            % Matrix exponential or scalar exponential
            % If x is a matrix, computes matrix exponential using matrixexp
            % Otherwise, computes scalar exponential
            if ~isscalar(x)
                error('Operations with matrices of MaximaExpressions are not supported. Please consider using symbolic matrices instead.');
            end

            if isa(x, 'MaximaExpression') && x.isMatrix
                x.validateMaximaInstance();
                if x.matrixRows ~= x.matrixCols
                    error('Matrix exponential requires a square matrix.');
                end
                cmd = ['matrixexp(', x.identifier, ')'];
                id = x.maximaInstance.sendNoWait(cmd);
                y = MaximaExpression.fromId(id, x.maximaInstance);
                y.isMatrix = true;
                y.matrixRows = x.matrixRows;
                y.matrixCols = x.matrixCols;
            else
                y = MaximaExpression.funcOp('exp', x);
            end
        end

        function y = log(x)
            if ~isscalar(x)
                error('Operations with matrices of MaximaExpressions are not supported. Please consider using symbolic matrices instead.');
            end        
            y = MaximaExpression.funcOp('log', x);
        end

        function y = sqrt(x)
            if ~isscalar(x)
                error('Operations with matrices of MaximaExpressions are not supported. Please consider using symbolic matrices instead.');
            end
            y = MaximaExpression.funcOp('sqrt', x);
        end

        function y = abs(x)
            if ~isscalar(x)
                error('Operations with matrices of MaximaExpressions are not supported. Please consider using symbolic matrices instead.');
            end
            y = MaximaExpression.funcOp('abs', x);
        end

        function y = ratsimp(x)
            if ~isscalar(x)
                error('Operations with matrices of MaximaExpressions are not supported. Please consider using symbolic matrices instead.');
            end
            y = MaximaExpression.funcOp('ratsimp', x);
        end

        function y = trigsimp(x)
            if ~isscalar(x)
                error('Operations with matrices of MaximaExpressions are not supported. Please consider using symbolic matrices instead.');
            end
            y = MaximaExpression.funcOp('trigsimp', x);
        end

        function y = simplify(x)
            % Calls trigsimp(ratsimp(x)) in Maxima
            if ~isscalar(x)
                error('Operations with matrices of MaximaExpressions are not supported. Please consider using symbolic matrices instead.');
            end
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
            % This should never be called obj being an array
            if isempty(obj.maximaInstance) || ~isvalid(obj.maximaInstance)
                error('MaximaExpression: The Maxima interface instance has been deleted. Expression is invalid.');
            end
        end

        function newObj = copy(obj)
            % Create a copy of the MaximaExpression
            % For matrices, creates a new matrix in Maxima using copymatrix
            % For other expressions, creates a shallow copy with the same identifier
            if ~isscalar(obj)
                error('Operations with matrices of MaximaExpressions are not supported. Please consider using symbolic matrices instead.');
            end
            
            obj.validateMaximaInstance();
            
            if obj.isMatrix
                % For matrices, use Maxima's copymatrix to create a true copy
                cmd = ['copymatrix(', obj.identifier, ')'];
                id = obj.maximaInstance.sendNoWait(cmd);
                newObj = MaximaExpression.fromId(id, obj.maximaInstance);
                % Copy matrix properties
                newObj.isMatrix = true;
                newObj.matrixRows = obj.matrixRows;
                newObj.matrixCols = obj.matrixCols;
            else
                % For non-matrices, shallow copy is sufficient
                newObj = MaximaExpression(obj.identifier, true, obj.maximaInstance);
                % Copy type flags
                newObj.isNumber = obj.isNumber;
                newObj.isSymbol = obj.isSymbol;
                newObj.isExpression = obj.isExpression;
            end
            
            % Note: cachedOutput is intentionally not copied to allow fresh evaluation
        end

        function y = matrixElement(obj, row, col)
            % Get a matrix element
            % row: row index (1-based)
            % col: column index (1-based)
            arguments
                obj (1,1) MaximaExpression
                row (1,1) double {mustBeInteger, mustBePositive}
                col (1,1) double {mustBeInteger, mustBePositive}
            end
            
            if ~obj.isMatrix
                error('matrixElement can only be called on matrix MaximaExpressions.');
            end
            
            if row < 1 || row > obj.matrixRows || col < 1 || col > obj.matrixCols
                error('Matrix index (%d, %d) out of bounds [1:%d, 1:%d].', row, col, obj.matrixRows, obj.matrixCols);
            end
            
            obj.validateMaximaInstance();
            cmd = [obj.identifier, '[', num2str(row), ',', num2str(col), ']'];
            id = obj.maximaInstance.sendNoWait(cmd);
            y = MaximaExpression.fromId(id, obj.maximaInstance);
        end

        function y = matrixTranspose(obj)
            arguments
                obj (1,1) MaximaExpression
            end 
            % Transpose a matrix
            
            if ~obj.isMatrix
                error('matrixTranspose can only be called on matrix MaximaExpressions.');
            end
            
            obj.validateMaximaInstance();
            cmd = ['transpose(', obj.identifier, ')'];
            id = obj.maximaInstance.sendNoWait(cmd);
            
            y = MaximaExpression.fromId(id, obj.maximaInstance);
            y.isMatrix = true;
            y.matrixRows = obj.matrixCols;
            y.matrixCols = obj.matrixRows;
        end

        function y = matrixDeterminant(obj)
            arguments
                obj (1,1) MaximaExpression
            end
            % Calculate matrix determinant
            
            if ~obj.isMatrix
                error('matrixDeterminant can only be called on matrix MaximaExpressions.');
            end
            
            if obj.matrixRows ~= obj.matrixCols
                error('Determinant requires a square matrix.');
            end
            
            obj.validateMaximaInstance();
            cmd = ['determinant(', obj.identifier, ')'];
            id = obj.maximaInstance.sendNoWait(cmd);
            y = MaximaExpression.fromId(id, obj.maximaInstance);
        end

        function y = matrixInverse(obj)
            arguments
                obj (1,1) MaximaExpression
            end
            
            if ~obj.isMatrix
                error('matrixInverse can only be called on matrix MaximaExpressions.');
            end
            
            if obj.matrixRows ~= obj.matrixCols
                error('Inverse requires a square matrix.');
            end
            
            obj.validateMaximaInstance();
            cmd = ['invert(', obj.identifier, ')'];
            id = obj.maximaInstance.sendNoWait(cmd);
            
            y = MaximaExpression.fromId(id, obj.maximaInstance);
            y.isMatrix = true;
            y.matrixRows = obj.matrixRows;
            y.matrixCols = obj.matrixCols;
        end

        function depends(var, varargin)
            arguments
                var (1,1) MaximaExpression
            end
            arguments (Repeating)
                varargin
            end
            % Declare dependencies between variables
            % var: MaximaExpression that is a symbol (dependent variable)
            % varargin: list, array or cell array of MaximaExpression symbols (independent variables)
            
            % Validate var is a symbol
            if ~var.isSymbol
                error('First argument to depends must be a symbol.');
            end
            
            if isempty(varargin)
                error('At least one independent variable must be supplied.')
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
                if ~isa(args{i}, 'MaximaExpression') || ~args{i}.isSymbol
                    error('All elements in indepVars must be MaximaExpression symbols.');
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
                    
                    if length(s(1).subs) ~= 2
                        error('Matrix indexing requires exactly 2 subscripts (row, col).');
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
                        
                        % Create MaximaExpression object and set matrix properties
                        y = MaximaExpression.fromId(id, obj.maximaInstance);
                        y.isMatrix = true;
                        y.matrixRows = numel(rows);
                        y.matrixCols = numel(cols);
                    else
                        % Get one element
                        y = obj.matrixElement(rows, cols);
                    end

                    % Chained indexing doesn't make sense for symbols: throw error
                    if length(s) > 1
                        error('Chained indexing not supported for MaximaExpression symbols.');
                    else
                        varargout = {y};
                    end
                    
                case '.'
                    % Handle dot notation (property access)
                    varargout = {builtin('subsref', obj, s)};
                    
                case '{}'
                    % Cell array indexing not supported
                    error('Cell array indexing {} not supported for MaximaExpression.');
                    
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
                    if isa(val, 'MaximaExpression') && val.isMatrix
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
                                dataStr{i,j} = MaximaExpression.toMaxima(val{i,j});
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
                    error('Cell array indexing {} not supported for MaximaExpression.');
                    
                otherwise
                    obj = builtin('subsasgn', obj, s, val);
            end
        end

        function idx = subsindex(obj)
            % Convert MaximaExpression to an index for use in subscripting
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
                error('Only numeric MaximaExpressions can be used as array indices.');
            end
        end

        function idx = end(obj, k, n)
            % Support the use of 'end' keyword in subscripting
            % k: dimension being indexed (1 for rows, 2 for columns)
            % n: total number of subscripts
            
            if isscalar(obj)
                % Single MaximaExpression object
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
                % Array of MaximaExpression objects
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

        function y = matrix(data, maxima)
            % Create a Maxima matrix from MATLAB data
            % data: MATLAB matrix of numbers, strings, or MaximaExpressions
            % maxima: optional MaximaInterface instance
            
            arguments
                data
                maxima = []
            end
            
            if isempty(maxima)
                maxima = MaximaInterface.getInstance();
            end
            
            if ~isnumeric(data) && ~isstring(data) && ~iscell(data) && ~isa(data, 'MaximaExpression')
                error('Matrix data must be MaximaExpression, numeric, string array, cell array, or contain MaximaExpressions.');
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
                    dataStr{i,j} = MaximaExpression.toMaxima(data{i,j});
                end
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
            
            % Create MaximaExpression object and set matrix properties
            y = MaximaExpression.fromId(id, maxima);
            y.isMatrix = true;
            y.matrixRows = rows;
            y.matrixCols = cols;
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
            
            if isempty(maxima)
                maxima = MaximaInterface.getInstance();
            end
            
            name = char(name);
            
            % Create symbolic matrix using genmatrix
            cmd = ['genmatrix(lambda([i,j], ', name, '[i,j]), ', num2str(rows), ', ', num2str(cols), ')'];
            id = maxima.sendNoWait(cmd);
            
            % Create MaximaExpression object and set matrix properties
            y = MaximaExpression.fromId(id, maxima);
            y.isMatrix = true;
            y.matrixRows = rows;
            y.matrixCols = cols;
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
            
            if isempty(maxima)
                maxima = MaximaInterface.getInstance();
            end
            
            cmd = ['zeromatrix(', num2str(rows), ', ', num2str(cols), ')'];
            id = maxima.sendNoWait(cmd);
            
            y = MaximaExpression.fromId(id, maxima);
            y.isMatrix = true;
            y.matrixRows = rows;
            y.matrixCols = cols;
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
            
            if isempty(maxima)
                maxima = MaximaInterface.getInstance();
            end
            
            cmd = ['genmatrix(1, ', num2str(rows), ', ', num2str(cols), ')'];
            id = maxima.sendNoWait(cmd);
            
            y = MaximaExpression.fromId(id, maxima);
            y.isMatrix = true;
            y.matrixRows = rows;
            y.matrixCols = cols;
        end

        function y = eye(n, maxima)
            % Create an identity matrix
            % n: matrix dimension (n x n)
            % maxima: optional MaximaInterface instance
            
            arguments
                n (1,1) double {mustBeInteger, mustBePositive}
                maxima = []
            end
            
            if isempty(maxima)
                maxima = MaximaInterface.getInstance();
            end
            
            cmd = ['ident(', num2str(n), ')'];
            id = maxima.sendNoWait(cmd);
            
            y = MaximaExpression.fromId(id, maxima);
            y.isMatrix = true;
            y.matrixRows = n;
            y.matrixCols = n;
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

            % Handle element-wise matrix operations
            if isa(a, 'MaximaExpression') && a.isMatrix
                y.isMatrix = true;
                y.matrixRows = a.matrixRows;
                y.matrixCols = a.matrixCols;
            elseif isa(b, 'MaximaExpression') && b.isMatrix
                y.isMatrix = true;
                y.matrixRows = b.matrixRows;
                y.matrixCols = b.matrixCols;
            end
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

            if isa(a, 'MaximaExpression') && a.isMatrix
                y.isMatrix = true;
                y.matrixRows = a.matrixRows;
                y.matrixCols = a.matrixCols;
            end
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

            if isa(arg, 'MaximaExpression') && arg.isMatrix
                y.isMatrix = true;
                y.matrixRows = arg.matrixRows;
                y.matrixCols = arg.matrixCols;
            end
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
