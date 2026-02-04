% Extended test for sym
addpath(fullfile(fileparts(mfilename("fullpath")), '..', 'toolbox'));
import MAMaS.*
clc
if exist('maxima', 'var')
    maxima.stopInstance;
end

clear all

maxima = MaximaInterface.getInstance(4); %#ok<NASGU>

% Basic symbols
x = sym('x');
y = sym('y');
t = sym('t');
assert(x.isSymbol && ~x.isNumber && ~x.isMatrix, 'Expected x to be a symbol.');

% Numeric input (scalar)
n1 = sym(3.5, 'd');
assert(n1.isNumber, 'Expected numeric scalar to be a number.');
assert(strcmp(n1.identifier, num2str(3.5, 16)), 'Numeric identifier mismatch.');

% Numeric input (string)
n2 = sym("2");
assert(n2.isNumber, 'Expected numeric string to be a number.');

% Basic expression
expr = sin(x) + cos(y)^2 + sym.pi();
result = char(expr);
assert(ischar(result), 'Expected char output from sym.');
assert(~isempty(result), 'Expected non-empty output from sym.');

% Ensure constant factory works
c = sym.const('e');
assert(strcmp(c.identifier, '%e'), 'Expected identifier %e for constant e.');

% Matrix creation via constructor
M = sym([1, 2; 3, 4]);
assert(M.isMatrix, 'Expected matrix flag to be true.');
assert(M.matrixRows == 2 && M.matrixCols == 2, 'Matrix dimensions mismatch.');

% Matrix element access (subsref)
e12 = M(1, 2);
assert(isa(e12, 'sym'), 'Expected matrix element to be a sym.');

% end() indexing
e1end = M(1, end);
assert(isa(e1end, 'sym'), 'Expected end-indexed element to be a sym.');

% Subscript assignment (subsasgn)
M(1, 1) = 5;
e11 = M(1, 1);
e11str = char(e11);
assert(ischar(e11str) && ~isempty(e11str), 'Expected element to evaluate after assignment.');

% Matrix copy should create a different identifier
Mc = M.copy();
assert(Mc.isMatrix, 'Expected copy to be a matrix.');
assert(Mc.matrixRows == M.matrixRows && Mc.matrixCols == M.matrixCols, 'Copied matrix dimensions mismatch.');
assert(~strcmp(M.identifier, Mc.identifier), 'Matrix copy should have a different identifier.');

% Matrix exponential
Me = expm(M);
assert(Me.isMatrix, 'Expected expm result to be a matrix.');
assert(Me.matrixRows == M.matrixRows && Me.matrixCols == M.matrixCols, 'expm result dimensions mismatch.');

% Matrix transpose
Mt = M.';  % or M'
assert(Mt.isMatrix, 'Expected transpose result to be a matrix.');
assert(Mt.matrixRows == M.matrixCols && Mt.matrixCols == M.matrixRows, 'Transpose dimensions should be swapped.');

% Transpose using method
Mt2 = transpose(M);
assert(Mt2.isMatrix, 'Expected transpose() result to be a matrix.');

% Conjugate transpose
Mt3 = M';
assert(Mt3.isMatrix, 'Expected ctranspose result to be a matrix.');

% depends wrapper
depends(y, x);
depends(y, [x, t]);

% Test diff (differentiation)
fprintf('\nTesting diff()...\n');
expr_diff = x^3 + 2*x^2 + x;
dx = diff(expr_diff, x);
assert(isa(dx, 'sym'), 'Expected diff result to be sym.');
dxx = diff(expr_diff, x, 2);  % Second derivative
assert(isa(dxx, 'sym'), 'Expected second derivative to be sym.');
fprintf('  diff() tests passed!\n');

% Test jacobian
fprintf('Testing jacobian()...\n');
f1 = x^2 + y;
f2 = x*y + x;
jac = jacobian([f1, f2], [x, y]);
assert(jac.isMatrix, 'Expected jacobian result to be a matrix.');
assert(jac.matrixRows == 2 && jac.matrixCols == 2, 'Jacobian dimensions should be 2x2.');
fprintf('  jacobian() tests passed!\n');

% Test subs (substitution)
fprintf('Testing subs()...\n');
expr_subs = x^2 + y^2;
result_subs = subs(expr_subs, x, 2);
assert(isa(result_subs, 'sym'), 'Expected subs result to be sym.');
% Substitute multiple variables
result_subs2 = subs(expr_subs, [x, y], [2, 3]);
assert(isa(result_subs2, 'sym'), 'Expected multi-subs result to be sym.');
fprintf('  subs() tests passed!\n');

% Test symvar (find symbolic variables)
fprintf('Testing symvar()...\n');
expr_symvar = x^2 + y*sin(t) + 5;
vars = symvar(expr_symvar);
assert(isstring(vars) || ischar(vars), 'Expected symvar to return string array or char.');
% Should find x, y, t (and possibly sin as a function)
fprintf('  Found variables: ');
disp(vars);
fprintf('  symvar() tests passed!\n');

fprintf('\nAll tests passed!\n');

clear maxima;
