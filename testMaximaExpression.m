% Short test for MaximaExpression

maxima = MaximaInterface.getInstance();

x = MaximaExpression('x');
y = MaximaExpression('y');

expr = sin(x) + cos(y)^2 + MaximaExpression.pi();

% Force evaluation/printing
result = char(expr);

assert(ischar(result), 'Expected char output from MaximaExpression.');
assert(~isempty(result), 'Expected non-empty output from MaximaExpression.');

% Ensure constant factory works
c = MaximaExpression.const('e');
assert(strcmp(c.identifier, '%e'), 'Expected identifier %e for constant e.');

clear maxima;
