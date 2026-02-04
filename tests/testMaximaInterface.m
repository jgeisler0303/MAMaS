clc
% Create Maxima interface
maxima = MaximaInterface();

% Mathematische Operationen
res1 = maxima.send('factor(x^2 - 1);');
disp(res1);

% Variablen speichern
maxima.send('a: 5; b: 7;');
res2 = maxima.send('a*b + 2;');
disp(res2);

% Maxima instance is automatically terminated on close
clear maxima;
