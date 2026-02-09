% isvector of MATLAB objects (not Maxima matrix)
function tf = isvector_matlab(obj)
    % Vector if MATLAB array with one dimension being 1
    % Also returns true for scalars (1x1 MATLAB objects)
    sz = builtin('size', obj);
    if prod(sz) == 1
        % Scalar is considered a vector
        tf = true;
    elseif prod(sz) == 0
        tf = false; % Empty array is not a vector
    else
        % Multiple elements - check if one dimension is 1
        tf = length(sz)==2 && (sz(1) == 1 || sz(2) == 1);
    end
end
