% isscalar MATLAB object
function tf = isscalar_matlab(obj)
    tf = prod(builtin('size', obj))==1;
end

