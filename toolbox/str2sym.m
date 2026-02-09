function s = str2sym(name)
    % STR2SYM Create a symbolic msym from string
    % Wrapper for MATLAB Symbolic Toolbox compatibility with function notation
    %
    % s = str2sym('x') - creates a msym 'x'
    % s = str2sym('phi(time)') - creates 'phi' with dependency on 'time'
    % s = str2sym('x(a,b,c)') - creates 'x' with dependencies on a, b, c
    
    name = char(name);
    
    % Check if there are parentheses indicating dependencies
    parenIdx = strfind(name, '(');
    
    if isempty(parenIdx)
        % Simple symbol without dependencies
        s = msym(name);
    else
        % Extract symbol name and dependencies
        symName = name(1:parenIdx(1)-1);
        
        % Extract content between parentheses
        closeParenIdx = strfind(name, ')');
        if isempty(closeParenIdx) || closeParenIdx(end) < parenIdx(1)
            error('Mismatched parentheses in symbol name: %s', name);
        end
        
        depsStr = name(parenIdx(1)+1:closeParenIdx(end)-1);
        
        % Parse dependencies (comma-separated)
        if isempty(strtrim(depsStr))
            % Empty parentheses - just create symbol
            s = msym(symName);
        else
            % Split by comma and create dependency symbols
            depNames = strsplit(depsStr, ',');
            depNames = strtrim(depNames); % Remove whitespace
            
            % Create the main symbol
            s = msym(symName);
            
            % Create dependency symbols
            deps = cell(size(depNames));
            for i = 1:length(depNames)
                deps{i} = msym(depNames{i});
            end
            
            % Declare dependencies
            depends(s, deps);
        end
    end
end
