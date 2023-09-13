function [] = checkClassProperties(O, classProperties)
%CheckCLASSPROPERTIES Checks if properties of an object are settable.
% Checks if the input "classProperties", which is a cell array of
% string/value pairs, defines a valid set of properties for an object.
%
% Author: Matt Dvorsky

arguments
    O;
    classProperties(:, 1) {mustBeA(classProperties, "cell")};
end

%% Check Input Class Properties
if mod(numel(classProperties), 2) ~= 0
    error("Parameter and value arguments must come in pairs.");
end
for ii = 1:2:numel(classProperties)
    if ~isprop(O, classProperties{ii})
        error("The name '%s' is not an accessible property " + ...
            "for an instance of the class '%s'. Properties " + ...
            "must be one of the following: \n\n\t%s", ...
            classProperties{ii}, class(O), ...
            join(string(properties(O)), sprintf("\n\t")));
    end
end

end

