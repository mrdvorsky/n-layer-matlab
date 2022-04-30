function [] = printStructureSetup(O, options)
%PRINTSTRUCTURESETUP Summary of this function goes here
%   Detailed explanation goes here

arguments
    O;
    options.ShowLimits(1, 1) {mustBeNumericOrLogical} = true;
end

%% Options
layerWidth = 60;
tableItemWidth = 20;
erFormat = "%g";
erpFormat = "%g";
thkFormat = "%g";

%% Create Helper Strings
sepString = string(repmat('_', 1, layerWidth));
backString = string(repmat('/', 1, layerWidth));

%% Print Structure
fprintf("\t\t\t\tMultilayer Structure Stackup\n");
for ii = 1:O.layerCount
    % Layer Separator
    fprintf("%s\n\n", sepString);
    
    % Value Specifiers
    thkString = sprintf("thk = %s mm", ...
        sprintf(thkFormat, O.thkGuess(ii)));
    
    erString = sprintf("er = %s - j%s", ...
        sprintf(erFormat, O.erGuess(ii)), ...
        sprintf(erpFormat, O.erpGuess(ii)));
    
    urString = "ur = 1";
    
    fprintf("\t%s", strjoin(pad([thkString, erString, urString], ...
        tableItemWidth)));
    fprintf("\n");
    
    if options.ShowLimits
        % Lower Limit Specifier
        thkLowerString = sprintf(" [%s,", ...
            sprintf(thkFormat, O.thkRange(1, ii)));
        
        erLowerString = sprintf(" [%s - j%s,", ...
            sprintf(erFormat, O.erRange(1, ii)), ...
            sprintf(erpFormat, O.erpRange(1, ii)));
        
        urLowerString = " [1,";
        
        fprintf("\t%s", strjoin(pad([thkLowerString, erLowerString, ...
            urLowerString], tableItemWidth)));
        fprintf("\n");
        
        % Upper Limit Specifier
        thkUpperString = sprintf("  %s]", ...
            sprintf(thkFormat, O.thkRange(2, ii)));
        
        erUpperString = sprintf("  %s - j%s]", ...
            sprintf(erFormat, O.erRange(2, ii)), ...
            sprintf(erpFormat, O.erpRange(2, ii)));
        
        urUpperString = "  1]";
        
        fprintf("\t%s", strjoin(pad([thkUpperString, erUpperString, ...
            urUpperString], tableItemWidth)));
        fprintf("\n");
    end
end

%% Backing Layer
fprintf("%s\n", sepString);
fprintf("%s\n", backString);

end
