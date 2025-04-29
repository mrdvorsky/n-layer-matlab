function [] = cndeMatlabTools_startup(options)
%Function that is run at Matlab startup.
% 
% Author: Matt Dvorsky

arguments
    options.CheckForUpdates(1, 1) logical = true;
end

%% Inputs and Paths
% The line below will be dynamically replaced on install.
libPath = "<LIBRARY_PATH>";
if contains(libPath, "<") && contains(libPath, ">")
    error("String replacement was not completed.");
end

%% Add to Path
addpath(fullfile(libPath, "install"));
addpath(genpath(fullfile(libPath, "source")));

%% Check for Updates
if options.CheckForUpdates
    cndeMatlabTools_checkForUpdates(AlwaysShowPopupWindow=false);
end

end

