function [] = nLayerMatlab_uninstall()
%This function "uninstalls" the n-layer-matlab library.
% Essentially, this function does the following:
%   - Removes the line "nLayerMatlab_startup(...)" from the
%       "startup.m" file, if the line is present.
%   - Deletes the "nLayerMatlab_startup.m" file, if it exists.
%   - Deletes the "nLayerMatlab_checkForUpdates.m" file, if it exists.
%   - Remove library folders from path.
%
% Author: Matt Dvorsky

%% Inputs and Paths
startupLineSearch = "nLayerMatlab_startup(";

libStartupFileName = "nLayerMatlab_startup.m";
libUpdateFileName = "nLayerMatlab_checkForUpdates.m";

libName = "n-layer-matlab";
libPath = fileparts(fileparts(mfilename("fullpath")));

startupFileLocation = userpath();
startupFilePath = fullfile(startupFileLocation, "startup.m");

%% Check if "startup.m" Exists at UserPath
if ~isempty(dir(startupFilePath))
    % "startup.m" exists. Search for relevant lines.
    startupFileTextLines = readlines(startupFilePath);
    isStartupLine = startsWith(startupFileTextLines, ...
        startupLineSearch);

    if any(isStartupLine)
        writelines(startupFileTextLines(~isStartupLine), ...
            startupFilePath, ...
            WriteMode="overwrite");
    end
end

%% Clean Path
pathAll = split(path(), pathsep());
pathAllLib = pathAll(startsWith(pathAll, libPath));
if ~isempty(pathAllLib)
    rmpath(pathAllLib{:});
end

%% Delete the Functions in UserPath Folder
delete(fullfile(startupFileLocation, strcat(libStartupFileName, "*")));
delete(fullfile(startupFileLocation, strcat(libUpdateFileName, "*")));

%% Popup Window
resp = questdlg(sprintf(...
    "Library '%s' was successfully uninstalled. The git repository " + ...
    "is still present on your computer and must be deleted manually. " + ...
    "Do you want to change directories so that you can delete it?", ...
    libName), ...
    sprintf("'%s' Uninstaller", libName), ...
    "Yes", "No", "Yes");

if strcmp(resp, "Yes")
    cd(fileparts(libPath));
end

end

