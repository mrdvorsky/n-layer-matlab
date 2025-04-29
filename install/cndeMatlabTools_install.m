function [] = cndeMatlabTools_install(options)
%This function "installs" the cnde-matlab-tools library.
% Essentially, this function does the following:
%   - Creates the "startup.m" file in the userpath directory, if it does
%       not already exist.
%   - Adds the line "cndeMatlabTools_startup(...);" to the "startup.m"
%       file, if it does not already exist.
%   - Copies "cndeMatlabTools_startup.m" to the userpath directory.
%   - Copies "cndeMatlabTools_checkForUpdates.m" to the userpath directory.
%   - Runs the "cndeMatlabTools_startup" function.
%
% Author: Matt Dvorsky

arguments
    options.ShowFinishedPopup(1, 1) logical = true;
end

%% Inputs and Paths
startupLineSearch = "cndeMatlabTools_startup(";
startupLineFull   = "cndeMatlabTools_startup(CheckForUpdates=true);";

libStartupFileName = "cndeMatlabTools_startup.m";
libUpdateFileName = "cndeMatlabTools_checkForUpdates.m";

libName = "cnde-matlab-tools";
libPath = fileparts(fileparts(mfilename("fullpath")));

startupFileLocation = userpath();
startupFilePath = fullfile(startupFileLocation, "startup.m");

%% Check if "startup.m" Exists at UserPath
if isempty(dir(startupFilePath))
    % Create empty file.
    writelines("", startupFilePath);
end

% "startup.m" exists. Search for relevant line.
startupFileTextLines = readlines(startupFilePath);
isStartupLine = startsWith(startupFileTextLines, ...
    startupLineSearch);

lineDoesNotExist = ~any(isStartupLine);
if lineDoesNotExist
    writelines(startupLineFull, ...
        startupFilePath, ...
        WriteMode="append");
end

%% Copy Files to UserPath Folder
libStartupFilePath = fullfile(libPath, ...
    "install", "installer-files", ...
    libStartupFileName);
libStartupLines = readlines(libStartupFilePath);
libStartupLines = strrep(libStartupLines, "<LIBRARY_PATH>", libPath);
writelines(libStartupLines, fullfile(startupFileLocation, libStartupFileName));


libUpdateFilePath = fullfile(libPath, ...
    "install", "installer-files", ...
    libUpdateFileName);
libUpdateLines = readlines(libUpdateFilePath);
libUpdateLines = strrep(libUpdateLines, "<LIBRARY_PATH>", libPath);
writelines(libUpdateLines, fullfile(startupFileLocation, libUpdateFileName));

%% Clean Path
pathAll = split(path(), pathsep());
pathAllLib = pathAll(startsWith(pathAll, libPath));
rmpath(pathAllLib{:});

%% Run Startup
cndeMatlabTools_startup(CheckForUpdates=false);

%% Finished Dialogue
if options.ShowFinishedPopup
    msgbox(sprintf("Library '%s' was successfully installed.", libName), ...
        sprintf("'%s Installer'", libName));
end

end

