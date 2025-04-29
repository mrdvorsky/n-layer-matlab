function [] = cndeMatlabTools_checkForUpdates(options)
%This function "updates" the cnde-matlab-tools library.
% Essentially, this function does the following:
%   - Fetch latest from remote git repo.
%   - Check if there is an update.
%   - If there is an update, ask the user if they want to update.
%       - If user has modified libary, ask if they want to revert all
%           changes, and then do so.
%       - Perform the update.
%       - Re-install (fix paths, and then copy update and startup files).
%
% Author: Matt Dvorsky

arguments
    options.AlwaysShowPopupWindow(1, 1) logical = true;
end

%% Inputs and Paths
% The line below will be dynamically replaced on install.
libPath = "<LIBRARY_PATH>";
if contains(libPath, "<") && contains(libPath, ">")
    error("String replacement was not completed.");
end

[~, libName] = fileparts(libPath);
updaterTitle = sprintf("'%s' Updater", libName);

%% Look for Git Repo
try
    repo = gitrepo(libPath);
catch ex
    warning("Failed to open git repository. Git may not be installed.");
    rethrow(ex);
end
fetch(repo);

%% Check Alternate Branch
if ~strcmp(repo.CurrentBranch.Name, "main")
    if options.AlwaysShowPopupWindow
        msgbox(sprintf("The '%s' libary is currently switched to " + ...
            "a different branch (%s). Update will not be performed. " + ...
            "Make sure to run the update function once you have " + ...
            "switched back.", ...
            libName, repo.CurrentBranch.Name), ...
            updaterTitle);
    end
    return;
end

%% Check for Remote Updates
updateLog = repo.log(Revisions="origin/main");

isThereAnUpdate = ~strcmp(updateLog.ID(1), ...
    repo.CurrentBranch.LastCommit.ID);

if ~isThereAnUpdate
    if options.AlwaysShowPopupWindow
        msgbox(sprintf("The '%s' libary is up to date (%s).", ...
            libName, cndeMatlabTools_getVersion()), ...
            updaterTitle);
    end
    return;
end

%% Ask to Update
resp = questdlg(sprintf(...
    "The '%s' library has updates. Do you want to update?", ...
    libName), ...
    updaterTitle, ...
    "Update", "Cancel", "Update");

if strcmp(resp, "Cancel")
    return;
end

%% Check for Local Changes to Library
changesTable = repo.status(IncludeIgnoredFiles=false, ...
    IncludeUntrackedFiles=false);

if height(changesTable) > 0
    resp = questdlg(sprintf(...
        "The '%s' library has been modified, and thus updating " + ...
        "will result in any changes being lost. Do you want to proceed?", ...
        libName), ...
        updaterTitle, ...
        "?? Revert all Changes ??", "Cancel", "Cancel");
    if strcmp(resp, "Cancel")
        return;
    end

    for ii = 1:height(changesTable)
        repo.discardChanges(changesTable.File{ii});
    end
end

%% Perform Update
pull(repo);
msgbox(sprintf("The '%s' libary was successfully updated to %s.", ...
    libName, cndeMatlabTools_getVersion()), ...
    updaterTitle);

%% Setup Paths by Reinstalling
cndeMatlabTools_install(ShowFinishedPopup=false);

end


