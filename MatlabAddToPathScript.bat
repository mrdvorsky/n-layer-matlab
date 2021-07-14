@ECHO OFF
SETLOCAL ENABLEDELAYEDEXPANSION

SET startupHeader=%%%% Amntl Matlab Utilities

SET scriptPath=%~dp0
SET startupPath=C:\Users\%username%\Documents\MATLAB\startup.m
SET copyString=addpath(genpath("%scriptPath%"));


IF EXIST %startupPath% (
	FIND /c "addpath(genpath(""%scriptPath%""))" %startupPath% > NUL
	IF !errorlevel! EQU 0 (
		ECHO Already installed.
	) ELSE (
		CALL :performInstallation
	)
) ELSE (
	CALL :performInstallation
)

PAUSE
EXIT /B %errorlevel%
:performInstallation
	ECHO.>>%startupPath%
	ECHO.>>%startupPath%
	ECHO %startupHeader%>>%startupPath%
	ECHO %copyString%>>%startupPath%
	ECHO | SET /p="%copyString%" | clip
	ECHO Installation successful. Restart Matlab or run the following command (automatically copied to clipboard): %copyString%