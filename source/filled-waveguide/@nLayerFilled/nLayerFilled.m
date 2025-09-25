classdef nLayerFilled < nLayerForward
    %Implementation of nLayerForward for filled waveguides. 
    % This class can be used to calculate the reflection and transmission
    % for a MUT in a filled waveguide measurement configuration. 
    % This class is meant to be subclassed, where the subclass defines the 
    % mode to be considered and defines properties of the waveguide
    % (cutoff, k, kz, etc.).
    %
    % Note: The units of all parameters should match that of the speed of
    % light specified by the "speedOfLight", "distanceUnitScale", and
    % "frequencyUnitScale" parameters (defaults are mm and GHz).
    %
    % Author: Trent Moritz

    %% Public Properties
    properties (GetAccess=public, SetAccess=protected)
        mode_kc0(1,1) {mustBeNonnegative, mustBeFinite} = 0; % Free space cutoff wave number (default plane wave)
        modeType(1,1) {mustBeMember(modeType,["TE","TM"])} = "TE"; % Type of mode ("TE", "TM").
    end
    
    properties (Dependent, Access=public, AbortSet)
        frequencyRange(1, :) {mustBeNonnegative, mustBeFinite}; % Operating frequency range of the object.
    end

    properties (Access=public)
        waveguidePort1er(1,1) {nLayer.mustBeErUrCallable} = 1; %Permittivity of port 1 filling. 
        waveguidePort1ur(1,1) {nLayer.mustBeErUrCallable} = 1; %Permeability of port 1 filling
        waveguidePort2er(1,1) {nLayer.mustBeErUrCallable} = 1; %Permittivity of port 2 filling.
        waveguidePort2ur(1,1) {nLayer.mustBeErUrCallable} = 1; %Permeability of port 2 filling. 
    end

    properties (Dependent, GetAccess=public)
        mode_fc0;       % Free space cutoff frequency. 
    end

    %% Class Functions
    methods (Access=public)
        [Smn] = calculate(self,f,er,ur,thk);
        [outputLabels] = getOutputLabels(self);
    end

    methods (Access=protected)
        

    end

    %% Class Constructor
    methods (Access=public)
        function self = nLayerFilled(mode_kc0, modeType, classProperties)
            arguments
                 mode_kc0(1,1) {mustBeNonnegative, mustBeFinite} = 0; 
                 modeType(1,1) {mustBeMember(modeType,["TE","TM"])} = "TE";
                 classProperties.?nLayerFilled;
            end

            if ~isempty(mode_kc0)
                self.mode_kc0 = mode_kc0;
            end

            if ~isempty(modeType)
                self.modeType = modeType;
            end

            % Set Class Parameter Values
            propPairs = namedargs2cell(classProperties);
            for ii = 1:2:numel(propPairs)
                self.(propPairs{ii}) = propPairs{ii + 1};
            end
        end
    end

    %% Class Setters
    methods
        function set.frequencyRange(self,newFreqRange)
            self.frequencyRange_private = [min(newFreqRange),max(newFreqRange)];
        end

        function set.mode_kc0(self,newModeKc0)
            self.mode_kc0 = newModeKc0;
        end

    end

    %% Class Getters
    methods
        function [freq_range] = get.frequencyRange(self)
            freq_range = self.frequencyRange_private;
        end

        function [fc0] = get.mode_fc0(self)
            fc0 = self.speedOfLight * self.mode_kc0 ./ (2*pi);
        end
    end
end