function [gam] = calculate_impl(O, f, er, ur, thk)
%CALCULATEGAMMA Calculate S11 for clamped waveguides. "Clamping" can be
%done with any of the open-ended nLayerForward objects.
% 
% STRUCTURE FORMAT: When defining er, ur, and thk for each layer start with
%   port 1 and list them in order going towards port 2. 
%
% Inputs:
%   f - Column vector of frequencies.
%   er - Array of complex relative permittivities for each layer. Every row
%       er(nL, :) contains the permittivity of each layer. Alternatively
%       you can define er for each frequency as er(nF, nL).
%   ur - Same as er, except for complex relative permeability.
%   thk - Vector of thicknesses for each layer. The length of thk is the
%       same as the number of columns in er and ur. Last element can be
%       inf to represent an infinite half-space.
% Outputs:
%   gam - Equivalent S-parameters for the clamped waveguide in the form of
%   s(nf, 2, 2).
%
% Author: Trent Moritz

arguments
    O;
    f(:,1);
    er(1, :);
    ur(1,:);
    thk(1,:);
end

magCond = O.magCond;
layers = numel(thk);
NL = O.NL;

%% Format parameters for PEC and PMC calculations
if mod(layers, 2) ~=0 %Odd number of layers, split middle layer
    er_pec = er(1:(layers+1)/2);
    ur_pec = ur( 1:(layers+1)/2);
    thk_pec = [thk(1:(layers+1)/2-1), thk{(layers+1)/2}/2];
    er_pmc = [er(1:(layers+1)/2), 1];
    ur_pmc = [ur( 1:(layers+1)/2), 1-1j*magCond];
    thk_pmc = [thk(1:(layers+1)/2-1), thk{(layers+1)/2}/2, inf];
else
    error("nLayerClamped has not been implemented for even layered objects yet.")
    %% If even number of layers it MUST be symmetric on either side
    % er_pec = er(1,1:(layers)/2);
    % ur_pec = ur(1, 1:(layers)/2);
    % thk_pec = thk(1:(layers)/2);
    % er_pmc = [er(1, 1:(layers)/2), 1];
    % ur_pmc = [ur(1, 1:(layers)/2), 1-1j*magCond];
    % thk_pmc = [thk(1:(layers)/2), inf];
end

pec = NL.calculate(f, er_pec, ur_pec, thk_pec);
pmc = NL.calculate(f, er_pmc, ur_pmc, thk_pmc);
s11 = (pec+pmc)/2;
s21 = (pmc-pec)/2;


%% Set Output
gam = reshape([s11, s21, s21, s11], length(f), 2,2); 

if ~isempty(O.outputIndices)
    gam = gam(:, O.outputIndices);
end

end

