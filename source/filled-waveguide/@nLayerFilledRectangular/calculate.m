function [gam] = calculate(self, f, er, ur, thk)
%CALCULATEGAMMA Calculate S11 for rectangular waveguide TEmn mode excitation.
% Computes the reflection coefficient of the rectangular waveguide TE10
% mode when looking into a multilayer structure defined by er, ur, thk at
% the frequencies defined by f.
% 
% STRUCTURE FORMAT: When defining er, ur, and thk for each layer start with
% port 1 and list them in order going towards port 2. 
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
%   gam - Equivalent S-parameters for the waveguide in the form of s(2, 2, nf).
%
% Author: Trent Moritz

arguments
    self nLayerFilledRectangular;

    f(:, 1);
    er;
    ur;
    thk;
end

%% Check Inputs
[er, ur, thk] = nLayer.validateStructure(er, ur, thk, ...
    CheckStructureValues=false);

%% Fix
er = cell2mat(er).';
ur = cell2mat(ur).';
thk = cell2mat(thk).';

%% Define Variables
a = self.waveguideA;
b = self.waveguideB;
m = self.modeTE_m;
n = self.modeTE_n;
c = self.speedOfLight;

k0 = (2*pi) .* f./c;
eta0 = 120*pi;
eta = eta0 .* sqrt(ur./er);
k = k0 .* sqrt(ur.*er);
kc = hypot(m.*pi./a, n.*pi./b);
beta = sqrt(k.^2 - kc.^2);
z = k.*eta./beta;
z_ref = k0.*eta0 ./ sqrt(k0.^2 - kc.^2);

%% Create ABCD Matrix
A = cos(beta.*thk);
B = (1j.*z) .* sin(beta.*thk);
C = (1j./z) .* sin(beta.*thk);
D = cos(beta.*thk);
abcd_final = eye(2, 2);
abcd_layers = [permute(A, [3,4,1,2]), permute(B, [3,4,1,2]); ... 
               permute(C, [3,4,1,2]), permute(D, [3,4,1,2])];

for ii = (1:length(thk))
    abcd_final = pagemtimes(abcd_final, abcd_layers(:, :, :, ii));
end

abcd_final = ipermute(abcd_final, [3,4,1,2]);

%% Redefine Matrix Elements for Conversion to S Parameters. 
A_final = abcd_final(:, :, 1, 1);
B_final = abcd_final(:, :, 1, 2);
C_final = abcd_final(:, :, 2, 1);
D_final = abcd_final(:, :, 2, 2);

%% Convert ABCD Parameters to S Parameters
S11 = (A_final + B_final./z_ref - C_final.*z_ref - D_final) ./ ...
      (A_final + B_final./z_ref + C_final.*z_ref + D_final); 
S12 = (2.*(A_final.*D_final - B_final.*C_final)) ./ ...
      (A_final + B_final./z_ref + C_final.*z_ref + D_final);
S21 = (2) ./ ...
      (A_final + B_final./z_ref + C_final.*z_ref + D_final);
S22 = (-A_final + B_final./z_ref - C_final.*z_ref + D_final) ./ ...
      (A_final + B_final./z_ref + C_final.*z_ref + D_final);

%% Set Output
gam = reshape([S11, S12, S21, S22], length(f), 2, 2);

if ~isempty(self.outputIndices)
    gam = gam(:, self.outputIndices);
end

end

