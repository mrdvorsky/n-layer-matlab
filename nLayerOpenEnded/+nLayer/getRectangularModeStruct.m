function [modeStruct] = getRectangularModeStruct(m, n, wgA, wgB, TE_TM)
%GETSPECTRUMRECTANGULAR Get function object defining waveguide spectrums.
% This function returns function objects

arguments
    m(1, 1) {mustBeNonnegative, mustBeInteger};
    n(1, 1) {mustBeNonnegative, mustBeInteger};
    wgA(1, 1) {mustBePositive};
    wgB(1, 1) {mustBePositive};
    TE_TM(1, 1) string {mustBeMember(TE_TM, ["TE", "TM"])};
end

%% Create Mode Spectrum Functions
kc = pi * hypot(m/wgA, n/wgB);
scale_All = pi ./ kc;
if strcmp(TE_TM, "TE")
    scaleX =  (n/wgB) .* scale_All;
    scaleY = -(m/wgA) .* scale_All;
else
    scaleX =  (m/wgA) .* scale_All;
    scaleY =  (n/wgB) .* scale_All;
end

Ex = @(kx, ky, ~, ~) scaleX .* rectSpectrum_cos(kx, wgA, m).*rectSpectrum_sin(ky, wgB, n);
Ey = @(kx, ky, ~, ~) scaleY .* rectSpectrum_sin(kx, wgA, m).*rectSpectrum_cos(ky, wgB, n);

%% Line Integrals
Nx = 1000;
Ny = 1000;

scale_All = 2 ./ kc ./ sqrt(wgA * wgB) ./ (2*pi) ...
    ./ sqrt(1 + (m == 0)) ./ sqrt(1 + (n == 0));

[xg(1, :), wx(1, :)] = fejer2(Nx, 0, wgA);
[yg(1, :), wy(1, :)] = fejer2(Ny, 0, wgB);

% [xg(1, :), wx(1, :)] = trap(Nx, 0, wgA);
% [yg(1, :), wy(1, :)] = trap(Ny, 0, wgB);

% [xg(1, :), wx(1, :)] = gaussLegendre(Nx, 0, wgA);
% [yg(1, :), wy(1, :)] = gaussLegendre(Ny, 0, wgB);

hertz_x = [  xg, 0*yg + wgA, flip(xg),           0*yg]  - 0.5*wgA;
hertz_y = [0*xg,   yg,          0*xg + wgB,   flip(yg)] - 0.5*wgB;
hertz_w = [wx, wy, flip(wx), flip(wy)];
hertz_ang = deg2rad([0*xg, 0*yg + 90, 0*xg + 180, 0*yg + 270]);

hertz_Ez = [];
hertz_Hz = [];
if strcmp(TE_TM, "TE")
    signTop   = -sign(mod(n, 2) - 0.5);
    signRight = -sign(mod(m, 2) - 0.5);
    hertz_Hz = scale_All * [...
             cos(m*pi/wgA * xg), ...
             cos(n*pi/wgB * yg) * signRight, ...
        flip(cos(m*pi/wgA * xg) * signTop), ...
        flip(cos(n*pi/wgB * yg))
        ];
else
    signTop   = sign(mod(n, 2) - 0.5);
    signRight = sign(mod(m, 2) - 0.5);
    hertz_Ez = scale_All * [...
             (n*pi/wgB) * sin(m*pi/wgA * xg), ...
             (m*pi/wgA) * sin(n*pi/wgB * yg) * signRight, ...
        flip((n*pi/wgB) * sin(m*pi/wgA * xg) * signTop), ...
        flip((m*pi/wgA) * sin(n*pi/wgB * yg))
        ];
end

%% Function Handles
% hertz_x = [-1,  1;  1, 1; 1, -1; -1,  1] .* 0.5 * wgA;
% hertz_y = [-1, -1; -1, 1; 1,  1;  1, -1] .* 0.5 * wgB;
% 
% if strcmp(TE_TM, "TE")
%     signTop   = -sign(mod(n, 2) - 0.5);
%     signRight = -sign(mod(m, 2) - 0.5);
%     hertz_Hz = scale_All * [...
%              cos(m*pi/wgA * xg), ...
%              cos(n*pi/wgB * yg) * signRight, ...
%         flip(cos(m*pi/wgA * xg) * signTop), ...
%         flip(cos(n*pi/wgB * yg))
%         ];
% else
%     signTop   = sign(mod(n, 2) - 0.5);
%     signRight = sign(mod(m, 2) - 0.5);
%     hertz_Ez = scale_All * {...
%              @(x) (n*pi/wgB) * sin(m*pi/wgA * xg), ...
%              (m*pi/wgA) * sin(n*pi/wgB * yg) * signRight, ...
%         flip((n*pi/wgB) * sin(m*pi/wgA * xg) * signTop), ...
%         flip((m*pi/wgA) * sin(n*pi/wgB * yg))
%         };
% end

%% Line Segments
Nx = m * 10;
Ny = n * 10;
pointsPerSeg = 5;
segOff = linspace(-0.5, 0.5, pointsPerSeg).';

scale_All = (-1j).^(m + n + 2) * -2 ./ kc ./ sqrt(wgA * wgB) ./ (2*pi) ...
    ./ sqrt(1 + (m == 0)) ./ sqrt(1 + (n == 0));

xg = ((1:Nx) - 0.5) ./ Nx;
yg = ((1:Nx) - 0.5) ./ Nx;

line_x = ([  xg, 0*yg + 1, flip(xg),         0*yg]  - 0.5) * wgA;
line_y = ([0*xg,   yg,        0*xg + 1,   flip(yg)] - 0.5) * wgB;

dx = wgA ./ Nx;
dy = wgB ./ Ny;
line_d = [0*xg + dx, 0*yg + dy, 0*xg + dx, 0*yg + dy];
line_ang = [0*xg + 0, 0*yg + 0.5*pi, 0*xg + pi, 0*yg - 0.5*pi];

line_Ez = [];
line_Hz = [];
if strcmp(TE_TM, "TE")
    signTop   = -sign(mod(n, 2) - 0.5);
    signRight = -sign(mod(m, 2) - 0.5);

    cosX = cos(xg(:).' + segOff(:));
    cosY = cos(yg(:).' + segOff(:));

    line_Hz = scale_All * [...
             cosX, ...
             cosY * signRight, ...
        flip(cosX * signTop, [2]), ...
        flip(cosY, [2])
        ];
else
    signTop   = sign(mod(n, 2) - 0.5);
    signRight = sign(mod(m, 2) - 0.5);

    sinX = (n*pi/wgB) * sin(m * pi * (xg(:).' + segOff(:) ./ Nx));
    sinY = (m*pi/wgA) * sin(n * pi * (yg(:).' + segOff(:) ./ Ny));

    line_Ez = scale_All * [...
                      sinX, ...
                      sinY * signRight, ...
            flip(flip(sinX * signTop, 2), 1), ...
            flip(flip(sinY, 2), 1)
        ].';
end

%% Create Mode Struct
symmetryY = "PMC";
if mod(m, 2) == 0
    symmetryY = "PEC";
end

symmetryX = "PMC";
if mod(n, 2) == 0
    symmetryX = "PEC";
end

modeStruct = nLayer.createModeStruct(TE_TM, ...
    sprintf("%s_{%d,%d}", TE_TM, m, n), ...
    ExSpec=Ex, EySpec=Ey, ...
    CutoffWavenumber=kc, ...
    ApertureWidth=hypot(wgA, wgB), ...
    SymmetryX=symmetryX, ...
    SymmetryY=symmetryY, ...
    Hertz_Ez=hertz_Ez, ...
    Hertz_Hz=hertz_Hz, ...
    Hertz_x=hertz_x, ...
    Hertz_y=hertz_y, ...
    Hertz_w=hertz_w, ...
    Hertz_ang=hertz_ang, ...
    Line_x=line_x, ...
    Line_y=line_y, ...
    Line_d=line_d, ...
    Line_ang=line_ang, ...
    Line_Hz=line_Hz, ...
    Line_Ez=line_Ez);

end




%% Integrals over Sin and Cos
function v = rectSpectrum_sin(k, a, m)
    if m == 0
        v = 0;
        return;
    end
    v = sqrt(a*pi) * (0.5/pi) ...
        .* (               sinc((0.5/pi) .* (a.*k - m.*pi)) ...
        + (-1).^(m + 1) .* sinc((0.5/pi) .* (a.*k + m.*pi)) );
end

function v = rectSpectrum_cos(k, a, m)
    if m == 0
        v = sqrt(0.5*a/pi) ...
            .* sinc((0.5/pi) .* (a.*k));
        return;
    end

    v = a .* sqrt(a/pi) .* (0.5/pi) .* k ./ m ...
        .* (               sinc((0.5/pi) .* (a.*k - m.*pi)) ...
        + (-1).^(m + 1) .* sinc((0.5/pi) .* (a.*k + m.*pi)) );
end




