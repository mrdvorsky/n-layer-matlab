function [S11, S12, S21, S22] = coaxialToCircular(freq,epsilon,mu,sampleThk,a,b,d,numModes,numIntervals)
%written by Joseph Filbert
%this code implements a coaxial to circular to coaxial waveguide mode
%matching solution that can be used to optimize the unknown MUT properties
%to match the theoretical S parameters.
% a: coaxial inner conductor radius
% b: coaxial outer conductor radius
% d: outer radius circular waveguide
cc = 299792458; mu0 = 4*pi*1e-7; eps0 = 1/mu0/cc^2;
numFreqPts = length(freq);
chi = zeroBesselJ(0,numModes);
kc = coaxZero(numModes,a,b)/a;
skj = zeros(numModes,numModes,numFreqPts);
Ami = zeros(numModes,numModes);
Iji = zeros(numModes,numModes);
Bm = zeros(numModes,1);
ya = zeros(1,numModes);
yb = zeros(1,numModes);
bnTemp = zeros(numModes,1);
anjTemp = zeros(numModes,numModes);
firstRow = zeros(1,numModes);
iiVec = 1:numModes;
bj = zeros(numModes+1,numFreqPts);
ci = zeros(numModes,numFreqPts);
Ib = sameModeCrossProductIntegralCircularWaveguide(chi,d);
Ia = sameModeCrossProductIntegralCoaxialWaveguide(kc,a,b,numIntervals);
for jj=1:numModes
    for ii=1:numModes
        Iji(jj,ii) = modeCrossProductIntegral(chi(ii),kc(jj),a,b,d,numIntervals);
    end %for
end %for
for ww=1:numFreqPts
    for kk=1:numModes
        yb(kk) = admittanceCircularWaveguide(freq(ww),epsilon,mu,chi(kk),d);
        ya(kk) = admittanceCoaxialWaveguide(freq(ww),eps0*complex(1,0),mu0*complex(1,0),kc(kk));
    end %for
    for incidentMode=1:numModes
        for mm=1:numModes
            Bm(mm) = sum(yb(mm)*yb(incidentMode)*Iji(:,mm).*Iji(:,incidentMode)./transpose(ya.*Ia));
            for ii=1:numModes
                Ami(mm,ii) = sum(yb(mm)*yb(ii)*Iji(:,mm).*Iji(:,ii)./transpose(ya.*Ia));
            end %for
        end %for
        Ami = Ami+diag(yb.*Ib);
        Bm(incidentMode) = Bm(incidentMode)-yb(incidentMode)*Ib(incidentMode);
        skj(incidentMode,:,ww) = Ami\Bm;
    end %for
    incidentMode = 1;
    betaWaveguideB = betaCircularWaveguide(freq(ww),epsilon,mu,chi,d);
    phaseDelayMatrix = exp(-1i*transpose(betaWaveguideB)*sampleThk).*exp(-1i*betaWaveguideB*sampleThk);
    skjPrime = skj(:,:,ww).*phaseDelayMatrix;
    for nn=1:numModes
        bnTemp(nn) = -yb(nn)*Iji(incidentMode,nn);
        for jj=1:numModes
            firstRow(jj) = yb(jj)*Iji(incidentMode,jj)-sum(skjPrime(jj,:).*yb.*Iji(incidentMode,:));
            temp = 0;
            for ii=iiVec(1:end~=incidentMode)
                temp = temp+(yb(jj)*Iji(ii,jj)-sum(skjPrime(jj,:).*yb.*Iji(ii,:)))*yb(nn)*Iji(ii,nn)/Ia(ii)/ya(ii);
            end %for
            anjTemp(nn,jj) = temp+skjPrime(jj,nn)*yb(nn)*Ib(nn);
        end %for
    end %for
    Bn = [ya(incidentMode)*Ia(incidentMode);bnTemp];
    anjTemp = anjTemp+diag(yb.*Ib);
    firstColumn = transpose(yb.*Iji(incidentMode,:));
    Anj = [ya(incidentMode)*Ia(incidentMode),firstRow;...
           firstColumn,-anjTemp];
    bj(:,ww) = Anj\Bn;
    bjPrime = bj(2:end,ww).*transpose(exp(-1i*betaWaveguideB*sampleThk));
    for ii=1:numModes       
        temp = 0;
        for jj=1:numModes
            temp = temp+bjPrime(jj)*(yb(jj)*Iji(ii,jj)-sum(skj(jj,:,ww).*yb.*Iji(ii,:)));
        end %for
        ci(ii,ww) = temp/ya(ii)/Ia(ii);
    end %for
end %for

S11(:, 1) = bj(1, :);
S22 = S11;
S21(:, 1) = ci(1, :);
S12 = S21;
S = pagetranspose(reshape([S11;S12;S21;S22],[2 2 numFreqPts]));

%% Local Functions
function out = betaCircularWaveguide(freq,epsilon,mu,chi,d)
    out = sqrt(epsilon*mu*(2*pi*freq)^2-(chi/d).^2);
    out(imag(out)>0)=conj(out(imag(out)>0));
end %function

function out = coaxZero(numZeros,a,b)
%this gives kcTimesA
    beta =  [2.413373114376057, -0.016609217793955];
    %a = 3.0397/2*1e-3;
    %b = 7/2*1e-3;
    p = (1:numZeros);
    xInitial = beta(1)*p+beta(2);
    f = @(x) besselj(0,x).*bessely(0,x*b/a)-besselj(0,x*b/a).*bessely(0,x);
    out = zeros(1,numZeros);
    for ii=1:numZeros
        out(ii) = fzero(f,xInitial(ii));
    end %for
    out = [0,out(1:end-1)]; %chop off the extra zero
end %function

function out = modeCrossProductIntegral(chi,kc,a,b,d,numIntervals)
    %a = 3.0397/2*1e-3;
    %b = 7/2*1e-3;
    %d = 15/2*1e-3;
    h = (b-a)/numIntervals;
    xVec = (a+h/2):h:(b-h/2);
    chiTimesRho = xVec*chi;
    kcTimesRho = xVec*kc;
    kcTimesB = kc*b;
    factor = besselj(0,kcTimesB)./bessely(0,kcTimesB);
    if kc==0
        out = h*sum(-besselj(1,chiTimesRho/d));
    else
        out = h*sum((besselj(1,kcTimesRho)-factor.*bessely(1,kcTimesRho)).*besselj(1,chiTimesRho/d).*xVec);
    end %if
end %function

function out = sameModeCrossProductIntegralCircularWaveguide(chi,d)
    %d = 15/2*1e-3;
    out = (d^2/2)*(besselj(1,chi)).^2;
end %function

function out = sameModeCrossProductIntegralCoaxialWaveguide(kc,a,b,numIntervals)
    %a = 3.0397/2*1e-3;
    %b = 7/2*1e-3;
    kc = kc(2:end); %strip out kc=0
    N = length(kc);
    h = (b-a)/numIntervals;
    xVec = (a+h/2):h:(b-h/2);
    xVecMatrix = repmat(xVec,N,1);
    kcTimesRho = xVecMatrix.*kc.';
    kcTimesB = (kc*b).'; %Nx1
    factor = besselj(0,kcTimesB)./bessely(0,kcTimesB);
    out = h*(besselj(1,kcTimesRho)-factor.*bessely(1,kcTimesRho)).^2*xVec.';
    out = [log(b/a);out].'; %add back kc=0
end %function

function out = zeroBesselJ(n,numModes)
%halley's method
    out = zeros(1,numModes);
    for ii=1:numModes
        xInitial = 1+sqrt(2)+(ii-1)*pi+n+n^0.4;
        out(ii) = fzero(@(x) besselj(n,x),xInitial);
    end %for
end %function

function out = admittanceCircularWaveguide(freq,epsilon,mu,chi,d)
    %d = 15/2*1e-3;
    w = 2*pi*freq;
    beta = sqrt(epsilon*mu*w^2-(chi/d)^2);
    beta(imag(beta)>0) = conj(beta(imag(beta)>0));
    out = w*epsilon/beta;
end %for

function out = admittanceCoaxialWaveguide(freq,epsilon,mu,kc)
    if kc==0
        out = sqrt(epsilon/mu);
    else
        w = 2*pi*freq;
        beta = sqrt(epsilon*mu*w^2-(kc)^2);
        beta(imag(beta)>0) = conj(beta(imag(beta)>0));
        out = w*epsilon/beta;
    end %if
end %function

end %function