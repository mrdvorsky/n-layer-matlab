NL1 = nLayerRectangular(3, 2, Band='x'); 
NL2 = nLayerRectangular(5, 4, Band='x'); 
NL3 = nLayerRectangular(1, 0, Band='x');
% NL(4) = nLayerRectangular(maxM, maxN, Band=wgBand); 
% NL(5) = nLayerRectangular(maxM, maxN, Band=wgBand); 
% NL(6) = nLayerRectangular(maxM, maxN, Band=wgBand); 

f1 = linspace(8.2, 12.4, 21).';  
f2 = linspace(8.2, 12.4, 21).';  
f3 = linspace(8.2, 12.4, 21).';  

er = [11-11j, 2-0.2j];
thk = [1, 2];

legendCell = ["1", "2", "3"];

nLayerViewer(er, thk, NL1, f1, NL2, f2);
