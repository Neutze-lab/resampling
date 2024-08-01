clear
close all


pdb = pdbread('/home/arpitha/bm07_processing_april_24/100K/anom_ref/Anomalous_refmac1.pdb');


margin = 5;
gridSpacing = 0.25;

X = [pdb.Model.Atom(:).X]';
Y = [pdb.Model.Atom(:).Y]';
Z = [pdb.Model.Atom(:).Z]';

pdbLimits = [min(X) max(X) min(Y) max(Y) min(Z) max(Z)];

% Limits with border
limits = ceil([pdbLimits(1)-margin pdbLimits(2)+margin pdbLimits(3)-margin pdbLimits(4)+margin pdbLimits(5)-margin pdbLimits(6)+margin]);

% Cell parameters
CELLWORK = [limits(2)-limits(1) limits(4)-limits(3) limits(6)-limits(5) 90 90 90];

%
GRIDWORK = CELLWORK(1:3)*1/gridSpacing;

XYZLIM = limits/0.25;

disp(CELLWORK)
disp(GRIDWORK)
disp(XYZLIM)

%%%
% 0.016 
% CELL WORK 43 55 78 90 90 90
% GRID WORK 172 220 312
% XYZLIM -144 28 -180 40 -136 176
% 0.760 & 1725
% CELL WORK 44 55 78 90 90 90
% GRID WORK 176 220 312
% XYZLIM -148 28 -180 40 -136 176
