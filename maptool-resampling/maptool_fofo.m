% maptool_fofo                              version 240724
% (C) 2024 Adams Vallejos
% 	https://doi.org/10.1063/4.000025
%
% -------------------------------------------------------------------------
% CLI USAGE
% matlab -nodisplay -nosplash -nodesktop

% INPUT

% -------------------------------------------------------------------------
clear
% clc

startTic = tic;

% SPHERE AND SIGMA SETTINGS
radius = 1.0; % Å
gridSpacing = 0.1; %Å, how dense grid within sphere
% rmsCutoff = 0.01; % exclude data below this sigma level
rmsCutoff = 2.0; % exclude data below this sigma level

% Files
% Main folder of MAPTOOL
mainPath = '/PATH/TO/resampling/output/';
dataSet = 'light/';
method = 'B/';
mapAnalysis = 'analysis_fofo_0';


mapsDirectory = [mainPath dataSet method'];


formatSpec= 'Working on: %s \n';
fprintf(formatSpec,mapsDirectory)

%%
for i = 1:100
    pdbFile{i} = sprintf('%03d/merge_trunc_cad_uniq_001_hetatm2atom.pdb', i-1);
    % mapFileName{i} = sprintf('%03d/merge_trunc_cad_uniq_fofo_map_coeffs', i-1);
    mapFileName{i} = sprintf('%03d/%s/merge_trunc_cad_uniq_fofo_map_coeffs', i-1, mapAnalysis);


    % Corresponding labels of maps above
    timePointLabels{i} = sprintf('%03d', i-1);
end


%  START CALCULATIONS
% ------------------------------------------------------------------------
% Number of atoms in pdb file
restingState = pdbread([mapsDirectory pdbFile{1}]);

% -------------------------------------------------------------------------
% Only ATOM label in pdb file, here HETATM is substituted by ATOM
atomicCoordinates = [...
    [restingState.Model.Atom.X]'...
    [restingState.Model.Atom.Y]'...
	[restingState.Model.Atom.Z]'...
];

% -------------------------------------------------------------------------

% Get number of atoms in the restingstate
numberOfAtoms = size(atomicCoordinates, 1);

% Calculate distance from sphere center
distanceWithinSphere = @(x,y,z) sqrt(x^2 + y^2 + z^2);

% -------------------------------------------------------------------------
% Generate spots of cubic grid with sphere inscribed
spots = linspace(-radius, radius, 2*(radius/gridSpacing) + 1);
numberOfSpots = length(spots);

% List grid spots within the sphere
sphereList= zeros(numberOfSpots^3, 7);
count = 0;
for i = 1 : numberOfSpots
   for j = 1 : numberOfSpots
        for k = 1 : numberOfSpots
            distance = distanceWithinSphere( spots(i) , spots(j) , spots(k) );
            if distance <= radius
               count = count + 1;
               sphereList(count,:) = [spots(i) spots(j) spots(k) i j k distance];
            end
        end
   end
end

sphereList = sphereList(1:count, :); 
numberOfPoints = size(sphereList, 1);
% -------------------------------------------------------------------------

% % PRECALCULATE ALL COORDINATES IN ALL SPHERES
% % (example: 1834 rows = atoms, 2109 columns = points) 
% X = repmat(atomicCoordinates(:,1), 1, numberOfPoints)+repmat(sphereList(:,1)', numberOfAtoms, 1);
% Y = repmat(atomicCoordinates(:,2), 1, numberOfPoints)+repmat(sphereList(:,2)', numberOfAtoms, 1);
% Z = repmat(atomicCoordinates(:,3), 1, numberOfPoints)+repmat(sphereList(:,3)', numberOfAtoms, 1);
% %
% % Reshape to a single column starting with first point for each atom, then second etc. 
% X=reshape(X, numberOfAtoms*numberOfPoints, 1);
% Y=reshape(Y, numberOfAtoms*numberOfPoints, 1);
% Z=reshape(Z, numberOfAtoms*numberOfPoints, 1);
%
% LOAD GRID FROM MAP
% (using first experimental XYZinfo file - here all files have the same grid)
XYZinfo = dlmread([mapsDirectory mapFileName{1} '_XYZinfo.dat']); 

%
% Dimensions of cartesian cell
cellDimensions = XYZinfo(1,1:3);


% Grid dimensions of cartesian map
gridPoints = XYZinfo(2,1:3); 
% Boundaries of cell axes50
axesBoundaries = XYZinfo(3,:);
% Order of axes
axesOrder = XYZinfo(4,1:3);

% Grid spacing
dX = cellDimensions(1)/gridPoints(1);
dY = cellDimensions(2)/gridPoints(2);
dZ = cellDimensions(3)/gridPoints(3);

% Re-arrange axes according to ordering of coordinates order from maps
axes = sortrows([axesOrder(1) axesBoundaries(1:2);
                 axesOrder(2) axesBoundaries(3:4);
                 axesOrder(3) axesBoundaries(5:6)]);

% Number of points on each side of the grid
sX = (axes(1,2) : axes(1,3)) * dX;
sY = (axes(2,2) : axes(2,3)) * dY;
sZ = (axes(3,2) : axes(3,3)) * dZ;

% 3D grid coordinates with dimensions from above
[gY, gZ, gX] = meshgrid(sY, sZ, sX);


% LOAD MAPS AND EXTRACT DENSITIES
% ------------------------------------------------------------------------
numberOfMaps = length(mapFileName);

% Read rms(sigma) from maps
rmsFromMaps = zeros(numberOfMaps, 2);

% Populate with map densities
normalizedDensities = zeros(numberOfMaps, numberOfAtoms, numberOfPoints);

for m = 1:numberOfMaps
    % Number of atoms in pdb file
    restingState = pdbread([mapsDirectory pdbFile{1}]);

    % -------------------------------------------------------------------------
    % Only ATOM label in pdb file, here HETATM is substituted by ATOM
    atomicCoordinates = [...
        [restingState.Model.Atom.X]'...
        [restingState.Model.Atom.Y]'...
	    [restingState.Model.Atom.Z]'...
    ];
    % PRECALCULATE ALL COORDINATES IN ALL SPHERES
    % (example: 1834 rows = atoms, 2109 columns = points)
    X = repmat(atomicCoordinates(:,1), 1, numberOfPoints)+repmat(sphereList(:,1)', numberOfAtoms, 1);
    Y = repmat(atomicCoordinates(:,2), 1, numberOfPoints)+repmat(sphereList(:,2)', numberOfAtoms, 1);
    Z = repmat(atomicCoordinates(:,3), 1, numberOfPoints)+repmat(sphereList(:,3)', numberOfAtoms, 1);
    %
    % Reshape to a single column starting with first point for each atom, then second etc.
    X=reshape(X, numberOfAtoms*numberOfPoints, 1);
    Y=reshape(Y, numberOfAtoms*numberOfPoints, 1);
    Z=reshape(Z, numberOfAtoms*numberOfPoints, 1);
    formatSpec= 'Reading %s\r';
    fprintf(formatSpec,mapFileName{m})

%   Load sigma info
    XYZinfo = dlmread([mapsDirectory mapFileName{m} '_XYZinfo.dat']);
    rmsFromMaps(m,1:2) = XYZinfo(5:6,1);

%   LOAD MAP AND CALCULATE DENSITY AT POINTS BY INTERPOLATION
%   reshape back to: rows = atoms, columns = points
    cartesianMap = h5read([mapsDirectory mapFileName{m} '_fullCell_cartesian.h5'],'/map');
    interpolatedDensities = interp3(gY,gZ,gX,cartesianMap,Y,Z,X);

    mapDensities = reshape(interpolatedDensities, numberOfAtoms, numberOfPoints);
%   Divide with sigma of original map
    normalizedDensities(m,:,:) = mapDensities/rmsFromMaps(m,2);
end
%%

% CALCULATE AVERAGE DENSITIES AND CORRELATIONS
%------------------------------------
% Calculate mean positive and negative densities
mapDensitiesAbsolute = normalizedDensities;
mapDensitiesAbsolute(abs(normalizedDensities) < rmsCutoff) = 0;

mapDensities_positive = mapDensitiesAbsolute;
mapDensities_positive(mapDensitiesAbsolute < 0) = 0;
meanpositiveDensities = mean(mapDensities_positive,3);

mapDensities_negative = mapDensitiesAbsolute;
mapDensities_negative(mapDensitiesAbsolute > 0) = 0;
meanNegativeDensities = mean(mapDensities_negative,3);

%% SAVE PEAKS
% save([mapsDirectory dataSet '_' mapAnalysis  '_plus.mat'], 'meanpositiveDensities' );
% save([mapsDirectory dataSet '_' mapAnalysis  '_minus.mat'], 'meanNegativeDensities' );


%% PREPARE PLOTS
%------------------------------------
fprintf('Preparing plots.\n')

% Plot settings
% RGB Colors normalized to 256
golden = [0.7, 0.5, 0];
purple = [0.5, 0.4, 1];
silver = [0.8, 0.8, 0.8];
fontSize = 12;
fontName = 'helvetica narrow';
lineWidth = 1; % Default value 0.5
labelCA = {'A', 'B', 'C', 'D', 'E', 'F','G'};
labelAll = {'A', 'B', 'C', 'D', 'E', 'F','G','RET'};

% Set fonts in all figures
set(0,'DefaultAxesFontName',fontName,'DefaultTextFontName',fontName);

% Avoid overlapping in plots        
scale_E = max(max(meanpositiveDensities-meanNegativeDensities));

% For plotting all atoms, find where helices start and stop
helix = [3   26;
         33  56;
         70  92;
         94 117;
         121 150;
         153 181;
         189 219;
         301 301];

% Get residues index from resting state
allResidues = [restingState.Model.Atom.resSeq]'; 

% Select residues from limits above
limitsHelix = zeros(size(helix));
for h = 1:size(helix,1)
    limitsHelix(h,1) = find(allResidues==helix(h,1), 1, 'first');
    limitsHelix(h,2) = find(allResidues==helix(h,2), 1, 'last');
end

% For plotting C alphas, find where helices start and stop
selectedCA = find(strcmp({restingState.Model.Atom.AtomName}','CA'));
residuesCA = [restingState.Model.Atom(selectedCA).resSeq]'; 

% PLOT
%------------------------------------
% figure('units','normalized','outerposition',[0 0 1 1],'name',['radius '...
%num2str(radius) ' Å, ' num2str(rmsCutoff) ' sigma,'])
figure('name',[dataSet  ' radius ' num2str(radius) ' Å, ' num2str(rmsCutoff) ' sigma,'])

% subplot(2,3,[2 3])
hold all
for n = 1:numberOfMaps
    plot(1:numberOfAtoms,...
         -meanpositiveDensities(n,:)/scale_E+(n),...
         'color', purple,...
         'LineWidth',lineWidth)
     
    plot(1:numberOfAtoms,...
         -meanNegativeDensities(n,:)/scale_E+(n),...
         'color', golden,...
         'LineWidth',lineWidth)
end

for i = 1:numberOfMaps
    line([1 numberOfAtoms], [(i) (i)],'color', silver)
    hold all
end 

for h = 1:length(limitsHelix)
    line([limitsHelix(h,1)-0.5 limitsHelix(h,1)-0.5],[0 numberOfMaps+1],'color', 'k')
    line([limitsHelix(h,2)+0.5 limitsHelix(h,2)+0.5],[0 numberOfMaps+1],'color', 'k')
end 
set(gca,'TickDir','out','Ytick',1:numberOfMaps,'Yticklabel', timePointLabels)
set(gca,'XTickLabel',labelAll,'XTick', mean(limitsHelix,2))
ylim([0 numberOfMaps+1])
xlim([1 numberOfAtoms])
title('All atoms (protein only)')
set(gca,'Ydir','reverse', 'XAxisLocation', 'top', 'box','on','FontSize',fontSize)

%% Plot averaged 
figure('name',[dataSet  ' radius ' num2str(radius) ' Å, ' num2str(rmsCutoff) ' sigma', 'Averaged'])
plot(mean(meanNegativeDensities))
hold on
plot(mean(meanpositiveDensities))

for h = 1:length(limitsHelix)
     line([limitsHelix(h,1)-0.5 limitsHelix(h,1)-0.5],[-1 +1],'color', 'k')
     line([limitsHelix(h,2)+0.5 limitsHelix(h,2)+0.5],[-1 +1],'color', 'k')
end 
%set(gca,'TickDir','out','Ytick',-1:+1,'Yticklabel', timePointLabels)
set(gca,'XTickLabel',labelAll,'XTick', mean(limitsHelix,2))
ylim([-0.2 +0.2])
xlim([1 numberOfAtoms])
title('Averaged resampling')
set(gca,'Ydir','reverse', 'XAxisLocation', 'top', 'box','on','FontSize',fontSize)


%saveas('Averaged values', 'Averaged_resampling', [pwd '/subFolderName/myFig.fig'])
%% JOB FINISHED
fprintf('\r\033[KJob finished! ET %.1f min\n', toc(startTic)/60)
% clear all
