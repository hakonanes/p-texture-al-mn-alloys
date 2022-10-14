%% Macrotexture from AA3xxx cold rolled to a true strain of 3.0 with a
% pronounced P texture component
%
% Created 2019-10-24 by Håkon Wiik Ånes (hakon.w.anes@ntnu.no)

clear all
close all
home

%% Import pole figure data and create PoleFigure object
cs = crystalSymmetry('m-3m', [4.04 4.04 4.04], 'mineral', 'Al');

to_plot = 1;

res = '-r100';

% Select sample
s0 = '00_191023_e3_asdef';
s0Fname = '3xxx_asdef_R';
%s1 = '01_200814_e3_50c';
%s1Fname = '3xxx_e3_50c_R';
s1 = '01_200909_e3_50c';
s1Fname = '3xxx_e3_50c_2_R';
s2 = '02_200817_e3_75c';
s2Fname = '3xxx_e3_75c_R';
s3 = '03_200818_e3_100c';
s3Fname = '3xxx_e3_100c_R';
s4 = '04_200818_e3_125c';
s4Fname = '3xxx_e3_125c_R';
s5 = '05_200818_e3_150c';
s5Fname = '3xxx_e3_150c_R';
s6 = '06_200819_e3_175c';
s6Fname = '3xxx_e3_175c_R';
s7 = '07_200819_e3_200c';
s7Fname = '3xxx_e3_200c_R';
s8 = '08_200819_e3_225c';
s8Fname = '3xxx_e3_225c_R';
s9 = '09_200827_e3_250c';
s9Fname = '3xxx_e3_250c_R';
s10 = '10_200911_e3_275c';
s10Fname = '3xxx_e3_275c_R';
s11 = '11_200911_e3_300c';
s11Fname = '3xxx_e3_300c_R';
s12 = '12_200914_e3_325c';
s12Fname = '3xxx_e3_325c_R';
s13 = '13_200914_e3_350c_2';
s13Fname = '3xxx_e3_350c_R';
s14 = '14_200915_e3_375c';
s14Fname = '3xxx_e3_375c_R';
s15 = '15_200915_e3_400c_2';
s15Fname = '3xxx_e3_400c_R';

sample = s15;
fnamesPrefix = s15Fname;
disp(sample);

basedir = '/home/hakon/phd/data/p/macrotexture';
inpath = fullfile(basedir, sample, 'texeval');
outpath = fullfile(basedir, sample, 'mtex');

fnames = {
    fullfile(inpath, [fnamesPrefix '_pf111_uncorr.dat']),...
    fullfile(inpath, [fnamesPrefix '_pf200_uncorr.dat']),...
    fullfile(inpath, [fnamesPrefix '_pf220_uncorr.dat']),...
    fullfile(inpath, [fnamesPrefix '_pf311_uncorr.dat'])};

% Specimen symmetry
ss = specimenSymmetry('1'); % Triclinic
ssO = specimenSymmetry('222'); % Orthorhombic

% Plotting convention
setMTEXpref('xAxisDirection', 'north');
setMTEXpref('zAxisDirection', 'outOfPlane');

% Set annotations to highlight spatial reference frame
pfAnnotations = @(varargin) text([vector3d.X, vector3d.Y],...
    {'RD', 'TD'}, 'BackgroundColor', 'w', 'tag', 'axesLabels', varargin{:});
setMTEXpref('pfAnnotations', pfAnnotations);

h = {Miller(1, 1, 1, cs), Miller(2, 0, 0, cs), Miller(2, 2, 0, cs),...
    Miller(3, 1, 1, cs)};

% Load pole figures separately
columnNames = {'Polar Angle', 'Azimuth Angle', 'Intensity'};
pf1 = loadPoleFigure_generic(fnames{1}, 'ColumnNames', columnNames);
pf2 = loadPoleFigure_generic(fnames{2}, 'ColumnNames', columnNames);
pf3 = loadPoleFigure_generic(fnames{3}, 'ColumnNames', columnNames);
pf4 = loadPoleFigure_generic(fnames{4}, 'ColumnNames', columnNames);

% Construct pole figure object of the four pole figures
intensities = {
    pf1.intensities,...
    pf2.intensities,...
    pf3.intensities,...
    pf4.intensities};
pfs = PoleFigure(h, pf1.r, intensities, cs, ss);

%% Plot pole figures of raw, corrected data

if to_plot
    figure;
    plot(pfs, 'upper', 'projection', 'eangle', 'minmax')
    mtexColorbar('location', 'southOutside')
    export_fig(fullfile(outpath, 'pfs_uncorr.png'), res);
end

%% Calculate the ODF using default settings
odf = calcODF(pfs, 'silent');

%% Set correct specimen symmetry for calculation of texture strength
odf.SS = ssO;
textureIndex = odf.textureindex
entropy = odf.entropy
odfMax = odf.max

%% Define ideal texture components and spread acceptance angle
br = orientation.byMiller([0 1 1], [2 -1 1], cs, ssO);
cu = orientation.byMiller([1 1 2], [1 1 -1], cs, ssO);
s = orientation.byMiller([1 2 3], [6 3 4], cs, ssO);
cube = orientation.byEuler(0, 0, 0, cs, ssO);
cubeND = orientation.byMiller([0 0 1], [3 1 0], cs, ssO);
p = orientation.byMiller([0 1 1], [-5 -6 6], cs, ssO);

% Component parameters
comps = {br, cu, s, cube, cubeND, p};
comp_colors = {'g', 'b', 'm', 'r', 'orange', 'c'};
comp_markers = {'d', '^', 'p', 's', 's', '>'};

%% ODF in PFs with specified contour levels

if to_plot
    levelsPF = [0, 1, 2, 3, 4, 5];

    odf.SS = ss;
    figure
    plotPDF(odf, h, 'upper', 'projection', 'eangle', 'contourf', levelsPF)
    mtexColorMap white2black
    mtexColorbar
    export_fig(fullfile(outpath, 'odf_pfs.png'), res);
end

%% ODF in PFs with specified contour levels
if to_plot
    levelsPF = [0, 1, 2, 3, 4, 5];

    odf.SS = ss;
    figure
    plotPDF(odf, h, 'upper', 'projection', 'eangle', 'contourf', levelsPF)
    mtexColorMap white2black
    hold on
    for i=1:length(comps)
        annotate(comps{i}.symmetrise, 'marker', comp_markers{i},...
            'markerfacecolor', comp_colors{i});
    end
    mtexColorbar('title', 'Multiples of Random Density (MRD)')
    hold off
    export_fig(fullfile(outpath, 'odf_pfs_annotated.png'), res);
end

%% Plot ODF in most relevant Euler space phi2 sections

if to_plot
    levelsODF = [0, 1, 2, 3, 4, 8, 12];

    odf.SS = ssO;
    figure
    plot(odf, 'phi2', [0 45 65]*degree, 'contourf', levelsODF)
    mtexColorMap white2black
    hold on
    for i=1:length(comps)
        annotate(comps{i}.symmetrise, 'marker', comp_markers{i},...
            'markerfacecolor', comp_colors{i});
    end
    mtexColorbar('title', 'Multiples of Random Density (MRD)')
    hold off
    export_fig(fullfile(outpath, 'odf_sections.png'), res);
end

%% Plot inverse pole figure

if to_plot
    figure
    plotIPDF(odf, [xvector, yvector, zvector], 'contourf') % contoured
    mtexColorMap white2black
    hold on
    for i=1:length(comps)
        annotate(comps{i}.symmetrise, 'marker', comp_markers{i},...
            'markerfacecolor', comp_colors{i}, 'markersize', 20);
    end
    mtexColorbar('title', 'Multiples of Random Density (MRD)')
    hold off
    export_fig(fullfile(outpath, 'odf_ipfs.png'), res);
end

%% Plot intensity along beta fibre from Cu to Brass and write results to
% file

cu_fiber = orientation.byEuler([90 35 45] * degree, cs, ssO);
br_fiber = orientation.byEuler([35 45 90] * degree, cs, ssO);

odf.SS = ssO;
f = fibre(cu_fiber, br_fiber, cs, ssO);

% generate list from fibres and evalute ODF at specific orientations
fibreOris = f.orientation;
evalOris = [];
evalIndex = [1 84 167 254 346 446 556 680 824 1000];
evalValues = zeros(1, 10);
for i=1:10
    ori = fibreOris(evalIndex(i));
    evalOris = [evalOris ori];
    evalValues(i) = eval(odf, ori);
end

if to_plot
    figure
    plot(evalOris.phi2/degree, evalValues, '-o')
    xlabel('\phi_2 \rightarrow', 'interpreter', 'tex')
    ylabel('Orientation density f(g)', 'interpreter', 'tex')
    xlim([45 90])
    export_fig(fullfile(outpath, 'fibre_beta.png'), res);
end

% Write fibre data to a csv file for further analysis
datafname = fullfile(outpath, 'data_fibre_beta.csv');

% Write header to file
fid = fopen(datafname, 'w');
fprintf(fid, '%s\r\n', 'phi1,Phi,phi2,fibreValue');
fclose(fid);

% Write Euler angles and intensities to file
dlmwrite(datafname, [(evalOris.phi1/degree)' (evalOris.Phi/degree)'...
    (evalOris.phi2/degree)' evalValues'], '-append')

%% Plot intensity along fibre from Cube to Goss and write results to file
f = fibre(cube, goss, cs, ssO);

% Generate list from fibres and evalute ODF at specific orientations
fibreOris = f.orientation;
evalOris = [];
evalIndex = [1 111 222 333 444 555 666 777 888 1000];
evalValues = zeros(1, 10);
for i=1:10
    ori = fibreOris(evalIndex(i));
    evalOris = [evalOris ori];
    evalValues(i) = eval(odf, ori);
end

if to_plot
    figure
    plot(evalOris.Phi/degree, evalValues, '-o')
    xlabel('\Phi \rightarrow', 'interpreter', 'tex')
    ylabel('Orientation density f(g)', 'interpreter', 'tex')
    xlim([min([f.o1.Phi f.o2.Phi]) max([f.o1.Phi f.o2.Phi])] / degree)
    export_fig(fullfile(outpath, 'fibre_cube_goss.png'), res);
end

% Write fibre data to a csv file for further analysis
datafname = fullfile(outpath, 'data_fibre_cube_goss.csv');

% Write header to file
fid = fopen(datafname, 'w');
fprintf(fid, '%s\r\n', 'phi1,Phi,phi2,fibreValue');
fclose(fid);

% Write Euler angles and intensities to file
dlmwrite(datafname, [(evalOris.phi1/degree)' (evalOris.Phi/degree)'...
    (evalOris.phi2/degree)' evalValues'], '-append')

%% plot intensity along fibre from Cube to ND-rotated Cube
f = fibre(cube, orientation.byEuler([45 0 0] * degree, cs, ssO), cs, ssO);

% Generate list from fibres and evalute ODF at specific orientations
fibreOris = f.orientation;
evalOris = [];
evalIndex = [1 111 222 333 444 555 666 777 888 1000];
evalValues = zeros(1, 10);
for i=1:10
    ori = fibreOris(evalIndex(i));
    evalOris = [evalOris ori];
    evalValues(i) = eval(odf, ori);
end

if to_plot
    figure
    plot(evalOris.phi1/degree, evalValues, '-o')
    xlabel('\phi_1 \rightarrow', 'interpreter', 'tex')
    ylabel('Orientation density f(g)', 'interpreter', 'tex')
    xlim([min([f.o1.phi1 f.o2.phi1]) max([f.o1.phi1 f.o2.phi1])] / degree)
    export_fig(fullfile(outpath, 'fibre_cube_cubend.png'), res);
end

% Write fibre data to csv file for further analysis in Python
datafname = fullfile(outpath, 'data_fibre_cube_cubeND45.csv');

% Write header to file
fid = fopen(datafname, 'w');
fprintf(fid, '%s\r\n', 'phi1,Phi,phi2,fibreValue');
fclose(fid);

% Write Euler angles and intensities to file
dlmwrite(datafname, [(evalOris.phi1/degree)' (evalOris.Phi/degree)'...
    (evalOris.phi2/degree)' evalValues'], '-append')

%% Plot intensity along fibre from Goss to P
goss_fiber = orientation.byEuler([0 45 0] * degree, cs, ssO);
f = fibre(goss_fiber, p, cs, ssO);

% Generate list from fibres and evalute ODF at specific orientations
fibreOris = f.orientation;
evalOris = [];
evalIndex = [1 111 222 333 444 555 666 777 888 1000];
evalValues = zeros(1, 10);
for i=1:10
    ori = fibreOris(evalIndex(i));
    evalOris = [evalOris ori];
    evalValues(i) = eval(odf, ori);
end

if to_plot
    figure
    plot(evalOris.phi1/degree, evalValues, '-o')
    xlabel('\phi_1 \rightarrow', 'interpreter', 'tex')
    ylabel('Orientation density f(g)', 'interpreter', 'tex')
    xlim([min([f.o1.phi1 f.o2.phi1]) max([f.o1.phi1 f.o2.phi1])] / degree)
    export_fig(fullfile(outpath, 'fibre_goss_p.png'), res);
end

% Write fibre data to csv file for further analysis in Python
datafname = fullfile(outpath, 'data_fibre_goss_p.csv');

% Write header to file
fid = fopen(datafname, 'w');
fprintf(fid, '%s\r\n', 'phi1,Phi,phi2,fibreValue');
fclose(fid);

% Write Euler angles and intensities to file
dlmwrite(datafname, [(evalOris.phi1/degree)' (evalOris.Phi/degree)'...
    (evalOris.phi2/degree)' evalValues'], '-append')

%%
close all

%% Calculate volume fractions Mi
odf.SS = ssO;
spread = 15*degree;

Mbr = 100*volume(odf, br, spread)
Mcu = 100*volume(odf, cu, spread)
Ms = 100*volume(odf, s, spread)
Mcube = 100*volume(odf, cube, spread)
McubeND = 100*volume(odf, cubeND, spread)
Mp = 100*volume(odf, p, spread)