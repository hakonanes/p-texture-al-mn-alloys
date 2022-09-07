% Estimate geometrically necessary disloations prior to grain and texture
% analysis
%
% Håkon Wiik Ånes (hakon.w.anes@ntnu.no)
% 2022-09-01

% MTEX configuration
plotx2east
plotzIntoPlane

% Crystal and specimen symmetry
cs = {'notIndexed',...
    crystalSymmetry('m-3m', [4.04 4.04 4.04], 'mineral', 'al')};
cs_al = cs{2};

% Directory and file names
sample = '0s';
dset_no = '2';
dir_data = fullfile('/home/hakon/phd/data/p/prover', sample, dset_no);
dir_kp = fullfile(dir_data, 'kp');
dir_mtex = fullfile(dir_data, 'mtex');
fname_ori = 'xmap_refori2.ang';

%% Read orientation data
%ebsd = EBSD.load(fullfile(dir_kp, fname_ori), cs, 'columnNames', ...
%    {'phi1', 'Phi', 'phi2', 'x', 'y', 'iq', 'ci', 'Phase',...
%    'detector_signal', 'fit', 'n_particles', 'n_pixels'}, 'radians');
%
%% Align kikuchipy's to MTEX' crystal reference frame
%rot_tsl2mtex = rotation.byAxisAngle(xvector - yvector, 180 * degree);
%ebsd = rotate(ebsd, rot_tsl2mtex, 'keepXY');
%
%% Get gridify shape in
%size_gridify_in = size(ebsd.gridify);
%
%% Remove unnecessary fields
%ebsd.prop = rmfield(ebsd.prop, {'detector_signal', 'ci', 'iq', 'fit'});
%
%% Fix not getting non-indexed points from .ang file properly
%ebsd.phaseMap(1) = -1;
%
%% Set phase of particles to 'notIndexed'
%ebsd(ebsd.n_particles > 0).phase = -1;
%
%% Misorientation angle threshold (mat) for grain reconstruction
%mat = 1 * degree;
%
%% Grain reconstruction
%[grains, ebsd.grainId, ebsd.mis2mean] = calcGrains(ebsd, 'angle',...
%    mat, 'boundary', 'tight', 'unitCell');
%
%% Assign small Al grains to surrounding grains after a second
%% reconstruction
%ebsd2 = ebsd(grains(grains.grainSize < 5));
%ebsd2 = ebsd2('al');
%ebsd(ismember(ebsd.id, ebsd2.id)) = [];
%[grains2, ebsd.grainId, ebsd.mis2mean] = calcGrains(ebsd, 'angle',...
%    mat, 'boundary', 'tight', 'unitCell');
%
%% Add equivalent circular diameter (ECD)
%grains2.prop.ecd = 0.816 * 2 * grains2.equivalentRadius;
%
%% Smooth grain boundaries
%grains2 = smooth(grains2, 5);
%
%% Remove unnecessary variables
%clear grains ebsd2
%
%% Gridify
%ebsd = ebsd.gridify;
%
%% Estimate geometrically necessary dislocations
%% Define FCC dislocation system
%dS = dislocationSystem.fcc(cs_al);
%
%% Line edge energies
%dS(dS.isEdge).u = 1;
%dS(dS.isScrew).u = 1 - 0.347;
%
%% Get lattice curvature
%kappa = ebsd.curvature;
%
%% Rotate the dislocation tensors into the specimen reference frame
%dSRot = ebsd.orientations * dS;
%
%% Save intermediate results
%save(fullfile(dir_mtex, 'kappa'), 'kappa');
%save(fullfile(dir_mtex, 'dSRot'), 'dSRot');

% Load intermediate results
kappa = load(fullfile(dir_mtex, 'kappa.mat'));
dSRot = load(fullfile(dir_mtex, 'dSRot.mat'));
kappa = kappa.kappa;
dSRot = dSRot.dSRot;

% Fit dislocation tensor to the dislocation density tensor in each pixel
[rho, factor] = fitDislocationSystems(kappa, dSRot);

% Total dislocation energy
gnd = factor * sum(abs(rho .* dSRot.u), 2);

% Write GND density to file
save(fullfile(dir_mtex, 'gnd.mat'), 'gnd');

% Write figure of GND density to file
%figure
%plot(ebsd, gnd)
%mtexColorMap LaboTeX
%mtexColorbar('title', 'GND')
%caxis([1e13 1e15])
%hold on
%plot(grains2.boundary)
%plot(ebsd('notIndexed'), 'facecolor', 'k')
%legend('hide')
%export_fig(fullfile(dir_mtex, 'gnd.png'), '-r100')
