% Correlated BSE/EBSD analysis of Al and particles in an Al-Mn alloy
%
% Håkon Wiik Ånes (hakon.w.anes@ntnu.no), 2022-11-24
% Norwegian University of Science and Technology (NTNU)

%clear variables
%close all

% Generate figures?
to_plot = 1;
plot_sanity_check = 0;

% MTEX configuration
plotx2east
plotzIntoPlane

% export_fig configuration
res = '-r200';

% Crystal and specimen symmetry
cs = {'notIndexed',...
    crystalSymmetry('m-3m', [4.04 4.04 4.04], 'mineral', 'al')};
cs_al = cs{2};
ssO = specimenSymmetry('orthorhombic');

% Directory and file names
sample = '0s'; % 0s, 175c, 300c, 325c
dset_no = '3';
dir_data = fullfile('/home/hakon/phd/data/p/prover', sample, dset_no);
disp(dir_data)
dir_kp = fullfile(dir_data, 'kp');
dir_mtex = fullfile(dir_data, 'mtex');
fname_ori = 'xmap_refori2.ang';

% Read orientation data
ebsd = EBSD.load(fullfile(dir_kp, fname_ori), cs, 'columnNames', ...
    {'phi1', 'Phi', 'phi2', 'x', 'y', 'iq', 'ci', 'Phase',...
    'detector_signal', 'fit', 'n_particles', 'n_pixels'}, 'radians');

% Align kikuchipy's to MTEX' crystal orientation reference frame
rot_tsl2mtex = rotation.byAxisAngle(xvector - yvector, 180 * degree);
ebsd = rotate(ebsd, rot_tsl2mtex, 'keepXY');

% Remove unnecessary properties
ebsd.prop = rmfield(ebsd.prop, {'detector_signal', 'ci', 'iq', 'fit'});

% Step sizes
dx = ebsd.gridify.dx; % EBSD
dx_bse = 0.025; % Upscaled BSE image pixel size

% Fix not getting non-indexed points from .ang file properly
ebsd.phaseMap(1) = -1;

% Set phase of particles to 'notIndexed'
ebsd(ebsd.n_particles > 0).phase = -1;

% Misorientation angle threshold (mat) for grain reconstruction
mat = 1 * degree;

% How high is a high angle grain boundary (HAB)?
hab = 15; % Degrees

% Particle classification threshold in microns
dispersoid_threshold_min = 0.03;
dispersoid_threshold_max = 0.24;
constituent_particle_threshold = 0.8;

% Orientation color keys
om_al = ipfHSVKey(cs_al);
om_al.CS2 = ssO;

%% (Sub)grain reconstruction
ebsd2 = ebsd;

[grains, ebsd2.grainId, ebsd2.mis2mean] = calcGrains(ebsd2, 'angle',...
    mat, 'boundary', 'tight', 'unitCell');

% Remove small Al (sub)grains
ebsd3 = ebsd2(grains(grains.grainSize < 5));
ebsd3 = ebsd3('al');
ebsd2(ismember(ebsd2.id, ebsd3.id)) = [];
[grains2, ebsd2.grainId, ebsd2.mis2mean] = calcGrains(ebsd2, 'angle',...
    mat, 'boundary', 'tight', 'unitCell');

% Add equivalent circular diameter (ECD)
grains2.prop.ecd = 0.816 * 2 * grains2.equivalentRadius;

% Add actual size (ECD) of detected particles from unbinned particle map.
% This is zero for Al subgrains.
n_particle_pixels = grainMean(ebsd2, ebsd2.n_pixels, grains2, @sum);
grains2.prop.particle_ecd = 0.816 * 2 * dx_bse *...
    sqrt(n_particle_pixels / pi);

% Smooth subgrain boundaries
grains2 = smooth(grains2, 5);

% Remove unnecessary variables
clear grains ebsd3

%% Boundaries and misorientation angles
gb2 = grains2.boundary;
mori = gb2.misorientation.angle ./ degree;
bin_edges = [mat hab max(mori)];
[~, ~, gb2Id] = histcounts(mori, 'NumBins', 2, 'BinEdges', bin_edges);

%% Plot orientation maps
if to_plot
    directions = {xvector, yvector, zvector};
    titles = {'nd', 'rd', 'td'};
    for i=1:length(directions)
        om_al.inversePoleFigureDirection = directions{i};
        figure
        plot(ebsd('al'), om_al.orientation2color(ebsd('al').orientations),...
            'micronBar', 'off')
        hold on
        plot(ebsd('notIndexed'), 'facecolor', 'k')
        legend('hide')
        hold off
        pause(1)
        export_fig(fullfile(dir_mtex, ['maps_om_ipf_' titles{i} '.png']), res)
%        close(gcf)
    end
end

%% Plot RD orientation map with grain boundaries overlayed
if to_plot
    om_al.inversePoleFigureDirection = yvector;
    figure
    plot(ebsd('al'), om_al.orientation2color(ebsd('al').orientations),...
        'micronBar', 'off')
    hold on
    plot(ebsd('notIndexed'), 'facecolor', 'k')
    plot(gb2(gb2Id == 1), 'linecolor', [0.7 0.7 0.7], 'linewidth', 1)
    plot(gb2(gb2Id == 2), 'linecolor', [0 0 0], 'linewidth', 1);
    legend('hide')
    hold off
    export_fig(fullfile(dir_mtex, 'maps_om_ipf_rd_gb.png'), res)
    close(gcf)
end

%% Ideal orientations
% Ideal grain texture components, rotated to match the global reference
% frame:
%   X (east)         = ND
%   Y (north)        = RD
%   Z (out of plane) = TD
rot_scan2global = rotation.byMatrix([0 0 1; 1 0 0; 0 1 0]);

br = rot_scan2global * orientation.byMiller([0 1 1], [2 -1 1], cs_al, ssO);
cu = rot_scan2global * orientation.byMiller([1 1 2], [1 1 -1], cs_al, ssO);
s = rot_scan2global * orientation.byMiller([1 2 3], [6 3 4], cs_al, ssO);
cube = rot_scan2global * orientation.byEuler(0, 0, 0, cs_al, ssO);
cubend = rot_scan2global * orientation.byMiller([0 0 1], [3 1 0], cs_al, ssO);
p = rot_scan2global * orientation.byMiller([0 1 1], [-5 -6 6], cs_al, ssO);

ideal_oris = {br, cu, s, cube, cubend, p};
ideal_colors = {'m', 'b', 'g', 'r', [1 0.55 0], 'c'};
ideal_markers = {'d', '^', 'p', 's', 's', '>'};
ideal_oris_labels = {'br', 'cu',  's', 'cube', 'cubend', 'p'};
n_ideal = length(ideal_oris);

%% Plot inverse pole figure key with components annotated
if 0
    om_al.inversePoleFigureDirection = yvector; % RD
    figure
    plot(om_al)
    hold on
    ms = 20;
    for i=1:length(ideal_oris)
        annotate(ideal_oris{i}, 'marker', ideal_markers{i}, 'markersize',...
            ms, 'markerfacecolor', ideal_colors{i});
    end
    export_fig(fullfile(dir_mtex, 'ipf_annotated.png'), '-r600')
    close(figure)
end

%% Assign a texture component to each grain (without overlap)
mori_threshold_deg = 15;

% Mean orientations of grains
grain_oris = grains2.meanOrientation;
grain_oris.SS = ssO;
n_grains = length(grains2);

% Get misorientation angle between each grain and component
grain_ideal_mangle = zeros(n_ideal, n_grains);
for i=1:n_ideal
    grain_ideal_mangle(i, :) = angle(grain_oris, ideal_oris{i}) / degree;
end

% Get index of minimum misorientation angle
[~, idx] = nanmin(grain_ideal_mangle);

% Assign ideal orientation and ideal orientation ID to each grain
default_ori = orientation.byEuler(1 * degree, 0, 0, cs_al, ssO);
ideal_ori = repmat(default_ori, n_grains, 1);
ideal_ori_id = zeros(n_grains, 1);
for i=1:length(grains2)
    idx_i = idx(i);
    mangle = grain_ideal_mangle(idx_i, i);
    if ~isnan(mangle) && mangle <= mori_threshold_deg
        ideal_ori(i) = ideal_oris{idx_i};
        ideal_ori_id(i) = idx_i;
    end
end
grains2.prop.ideal_ori = ideal_ori;
grains2.prop.ideal_ori_id = ideal_ori_id;

%% Plot orientation map of higher fidelity to illustrate
% a multimodal dataset (from dataset #1 at 300 C)
%region = [40, 35, 15, 15];

% EBSD ROI
%ebsd_roi = ebsd(inpolygon(ebsd, region));

% Grains within ROI
%grains_roi = grains2(inpolygon(grains2, region));
%particles_roi = grains_roi('notIndexed');
%particles_const_roi = particles_roi(...
%    particles_roi.particle_ecd >= constituent_particle_threshold);
%particles_disp_roi = particles_roi(...
%    (particles_roi.particle_ecd >= dispersoid_threshold_min) &...
%    (particles_roi.particle_ecd <= dispersoid_threshold_max)...
%);

% Subgrain boundaries within ROI
%gb_roi = gb2(inpolygon(gb2, region));
%mori_roi = gb_roi.misorientation.angle ./ degree;
%bin_edges_roi = [mat hab max(mori_roi)];
%[~, ~, gb_roiId] = histcounts(mori_roi, 'NumBins', 2, 'BinEdges',...
%    bin_edges_roi);

%figure
%for i=1:n_ideal
%    grains_i = grains_roi(ismember(grains_roi.ideal_ori_id, i));
%    if ~isempty(grains_i)
%        if i < 4
%            color = 'g';
%        elseif i == 7
%            continue
%        else
%            color = ideal_colors{i};
%        end
%        plot(grains_i, 'facecolor', color, 'micronBar', 'off', 'facealpha', 0.5)
%    end
%    hold on
%end
%plot(particles_roi, 'facecolor', 'k')
%hold on
%plot(particles_const_roi.boundary, 'linecolor', 'r', 'linewidth', 4)
%hold on
%plot(particles_disp_roi, 'facecolor', 'r')
%hold on
%plot(gb_roi(gb_roiId == 1), 'linecolor', [0.7 0.7 0.7], 'linewidth', 2)
%hold on
%plot(gb_roi(gb_roiId == 2), 'linecolor', [0 0 0], 'linewidth', 2)
%legend('hide')
%export_fig(fullfile(dir_mtex, 'maps_roi_grains_ideal_particles.png'), '-r300')

%% Plot of grains with grain boundaries per component
if to_plot
    figure
    for i=1:n_ideal
        grains_i = grains2(ismember(grains2.ideal_ori_id, i));
        if ~isempty(grains_i)
            plot(grains_i, 'facecolor', ideal_colors{i}, 'micronBar', 'off')
        end
        hold on
    end
    plot(ebsd('notIndexed'), 'facecolor', 'k')
    legend('hide')
    export_fig(fullfile(dir_mtex, 'maps_grains_ideal_particles.png'), res)

    pause(1)

    plot(gb2(gb2Id == 1), 'linecolor', [0.7 0.7 0.7], 'linewidth', 1)
    plot(gb2(gb2Id == 2), 'linecolor', [0 0 0], 'linewidth', 1)
    export_fig(fullfile(dir_mtex, 'maps_grains_ideal_particles_gb.png'), res)

    close(gcf)
end

%% ---------------------------------- DISPERSOIDS CLOSE TO GRAIN BOUNDARIES
distance_threshold = dx; % um
distance_considered = 2; % um. Sufficiently large, typically ECD dependent
pad = round(distance_considered / dx);

% Grain boundaries of Al-Al not on the map edges
gb2_al = gb2('al', 'al');
gb2_al = gb2_al(~any(gb2_al.grainId == 0, 2));
gb2_idx = 1:size(gb2_al);

% Extract dispersoids to loop over
dispersoid_condition = (grains2.phase == -1) & ...
    (grains2.particle_ecd <= dispersoid_threshold_max) & ...
    (grains2.particle_ecd >= dispersoid_threshold_min);
grains_particles = grains2(dispersoid_condition);

% Number of dispersoids
n = length(grains_particles);

% Indices of boundaries each particle is within the distance threshold to.
% 300 deemed a sufficiently large number.
boundary_idx = repmat(zeros(1, 600, 'int64'), [n, 1]);

% Keep track of minimum distance to grain boundary for each particle
grains2.prop.min_distance_to_gb = -ones(size(grains2));

h = waitbar(0, 'Finding boundary segments particles are close to');
for i=1:n
    waitbar(i / n)

    particle = grains_particles(i);

    % Extract coordinates of all particle boundary segment
    gb_particle = particle.boundary;
    gb_particle_midpoint = gb_particle.midPoint;
    gb_particle_midpoint_x = gb_particle_midpoint(:, 1);
    gb_particle_midpoint_y = gb_particle_midpoint(:, 2);

    % Calculate extent around particle
    x_min = min(gb_particle_midpoint_x);
    x_max = max(gb_particle_midpoint_x);
    x_extent = round((x_max - x_min) / dx);
    y_min = min(gb_particle_midpoint_y);
    y_max = max(gb_particle_midpoint_y);
    y_extent = round((y_max - y_min) / dx);
    extent = [x_min - distance_considered, y_min - distance_considered,...
        (x_extent + 2 * pad - 1) * dx, (y_extent + 2 * pad - 1) * dx];

    % Extract coordinates of grain boundary segments of interest
    % surrounding particle
    roi = inpolygon(gb2_al, extent);
    gb_of_interest = gb2_al(logical(roi));
    gb_of_interest_midpoint = gb_of_interest.midPoint;
    gb_of_interest_midpoint_x = gb_of_interest_midpoint(:, 1);
    gb_of_interest_midpoint_y = gb_of_interest_midpoint(:, 2);

    % Calculate distance to particle boundary coordinates for each grain
    % boundary coordinate
    distances = sqrt(...
        (gb_particle_midpoint_x - gb_of_interest_midpoint_x').^2 +...
        (gb_particle_midpoint_y - gb_of_interest_midpoint_y').^2);

    % Extract minimum distances for all grain boundary segments
    min_distance = squeeze(min(distances, [], 1))';

    % Extract indices of all boundaries within ROI within distance
    % threshold
    mask = logical(min_distance <= distance_threshold);
    idx_of_interest = gb2_idx(roi);
    boundary_idx(i, 1:sum(mask)) = idx_of_interest(mask);

    % Assign minimum distance to a grain boundary to the particle
    if ~isempty(min_distance)
        grains2(particle.id).min_distance_to_gb = min(min_distance);
    end
end

close(h)

% Trim array of grain boundary indices to non-zero elements by finding
% the index of max non-zero boundary index across all particles
[~, c] = find(boundary_idx);
max_c = max(c);
boundary_idx2 = boundary_idx(:, 1:max_c);

% Fill sub-array for dispersoids into full array for all grains
boundary_idx_all = repmat(zeros(1, max_c), [length(grains2), 1]);
boundary_idx_all(dispersoid_condition, :) = boundary_idx2;

% Assign boundary indices to all grains
grains2.prop.boundary_idx = boundary_idx_all;

% Assign number of dispersoid particles per boundary segment
boundary_idx_nz = nonzeros(grains2.boundary_idx);
[unique_boundary_idx, ~, ic] = unique(boundary_idx_nz);
particle_counts_present = accumarray(ic, 1);
particle_counts_all = zeros(size(gb2));
particle_counts_all(ismember(gb2_idx, unique_boundary_idx)) =...
    particle_counts_present;
gb2.prop.n_particles_close = particle_counts_all;

%% Assign unbinned size (ECD) of dispersoid particles per boundary segment
particles_close_size = zeros(size(gb2_al)); % Al GBs
particle_sizes = grains_particles.particle_ecd; % Unbinned particle ECD

for i=1:n % Loop over dispersoids
    gb_idx_i = nonzeros(boundary_idx2(i, :));
    if ~isempty(gb_idx_i)
        for j=1:size(gb_idx_i)
            k = gb_idx_i(j);
            particles_close_size(k) = particles_close_size(k) +...
                particle_sizes(i);
        end
    end
end

particles_close_size_all = zeros(size(gb2)); % All GBs
particles_close_size_all(gb2.isIndexed) = particles_close_size;
gb2.prop.particles_close_size = particles_close_size_all;

%% Get fraction of HAB for each grain
h = waitbar(0, 'Calculating fraction of HAB around each grain');
grain_ids = grains2('indexed').id;
numGrains = length(grains2('indexed'));
xhab = zeros(numGrains, 1);
for i = 1:numGrains
    waitbar(i/numGrains)
    id = grain_ids(i);
    % Create logical vector
    try
        isHAB = grains2(id).boundary('indexed', 'indexed').misorientation...
            .angle./degree > hab;
        % Create new property of HAB fraction (nnz = number of non-zero)
        xhab(i) = nnz(isHAB)/length(isHAB);
    catch ME
        fprintf('Grain %i is not here!\n', id);
    end
end
close(h)

xhab_all = zeros(length(grains2), 1);
xhab_all(grains2.isIndexed) = xhab;
grains2.prop.xhab = xhab_all;

%% Set whether a grain is recrystallized
%grains2_data = csvread(fullfile(dir_mtex, 'grains.txt'), 1, 0);
%grains2.prop.xhab = grains2_data(:, 11);

grains2.prop.is_rx = zeros(length(grains2), 1);
grains2(...
    grains2.phase==1 &...
    grains2.ecd > 4 &...  % 3 is too low
    grains2.GOS < 1 * degree &...  % 0.5 too low
    grains2.xhab > 0.5...  % 0.75 is too high
).is_rx = 1;

if to_plot
    figure
    plot(grains2('indexed'), grains2('indexed').is_rx)
    hold on
    plot(grains2('notIndexed'), 'FaceColor', 'k')
    hold off
    export_fig(fullfile(dir_mtex, 'maps_grains_is_rx.png'), res)
end

%% Set whether grain boundaries are of recrystallized grains
rx_ids = grains2(grains2.is_rx == 1).id;
gb2.prop.is_rx = ismember(gb2.grainId(:, 1), rx_ids) |...
    ismember(gb2.grainId(:, 2), rx_ids);

%% Plots as a sanity check of particle locations and boundaries of interest
gb2_al = gb2('al', 'al');
gb3 = gb2_al(ismember(gb2_idx, grains2.boundary_idx));

if plot_sanity_check
    % Check that grain boundary indices are correct by correlating particle
    % locations with boundary locations
    figure
    plot(grains_particles, 'facecolor', 'r')
    hold on
    plot(gb3('al', 'al'), 'linewidth', 2)

    % Check that number of particles per boundary is correct
    figure
    plot(grains_particles, 'facecolor', 'r')
    hold on
    plot(gb2_al, gb2_al.n_particles_close)

    % Check assigned minimum distance to grain boundary by ensuring that
    % particles on boundaries and particles within boundaries have correct
    % hues
    figure
    plot(grains2('notIndexed'), grains2('notIndexed').min_distance_to_gb)
    hold on
    plot(gb2_al)
end

%% ---------- SPECIAL BOUNDARIES: 40 deg<111> (CSL7) and 60 deg<111> (CSL3)
% Number of dispersoids per boundary of interest
gb_component1 = zeros([1, length(gb2)]);
gb_component2 = zeros([1, length(gb2)]);
gb2_ids = gb2.grainId;
gb2_ids1 = gb2_ids(:, 1);
gb2_ids2 = gb2_ids(:, 2);
for i=1:n_ideal
    comp_ids = grains2(grains2.ideal_ori == ideal_oris{i}).id;
    gb_component1(ismember(gb2_ids1, comp_ids)) = i;
    gb_component2(ismember(gb2_ids2, comp_ids)) = i;
end
gb2.prop.component1 = gb_component1;
gb2.prop.component2 = gb_component2;

% Special boundary misorientation
gb2.prop.is_csl3 = angle(CSL(3, cs_al), gb2.misorientation) <= 15 * degree;
gb2.prop.is_csl7 = angle(CSL(7, cs_al), gb2.misorientation) <= 15 * degree;

%% Store misorientation of all boundaries and boundaries with particles
gb2_al = gb2('al', 'al');
export(gb2_al.misorientation, fullfile(dir_mtex, 'mori_gb_all.txt'), 'quaternion', 'radians');
mori_particles = gb2_al(gb2_al.n_particles_close > 0).misorientation;
export(mori_particles, fullfile(dir_mtex, 'mori_gb_with_particles.txt'), 'quaternion', 'radians');

%% Plot special boundaries
if to_plot
    pause(1)
    figure
    for i=1:n_ideal
        grains_i = grains2(ismember(grains2.ideal_ori_id, i));
        if ~isempty(grains_i)
            plot(grains_i, 'facecolor', ideal_colors{i}, 'facealpha', 0.25)
        end
        hold on
    end
    plot(ebsd2('notIndexed'), 'facecolor', 'k', 'micronBar', 'off')
    plot(gb2, 'linecolor', [0.5, 0.5, 0.5], 'linewidth', 1)
    plot(gb2(gb2.is_csl7), 'linecolor', 'r', 'linewidth', 1)
    legend('hide')
    export_fig(fullfile(dir_mtex, 'grains_ideal_csl7.png'), res)
    close(figure)

    figure
    for i=1:n_ideal
        grains_i = grains2(ismember(grains2.ideal_ori_id, i));
        if ~isempty(grains_i)
            plot(grains_i, 'facecolor', ideal_colors{i}, 'facealpha', 0.25)
        end
        hold on
    end
    plot(ebsd2('notIndexed'), 'facecolor', 'k', 'micronBar', 'off')
    plot(gb2, 'linecolor', [0.5, 0.5, 0.5], 'linewidth', 1)
    plot(gb2(gb2.is_csl3), 'linecolor', 'r', 'linewidth', 1)
    legend('hide')
    export_fig(fullfile(dir_mtex, 'grains_ideal_csl3.png'), res)
    close(figure)
end

%% ----------------- MISORIENTATIONS OF GRAINS AROUND CONSTITUENT PARTICLES
grains_particle = grains2('notIndexed');

% Whether Al grains neighbor a particle
pairs1 = grains_particle.neighbors('full');
pairs1 = unique(pairs1);
grains2.prop.at_particle = (ismember(grains2.id, pairs1)) &...
    (grains2.phase == 1);

% Whether Al grains neighbor a constituent particle
grains_particle_constituent = grains_particle(...
    grains_particle.particle_ecd >= constituent_particle_threshold);
pairs1 = grains_particle_constituent.neighbors('full');
pairs1 = unique(pairs1);
grains2.prop.at_constituent_particle = (ismember(grains2.id, pairs1)) &...
    (grains2.phase == 1);

% Whether grain boundary segments neighbor a constituent particle
gb_is_at_constituent_particle = all(ismember(gb2_ids,...
    grains2(grains2.at_constituent_particle).id), 2);
gb2.prop.at_constituent_particle = gb_is_at_constituent_particle;

%% Misorientation angles of Al grains around constituent particles
gb_constituent = grains2(grains2.at_constituent_particle).boundary('al', 'al');
gb_other = grains2(~grains2.at_constituent_particle).boundary('al', 'al');
ids_between = grains2(grains2.at_constituent_particle).id;
mask2d = ismember(gb_constituent.grainId, ids_between);
mask1d = all(mask2d, 2);
gb_between = gb_constituent(mask1d);
gb_outwards = gb_constituent(~mask1d);

%% Misorientation angles of grains at constituent particles
fid = fopen(fullfile(dir_mtex, 'gb_at_constituent_mori_angles.txt'), 'w+');
fprintf(fid, '%.5f\n', gb_constituent.misorientation.angle');
fclose(fid);

% Misorientation angles between grains at constituent particles
fid = fopen(fullfile(dir_mtex, 'gb_at_constituent_between_mori_angles.txt'), 'w+');
fprintf(fid, '%.5f\n', gb_between.misorientation.angle');
fclose(fid);

% Misorientation angles outward of grains at constituent particles
fid = fopen(fullfile(dir_mtex, 'gb_at_constituent_outward_mori_angles.txt'), 'w+');
fprintf(fid, '%.5f\n', gb_outwards.misorientation.angle');
fclose(fid);

% Misorientation angles of grains NOT at constituent particles
fid = fopen(fullfile(dir_mtex, 'gb_not_at_constituent_mori_angles.txt'), 'w+');
fprintf(fid, '%.5f\n', gb_other.misorientation.angle');
fclose(fid);

%% Write all relevant data to file
fid = fopen(fullfile(dir_mtex, 'grains.txt'), 'w+');
fprintf(fid, ['#id,phase,size,particle_ecd,ideal,gos,gam,'...
    'at_particle,at_constituent_particle,dist_to_gb,xhab,is_rx\n']);
dataMat = [...
    grains2.id,...
    grains2.phase,...
    grains2.grainSize,...
    grains2.particle_ecd,...
    grains2.ideal_ori_id,...
    grains2.GOS,...
    ebsd2.grainMean(ebsd2.KAM),...
    grains2.at_particle,...
    grains2.at_constituent_particle,...
    grains2.min_distance_to_gb,...
    grains2.xhab,...
    grains2.is_rx,...
];
fprintf(fid, '%i,%i,%i,%.5f,%i,%.5f,%.5f,%i,%i,%.5f,%.5f,%i\n', dataMat');
fclose(fid);

% Grain neighbors
fid = fopen(fullfile(dir_mtex, 'grain_neighbors.txt'), 'w+');
fprintf(fid, '%i,%i\n', grains2.neighbors');
fclose(fid);

% Grain boundaries
al_ids = grains2('al').id;
gb2.prop.is_al = ismember(gb2.grainId(:, 1), al_ids) .*...
    ismember(gb2.grainId(:, 2), al_ids);
gb3 = gb2(gb2.is_al == 1);
fid = fopen(fullfile(dir_mtex, 'grain_boundaries.txt'), 'w+');
fprintf(fid, ['id1,id2,angle,a,b,c,d,is_csl3,is_csl7,is_rx,'...
    'n_dispersoids_close,dispersoids_close_size,'...
    'at_constituent_particle,component1,component2,length\n']);
dataMat = [...
    gb3.grainId(:, 1),...
    gb3.grainId(:, 2),...
    gb3.misorientation.angle,...
    gb3.misorientation.a,...
    gb3.misorientation.b,...
    gb3.misorientation.c,...
    gb3.misorientation.d,...
    gb3.is_csl3,...
    gb3.is_csl7,...
    gb3.is_rx,...
    gb3.n_particles_close,...
    gb3.particles_close_size,...
    gb3.at_constituent_particle,...
    gb3.component1',...
    gb3.component2',...
    gb3.segLength,...
];
fprintf(fid, '%i,%i,%.10f,%.5f,%.5f,%.5f,%.5f,%i,%i,%i,%i,%.5f,%i,%i,%i,%10f\n', dataMat');
fclose(fid);

close all
