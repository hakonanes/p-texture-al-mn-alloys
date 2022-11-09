% Misorientation axis and angle distribution of sub-boundaries w/o
% dispersoids at them, per component
%
% Håkon Wiik Ånes (hakon.w.anes@ntnu.no)
% 2022-11-08

clear variables
close all

%% MTEX configuration
plotx2east
plotzIntoPlane
setMTEXpref('FontSize', 20)

% export_fig configuration
res = '-r200';

% Other parameters
c_range = [0 2.7];
c_range_diff = [-1 1];
comp_names = {'random', 'b', 'c', 's', 'cube', 'cubend', 'p'};

% Crystal and specimen symmetry
cs = crystalSymmetry('m-3m', [4.04 4.04 4.04], 'mineral', 'al');

% Directory and file names
sample = '0s'; % 0s, 175c, 300c, 325c
dir_sample = fullfile('/home/hakon/phd/data/p/prover', sample);
fname = 'grain_boundaries.txt';

comp1 = [];
comp2 = [];
mori = [];
is_rx = [];
gb_length = [];
n_dispersoids_close = [];
for i=1:3
    dir_mtex = fullfile(dir_sample, num2str(i), 'mtex');

    % Data matrix contents:
    % (1, 2): (id1, id2)
    % 3: angle
    % (4, 5, 6, 7): (a, b, c, d)
    % (8, 9): (is_csl3, is_csl7)
    % 10: is_rx
    % 11: n_dispersoids_close
    % 12: dispersoids_close_size
    % 13: at_constituent_particle
    % (14, 15): (component1, component2)
    % 16: length
    A = csvread(fullfile(dir_mtex, fname), 1);

    % Texture component
    comp1 = [comp1; A(:, 14)];
    comp2 = [comp2; A(:, 15)];
    
    % Misorientation
    mori = [mori; orientation([A(:, 4), A(:, 5), A(:, 6), A(:, 7)], cs, cs)];

    % Number of dispersoids per boundary length
    is_rx = [is_rx; A(:, 10)];
    gb_length = [gb_length; A(:, 16)];
    n_dispersoids_close = [n_dispersoids_close; A(:, 11)];
end

% Only sub-boundaries
mori = mori(is_rx == 0);
gb_length = gb_length(is_rx == 0);
n_dispersoids_close = n_dispersoids_close(is_rx == 0);

%% Extract relevant quantities
% Sub-boundaries with particles on them
mask = n_dispersoids_close > 0;
mori_disp = mori(mask);
disp_per_length = n_dispersoids_close(mask) ./ gb_length(mask);
comp1_disp = comp1(mask);
comp2_disp = comp2(mask);

% Low angle sub-boundaries only
% All
mask_lagb = mori.angle < 15 * degree;
mori_lagb = mori(mask_lagb);
comp1_lagb = comp1(mask_lagb);
comp2_lagb = comp2(mask_lagb);
% With dispersoids
mask_disp_lagb = mori_disp.angle < 15 * degree;
mori_disp_lagb = mori_disp(mask_disp_lagb);
disp_per_length_lagb = disp_per_length(mask_disp_lagb);
comp1_disp_lagb = comp1_disp(mask_disp_lagb);
comp2_disp_lagb = comp2_disp(mask_disp_lagb);

% Misorientation axes
% All
v_mori = mori.axis;
v_mori.CS = cs;
v_mori_lagb = mori_lagb.axis;
v_mori_lagb.CS = cs;
% With dispersoids
v_mori_disp = mori_disp.axis;
v_mori_disp.CS = cs;
v_mori_disp_lagb = mori_disp_lagb.axis;
v_mori_disp_lagb.CS = cs;

%% Axis density of all GBs and those with particles on them
for i=1:7
    % All
    comp_mask_i = ismember(comp1, i) | ismember(comp2, i);
    v_mori_i = v_mori(comp_mask_i);
    % With dispersoids
    comp_disp_mask_i = ismember(comp1_disp, i) | ismember(comp2_disp, i);
    v_mori_disp_i = v_mori_disp(comp_disp_mask_i);

    % With dispersoids and all
    figure
    h1 = plot(v_mori_disp_i, 'fundamentalSector', 'contourf');
    mtexTitle('With dispersoids')
    nextAxis
    h2 = plot(v_mori_i, 'fundamentalSector', 'contourf', 'colorrange',...
        c_range);
    mtexTitle('All')
    mtexColorMap inferno
    mtexColorbar('title', 'MRD')
    export_fig(fullfile(dir_sample,...
        ['mori_disp_axis_distributions_' comp_names{i} '.png']), res)

    % Difference
    figure
    plot(v_mori_i, 'fundamentalSector')
    hold on
    v_grid = plotS2Grid(cs.fundamentalSector);
    plot(v_grid, h1.ZData - h2.ZData, 'fundamentalSector', 'contourf',...
        'colorrange', c_range_diff);
    hold off
    mtexColorMap blue2red
    mtexColorbar('title', 'MRD')
    export_fig(fullfile(dir_sample, ...
        ['mori_disp_axis_distributions_diff_' comp_names{i} '.png']), res)
end

%% Plot angle histograms and differences between them for all GB
bin_edges = linspace(0, maxAngle(cs), 20);
bin_midpoints = bin_edges(2:end) - 0.5 * bin_edges(2);
bin_midpoints_deg = bin_midpoints / degree;
bar_width = 0.5 * bin_edges(2) / degree;

for i=1:7
    % All
    comp_mask_i = ismember(comp1, i) | ismember(comp2, i);
    mori_i = mori(comp_mask_i);
    mori_angles_dens_i = histcounts(mori_i.angle, bin_edges,...
        'Normalization', 'probability');
    % With dispersoids
    comp_disp_mask_i = ismember(comp1_disp, i) | ismember(comp2_disp, i);
    mori_disp_i = mori_disp(comp_disp_mask_i);
    mori_disp_angles_dens_i = histcounts(mori_disp_i.angle, bin_edges,...
        'Normalization', 'probability');

    % W/o dispersoids
    figure
    plotAngleDistribution(mori_disp_i)
    hold on
    h = plotAngleDistribution(mori_i);
    h(1).FaceColor = 'r';
    h(2).FaceColor = 'b';
    legend('With dispersoids', 'All')
    ylim([0 25])
    export_fig(fullfile(dir_sample,...
        ['mori_angle_distributions_' comp_names{i} '.png']), res)

    % Difference for all GB
    mori_angle_diff_i = mori_disp_angles_dens_i - mori_angles_dens_i;

    figure
    for j=1:(length(bin_edges) - 1)
        mori_angle_diff_j = mori_angle_diff_i(j);
        if mori_angle_diff_j > 0
            bar_color = 'red';
        else
            bar_color = 'blue';
        end
        bar(bin_midpoints_deg(j), mori_angle_diff_j, bar_width,...
            'FaceColor', bar_color)
        hold on
    end
    xlabel('Misorientation angle [deg]')
    ylabel('Frequency difference (%)')
    ylim([-0.15 0.15])
    set(gcf, 'color', 'w');
    export_fig(fullfile(dir_sample,...
        ['mori_angle_distributions_diff_' comp_names{i} '.png']), res)
end

%% Axis density of all LAGBs and those with particles on them
figure
h1 = plot(v_mori_disp_lagb, 'fundamentalSector', 'contourf');
mtexTitle('With dispersoids')
nextAxis
h2 = plot(v_mori_lagb, 'fundamentalSector', 'contourf', 'colorrange', [0.5 2]);
mtexTitle('All')
mtexColorMap inferno
mtexColorbar('title', 'MRD')
export_fig(fullfile(dir_sample, 'mori_disp_axis_distributions_lagb.png'), res)

figure
plot(v_mori_lagb, 'fundamentalSector')
hold on
v_grid = plotS2Grid(cs.fundamentalSector);
%plot(v_grid, h1.ZData - h2.ZData, 'fundamentalSector', 'contourf', 'colorrange', c_range);
plot(v_grid, h1.ZData - h2.ZData, 'fundamentalSector', 'contourf');
hold off
mtexTitle('With dispersoids - All')
mtexColorMap blue2red
mtexColorbar('title', 'MRD')
export_fig(fullfile(dir_sample, 'mori_disp_axis_distributions_lagb_diff.png'), res)