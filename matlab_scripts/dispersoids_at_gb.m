% Dispersoids at (sub)grain boundaries
%
% Håkon Wiik Ånes (hakon.w.anes@ntnu.no)
% 2022-11-01

clear variables
close all

%% MTEX configuration
plotx2east
plotzIntoPlane

% export_fig configuration
res = '-r200';

% Crystal and specimen symmetry
cs = crystalSymmetry('m-3m', [4.04 4.04 4.04], 'mineral', 'al');

% Directory and file names
sample = '325c'; % 0s, 175c, 300c, 325c
dir_sample = fullfile('/home/hakon/phd/data/p/prover', sample);
fname = 'grain_boundaries.txt';

mori = [];
gb_length = [];
n_dispersoids_close = [];
for i=1:3
    dir_mtex = fullfile(dir_sample, num2str(i), 'mtex');

    % Data matrix contents:
    % (1, 2): (id1, id2)
    % 3: angle
    % (4, 5, 6, 7): (a, b, c, d)
    % (8, 9): (is_csl3, is_csl7)
    % 10: n_dispersoids_close
    % 11: dispersoids_close_size
    % 12: at_constituent_particle
    % (13, 14): (component1, component2)
    % 15: length
    A = csvread(fullfile(dir_mtex, fname), 1);

    % Misorientation
    mori = [mori; orientation([A(:, 4), A(:, 5), A(:, 6), A(:, 7)], cs, cs)];

    % Number of dispersoids per boundary length
    gb_length = [gb_length; A(:, 15)];
    n_dispersoids_close = [n_dispersoids_close; A(:, 10)];
end

%% Extract relevant quantities
mask = n_dispersoids_close > 0;

% Keep only data of boundaries with particles on them
mori_disp = mori(mask);
disp_per_length = n_dispersoids_close(mask) ./ gb_length(mask);

% Misorientation axes
v_mori = mori.axis;
v_mori_disp = mori_disp.axis;
v_mori.CS = cs;
v_mori_disp.CS = cs;

%% Plot axis density of all boundaries and those with particles on them
figure
h1 = plot(v_mori_disp, 'fundamentalSector', 'contourf');
mtexTitle('With dispersoids')
nextAxis
h2 = plot(v_mori, 'fundamentalSector', 'contourf', 'colorrange', [0.5 2]);
mtexTitle('All')
mtexColorMap inferno
mtexColorbar('title', 'MRD')
export_fig(fullfile(dir_sample, 'mori_disp_axis_distributions.png'), res)

figure
plot(v_mori, 'fundamentalSector')
hold on
v_grid = plotS2Grid(cs.fundamentalSector);
plot(v_grid, h1.ZData - h2.ZData, 'fundamentalSector', 'contourf', 'colorrange', [-0.1 0.1]);
hold off
mtexTitle('With dispersoids - All')
mtexColorMap blue2red
mtexColorbar('title', 'MRD')
export_fig(fullfile(dir_sample, 'mori_disp_axis_distributions_diff.png'), res)

%% Plot angle histograms and differences between them
bin_edges = linspace(0, maxAngle(cs), 20);
bin_midpoints = bin_edges(2:end) - 0.5 * bin_edges(2);
bar_width = 0.5 * bin_edges(2) / degree;

mori_angles_dens = histcounts(mori.angle, bin_edges, 'Normalization', 'probability');
mori_disp_angles_dens = histcounts(mori_disp.angle, bin_edges, 'Normalization', 'probability');

figure
plotAngleDistribution(mori_disp)
hold on
h = plotAngleDistribution(mori);
h(1).FaceColor = 'r';
h(2).FaceColor = 'b';
legend('With dispersoids', 'All')
ylim([0 25])
export_fig(fullfile(dir_sample, 'mori_angle_distributions.png'), res)

%% Differences
mori_angle_diff = mori_disp_angles_dens - mori_angles_dens;
bin_midpoints_deg = bin_midpoints / degree;

figure
for i=1:(length(bin_edges) - 1)
    mori_angle_diff_i = mori_angle_diff(i);
    if mori_angle_diff_i > 0
        bar_color = 'red';
    else
        bar_color = 'blue';
    end
    bar(bin_midpoints_deg(i), mori_angle_diff_i, bar_width, 'FaceColor', bar_color)
    hold on
end
xlabel('Misorientation angle [deg]')
ylabel('Frequency difference (%)')
ylim([-0.01 0.01])

set(gcf, 'color', 'w');
export_fig(fullfile(dir_sample, 'mori_angle_distributions_diff.png'), res)