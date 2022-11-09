% Misorientation axis and angle distribution of sub-boundaries w/o
% dispersoids at them
%
% Håkon Wiik Ånes (hakon.w.anes@ntnu.no)
% 2022-11-08

clear variables
close all

%% Read data
% MTEX configuration
plotx2east
plotzIntoPlane
setMTEXpref('FontSize', 40)

% export_fig configuration
res = '-r200';

% Other parameters
c_range = [-0.1 0.1];

% Crystal and specimen symmetry
cs = crystalSymmetry('m-3m', [4.04 4.04 4.04], 'mineral', 'al');

% Directory and file names
sample = '300c'; % 0s, 175c, 300c, 325c
dir_sample = fullfile('/home/hakon/phd/data/p/prover', sample);
fname = 'grain_boundaries.txt';

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

    % Misorientation
    mori = [mori; orientation([A(:, 4), A(:, 5), A(:, 6), A(:, 7)], cs, cs)];

    % Number of dispersoids per boundary length
    is_rx = [is_rx; A(:, 10)];
    gb_length = [gb_length; A(:, 16)];
    n_dispersoids_close = [n_dispersoids_close; A(:, 11)];
end

% Recrystallized boundaries only
mask_rx = is_rx == 1;
if any(mask_rx)
    mori_rx = mori(mask_rx);
    gb_length_rx = gb_length(mask_rx);
    n_dispersoids_close_rx = n_dispersoids_close(mask_rx);
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

% Low angle boundaries only
mori_lagb = mori(mori.angle < 15 * degree);
mori_disp_lagb = mori_disp(mori_disp.angle < 15 * degree);
disp_per_length_lagb = disp_per_length(mori_disp.angle < 15 * degree);
gb_length_lagb = gb_length(mori.angle < 15 * degree);

% Misorientation axes
v_mori = mori.axis;
v_mori_disp = mori_disp.axis;
v_mori_lagb = mori_lagb.axis;
v_mori_disp_lagb = mori_disp_lagb.axis;
v_mori.CS = cs;
v_mori_disp.CS = cs;
v_mori_lagb.CS = cs;
v_mori_disp_lagb.CS = cs;

%% Fraction of sigma7 sub-boundaries
frac_sigma7 = sum(angle(mori, CSL(7, cs)) < 15 * degree) / size(mori, 1);
frac_sigma7_disp = sum(angle(mori_disp, CSL(7, cs)) < 15 * degree) / size(mori_disp, 1);
disp(frac_sigma7_disp / frac_sigma7)

%% Axis density of all GBs and those with particles on them
mori_disp_dens = calcDensity(v_mori_disp, 'weights', disp_per_length);
mori_dens = calcDensity(v_mori, 'weights', gb_length);

figure
plot(mori_disp_dens, 'fundamentalSector', 'contourf');
mtexTitle('With dispersoids')
nextAxis
plot(mori_dens, 'fundamentalSector', 'contourf', 'colorrange', [0.5 2]);
mtexTitle('All')
mtexColorMap inferno
mtexColorbar('title', 'MRD')
export_fig(fullfile(dir_sample, 'mori_disp_axis_distributions.png'), res)

pause(0.5)

figure
plot(mori_disp_dens - mori_dens, 'fundamentalSector', 'contourf',...
    'colorrange', c_range)
mtexColorMap blue2red
if strcmp(sample, '325c')
    mtexColorbar('title', 'MRD difference')
end
export_fig(fullfile(dir_sample, 'mori_disp_axis_distributions_diff.png'), res)

%% Plot angle histograms and differences between them for all GB
figure
plotAngleDistribution(mori_disp)
hold on
h = plotAngleDistribution(mori);
h(1).FaceColor = 'r';
h(2).FaceColor = 'b';
legend('With dispersoids', 'All')
ylim([0 25])
export_fig(fullfile(dir_sample, 'mori_angle_distributions.png'), res)

%% Difference for all GB
max_omega = maxAngle(cs);
bins = linspace(-eps, max_omega + 0.01, 15);
bin_midpoints = 0.5*(bins(1:end - 1) + bins(2:end));

mori_angles_dens = histcounts(mori.angle, bins, 'Normalization',...
    'probability');
mori_disp_angles_dens = histcounts(mori_disp.angle, bins,...
    'Normalization', 'probability');
mori_angle_diff = 100 * (mori_disp_angles_dens - mori_angles_dens);

figure
for i=1:length(bin_midpoints)
    mori_angle_diff_i = mori_angle_diff(i);
    if mori_angle_diff_i > 0
        bar_color = 'red';
    else
        bar_color = 'blue';
    end
    bar(bin_midpoints(i) / degree, mori_angle_diff_i, 'FaceColor', bar_color)
    hold on
end
xlabel('Misorientation angle (degrees)')
ylabel('Frequency difference (%)')
ylim([-1.2 1.2])
set(gcf, 'color', 'w');
export_fig(fullfile(dir_sample, 'mori_angle_distributions_diff.png'), res)

%% Axis density of all LAGBs and those with particles on them
mori_disp_lagb_dens = calcDensity(v_mori_disp_lagb, 'weights',...
    disp_per_length_lagb);
mori_lagb_dens = calcDensity(v_mori_lagb, 'weights', gb_length_lagb);

figure
plot(mori_disp_lagb_dens, 'fundamentalSector', 'contourf');
mtexTitle('With dispersoids')
nextAxis
plot(mori_lagb_dens, 'fundamentalSector', 'contourf', 'colorrange', [0.5 2]);
mtexTitle('All')
mtexColorMap inferno
mtexColorbar('title', 'MRD')
export_fig(fullfile(dir_sample, 'mori_lagb_disp_axis_distributions.png'), res)

pause(0.5)

figure
plot(mori_disp_lagb_dens - mori_lagb_dens, 'fundamentalSector',...
    'contourf', 'colorrange', c_range)
mtexColorMap blue2red
if strcmp(sample, '325c')
    mtexColorbar('title', 'MRD difference')
end
export_fig(fullfile(dir_sample, 'mori_lagb_disp_axis_distributions_diff.png'), res)

%% Recrystallized boundaries
if any(mask_rx)
    % Get misorientations of boundaries with dispersoids
    mask_rx_disp = n_dispersoids_close_rx > 0;
    mori_rx_disp = mori_rx(mask_rx_disp);
    disp_per_length_rx = n_dispersoids_close_rx(mask_rx_disp) ./ gb_length_rx(mask_rx_disp);

    % Get misorientation axes w/o dispersoids
    v_mori_rx = mori_rx.axis;
    v_mori_rx_disp = mori_rx_disp.axis;
    v_mori_rx.CS = cs;
    v_mori_rx_disp.CS = cs;

    % Fraction of sigma7 boundaries
    frac_sigma7_rx = sum(angle(mori_rx, CSL(7, cs)) < 15 * degree) / size(mori_rx, 1);
    frac_sigma7_rx_disp = sum(angle(mori_rx_disp, CSL(7, cs)) < 15 * degree) / size(mori_rx_disp, 1);
    disp(frac_sigma7_rx_disp / frac_sigma7_rx)

    % Axis density of all GBs and those with dispersoids on them
    mori_rx_disp_dens = calcDensity(v_mori_rx_disp, 'weights',...
    disp_per_length_rx);
    mori_rx_dens = calcDensity(v_mori_rx, 'weights', gb_length_rx);

    figure
    plot(mori_rx_disp_dens, 'fundamentalSector', 'contourf');
    mtexTitle('With dispersoids')
    nextAxis
    plot(mori_rx_dens, 'fundamentalSector', 'contourf', 'colorrange', [0 2]);
    mtexTitle('All')
    mtexColorMap inferno
    mtexColorbar('title', 'MRD')
    export_fig(fullfile(dir_sample, 'mori_rx_disp_axis_distributions.png'), res)

    figure
    plot(mori_rx_disp_dens - mori_rx_dens, 'fundamentalSector',...
        'contourf', 'colorrange', [-0.5 0.5])
    mtexColorMap blue2red
    if strcmp(sample, '325c')
        mtexColorbar('title', 'MRD difference')
    end
    export_fig(fullfile(dir_sample, 'mori_rx_disp_axis_distributions_diff.png'), res)

    % Plot angle histograms and differences between them for all GB
    figure
    plotAngleDistribution(mori_rx_disp)
    hold on
    h = plotAngleDistribution(mori_rx);
    h(1).FaceColor = 'r';
    h(2).FaceColor = 'b';
    legend('With dispersoids', 'All')
    ylim([0 25])
    export_fig(fullfile(dir_sample, 'mori_rx_angle_distributions.png'), res)

    % Difference for all GB
    max_omega = maxAngle(cs);
    bins = linspace(-eps, max_omega + 0.01, 15);
    bin_midpoints = 0.5*(bins(1:end - 1) + bins(2:end));
    mori_rx_disp_dens = histcounts(mori_rx_disp.angle, bins,...
        'Normalization', 'probability');
    mori_rx_dens = histcounts(mori_rx.angle, bins, 'Normalization',...
        'probability');
    mori_rx_angle_diff = 100 * (mori_rx_disp_dens - mori_rx_dens);

    figure
    for i=1:length(bin_midpoints)
        mori_rx_angle_diff_i = mori_rx_angle_diff(i);
        if mori_rx_angle_diff_i > 0
            bar_color = 'red';
        else
            bar_color = 'blue';
        end
        bar(bin_midpoints(i) / degree, mori_rx_angle_diff_i, 'FaceColor', bar_color)
        hold on
    end
    xlabel('Misorientation angle (degrees)')
    ylabel('Frequency difference (%)')
    ylim([-6 6])
    set(gcf, 'color', 'w');
    export_fig(fullfile(dir_sample, 'mori_rx_angle_distributions_diff.png'), res)
end

%%
close all