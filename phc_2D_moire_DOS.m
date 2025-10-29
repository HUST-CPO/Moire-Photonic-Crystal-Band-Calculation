%% density of states of 2D moire photonic crystals based on calculated bands
% visualize the band structure
% calculate and visualize the density of states out of calculated band
% calculate and visualize the eigenmode out of eigenvectors

% Author: Wenfei Guo, Huazhong University of Science and Technology
% Date: 06/24/2024

close all
clear
clc

load('phc_2D_moire_data.mat');

%% plot the band structure
freq_max_plot = 0.6;
tmpx = 1:length(k0);
fig_2D_Moire_band = figure;
freq_temp = freq(1).real;
freq_temp = freq_temp(freq_temp < freq_max_plot);

if ~isempty(freq_temp)
    plot(tmpx(1), freq_temp, '.', 'linewidth', 1, 'color', 'red')
end

xlim([min(tmpx), max(tmpx)])
xticks([1, 100, 200, 300]);
xticklabels({'\Gamma', 'X', 'M', '\Gamma'});
hold on

for ii = 2:length(k0)
    freq_temp = freq(ii).real;
    freq_temp = freq_temp(freq_temp < freq_max_plot);

    if ~isempty(freq_temp)
        plot(tmpx(ii), freq_temp, '.', 'linewidth', 1, 'color', 'red')
    end

end

freq_temp = freq(1).edge;
freq_temp = freq_temp(freq_temp < freq_max_plot);

if ~isempty(freq_temp)
    plot(tmpx(1), freq_temp, '.', 'linewidth', 1, 'color', 'black')
end

for ii = 2:length(k0)
    freq_temp = freq(ii).edge;
    freq_temp = freq_temp(freq_temp < freq_max_plot);

    if ~isempty(freq_temp)
        plot(tmpx(ii), freq_temp, '.', 'linewidth', 1, 'color', 'black')
    end

end

hold off
title('Band structure of 2D Moire photonic crystal with edge states')
xlabel('wave vector')
ylabel('\omegaa/2\pi')
grid on
res = 600;
filename = 'fig_2D_Moire_band';
print(fig_2D_Moire_band, '-dpng', ['-r' num2str(res)], filename);

% plot the band structure without edge states
tmpx = 1:length(k0);
fig_2D_Moire_band_remove_edge = figure;
freq_temp = freq(1).real;
freq_temp = freq_temp(freq_temp < freq_max_plot);

if ~isempty(freq_temp)
    plot(tmpx(1), freq_temp, '.', 'linewidth', 1, 'color', 'red')
end

xlim([min(tmpx), max(tmpx)])
xticks([1, 100, 200, 300]);
xticklabels({'\Gamma', 'X', 'M', '\Gamma'});
hold on

for ii = 2:length(k0)
    freq_temp = freq(ii).real;
    freq_temp = freq_temp(freq_temp < freq_max_plot);

    if ~isempty(freq_temp)
        plot(tmpx(ii), freq_temp, '.', 'linewidth', 1, 'color', 'red')
    end

end

hold off
title('Band structure of 2D Moire photonic crystal')
xlabel('wave vector')
ylabel('\omegaa/2\pi')
grid on
res = 600;
filename = 'fig_2D_Moire_band_remove_edge';
print(fig_2D_Moire_band_remove_edge, '-dpng', ['-r' num2str(res)], filename);

%% calculate and visualize the density of state
step = 0.005; % step for discretizing the energy band
freq_max = max(abs(freq(1).real));
freq_min = min(abs(freq(1).real));

% Take the highest and lowest frequency of the band and calculate the range of the band
for ii = 1:length(k0)

    if max(abs(freq(ii).real)) > freq_max
        freq_max = max(abs(freq(ii).real));
    end

    if min(abs(freq(ii).real)) < freq_min
        freq_min = min(abs(freq(ii).real));
    end

    range_freq = freq_max - freq_min;
end

ndos = zeros(round(range_freq / step + 0.5) + 10, 1); % Store density of states
num_freq_points = round(range_freq / step + 0.5) + 10;
freq_points = linspace(freq_min, freq_max, num_freq_points);

% Count the number of eigenstates in each interval
for ii = 1:length(k0)
    index = round((abs(freq(ii).real) - freq_min) / step + 0.5);

    for jj = 1:length(index)

        ndos(index(jj)) = ndos(index(jj)) + 1;

    end

end

fig_2D_Moire_dos = figure; % Visualize the density of states
indices = freq_points <= freq_max_plot;
plot(freq_points(indices), ndos(indices), '.-', 'linewidth', 1);
% ylim([0, max(ndos(indices))])
xlabel('\omegaa/2\pi');
ylabel('Density of states');
title('Density of states of 2D Moire photonic crystal')
grid on
res = 600;
filename = 'fig_2D_Moire_dos';
print(fig_2D_Moire_dos, '-dpng', ['-r' num2str(res)], filename);

%% calculate the eigenmode with lowest frequency at k=0
num_points = 200;
x = linspace(-10 * a, 10 * a, num_points);
y = linspace(-10 * a, 10 * a, num_points);
[X, Y] = meshgrid(x, y); % Discretize the space
r = X + Y .* 1j;

% you can change the index in >eigenmode< to visualize any eigenstate
E_z_values = E_z(eigenmode(1).real(:, 1), q_k(1).q, r);

fig_2D_Moire_eigenmode_3D = figure;
surf(X, Y, abs(E_z_values), 'EdgeColor', 'none');
xlabel('X');
ylabel('Y');
zlabel('|E_z(r)|');
title('Eigenmode E_z with lowest frequency at k=0');
colorbar;

fig_2D_Moire_eigenmode = figure;
contourf(X, Y, abs(E_z_values), 50, 'EdgeColor', 'none');
xlabel('X');
ylabel('Y');
title('Eigenmode E_z with lowest frequency at k=0');
cb2 = colorbar;
set(get(cb2, 'title'), 'string', '|E_z|');
res = 600;
filename = 'fig_2D_Moire_eigenmode';
print(fig_2D_Moire_eigenmode, '-dpng', ['-r' num2str(res)], filename);
