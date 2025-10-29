%% calculate the band structure and DOS of 2D PhC using standard plane-wave expansion

close all
clear
clc

%% parameters
epsa = 2;
epsb = 1;
a = 1.0019;
R = 0.148;
num_band = 10;
i = sqrt(-1);
f = pi * (R ^ 2) / (a ^ 2);
a1 = a;
a2 = a * i;
b1 = 2 * pi / a;
b2 = (2 * pi / a) * i;
% n = input('please input n: ');
n = 3;

% number of plane waves for expansion
NumberofPW = (2 * n + 1) ^ 2;

%% diagonalization of Hamiltonian
% generate the G vectors
G = [];

for x = -n:n

    for y = -n:n
        G = [G, x * b1 + y * b2];
    end

end

% calculate the FFT of epsilon E(G-G') analytically for squared lattice
for x = 1:NumberofPW

    for y = x + 1:NumberofPW
        eps2(x, y) = (epsa - epsb) * 2 * f * besselj(1, abs(G(x) - G(y)) * R) ./ (abs(G(x) - G(y)) * R);
        eps2(y, x) = eps2(x, y);
    end

    eps2(x, x) = f * epsa + (1 - f) * epsb;
end

% generate the k points for band calculation
k1 = 2 * pi / a * 0.5 .* (0:0.01:1);
k2 = 2 * pi / a * (0.5 + (0.01:0.01:1) .* 0.5 * i);
k3 = 2 * pi / a * (0.5 + 0.5 * i) .* (0.99:-0.01:0.01);
k0 = [k1, k2, k3];

% generate the M matrix
eps2 = inv(eps2);

for ii = 1:length(k0)
    k = k0(ii);
    M = abs(k + G.') * abs(k + G) .* (eps2); % for TM wave
    % M = (real(k + G.') * real(k + G) + imag(k + G.') * imag(k + G)) .* (eps2); % for TE wave
    E = sort(abs(eig(M)));
    freq(:, ii) = sqrt(abs(E(1:num_band))) .* a ./ 2 ./ pi;
    fprintf('calculation of k=%f+%fi is finished\n', real(k), imag(k));
end

%% plot the band strcuture
freq_max_plot = 1;
tmpx = 1:length(k0);
fig_2D_band = figure;
plot(tmpx, freq, '.', 'linewidth', 1)
xlim([min(tmpx), max(tmpx)])
ylim([0, freq_max_plot])
xticks([1, 100, 200, 300]);
xticklabels({'\Gamma', 'X', 'M', '\Gamma'});
title('Band structure of 2D square lattice')
xlabel('wave vector')
ylabel('\omegaa/2\pi')
grid on
res = 600;
filename = 'fig_2D_band';
print(fig_2D_band, '-dpng', ['-r' num2str(res)], filename);

%% calculate the density of state
step = 0.015;
freq_min = min(min(freq));
freq_max = max(max(freq));
range_freq = freq_max - freq_min;

ndos = zeros(round(range_freq / step + 0.5) + 10, 1); % store the density of states
num_freq_points = round(range_freq / step + 0.5) + 10;
freq_points = linspace(freq_min, freq_max, num_freq_points);

for band = 1:num_band
    index = round((abs(freq(band, :)) - freq_min) / step);

    for ii = 1:length(index)

        if index(ii) == 0
            ndos(index(ii) + 1) = ndos(index(ii) + 1) + 1;
        else
            ndos(index(ii)) = ndos(index(ii)) + 1;
        end

    end

end

fig_2D_dos = figure;
plot(freq_points, ndos, '.-', 'linewidth', 1);
xlim([0, freq_max_plot]);
ylim([0, max(ndos(freq_points <= freq_max_plot))]);
xlabel('\omegaa/2\pi')
ylabel('Density of states');
title('Density of states of 2D square lattice')
grid on
res = 600;
filename = 'fig_2D_dos';
print(fig_2D_dos, '-dpng', ['-r' num2str(res)], filename);
