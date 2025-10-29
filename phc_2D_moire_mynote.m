%% program to calculate the band structure of 2D moire photonic crystal
% use coupling and energy cut-off

% Author: Wenfei Guo, Huazhong University of Science and Technology
% Date: 06/24/2024

tic
close all
clear
clc
diary 'log.txt'
diary on

%% parameter definition
theta = 45; % twist angle between layers in degree
theta = (theta / 180) * pi;
a = 1.0119; % lattice constant
R = 0.148; % radius of nanodisks
f = pi * (R ^ 2) / (a ^ 2); % filling fraction
i = sqrt(-1);
a1 = 2 * pi / a; % reciprocal vectors of the first layer
a2 = (2 * pi / a) * i;
b1 = (2 * pi / a) * cos(theta) + 1j * (2 * pi / a) * sin(theta); % reciprocal vectors of the second layer
b2 =- (2 * pi / a) * sin(theta) + 1j * (2 * pi / a) * cos(theta);

nc = 2; % threshold of coupling cut-off
ne = 1; % threshold to determine whether q is at the edge of the coupling cut-off in k space
qc_g1 = 2.1; % threshold of energy cut-off
qc = qc_g1 * abs(a1 + a2);
edge_state_crit = 1 / exp(1); % threshold to determine whether a state is a momentum edge state
num_band = 600; % number of calculated bands

eps_rad = 2; % permittivity of nanodisks
eps_bg = 1; % permittivity of the background material
A = 1 / eps_rad - 1 / eps_bg;
B = 1 / eps_bg;
C = 1 / eps_rad - 1 / eps_bg;

%% generate the k points for band calculation
delta_k = 0.01; % step to sweep the FBZ
k1 = a1 * 0.5 .* (0:delta_k:1); % Γ to X
k2 = a1 * 0.5 + a2 * 0.5 .* (delta_k:delta_k:1); % X to M
k3 = (a1 * 0.5 + a2 * 0.5) .* (1 - delta_k:-delta_k:delta_k); % M to Γ
k0 = [k1, k2, k3];
fig_k = figure; % visualize the wave vectors k
plot(real(k0), imag(k0), '.', 'linewidth', 1, 'color', 'red')
xlabel('k_x')
ylabel('k_y')
title('All wave vectors k in the FBZ of the first lattice for band calculation')
res = 600;
filename = 'fig_k';
print(fig_k, '-dpng', ['-r' num2str(res)], filename);

%% diagonalization of Hamiltonian
bar = waitbar(0, 'Computation started');

for ii = 1:length(k0)
    str = ['Computing...', num2str(100 * roundn(ii / length(k0), -4)), '%'];
    wtbar = waitbar(ii / length(k0), bar, str);

    % For each k, a Hamiltonian is constructed and diagonalized to calculate the eigenfrequency and eigenmode
    k = k0(ii);
    q = [];
    G2 = [];
    G1 = [];

    % the flag to show whether q is at the edge of the coupling cut-off in k space
    is_edge_q = [];

    % generate all the n1 and n2 which satisfy coupling cut-off
    n1_n2 = [];

    for x = -nc:1:nc

        for y = -nc:1:nc
            n1_n2 = [n1_n2, [x, y]'];
        end

    end

    % generate all the coupled wave vectors q for each initial k
    count = 1;

    for jj = 1:length(n1_n2(1, :))
        n1 = n1_n2(1, jj);
        n2 = n1_n2(2, jj);
        G2 = [G2, n1 * b1 + n2 * b2];

        % generate all the m1 and m2 which satisfy energy cut-off
        m1_m2 = [];

        for x = -floor(5 * qc_g1):1:floor(5 * qc_g1)

            for y = -floor(5 * qc_g1):1:floor(5 * qc_g1)
                temp = k + x * a1 + y * a2 + n1 * b1 + n2 * b2;

                if abs(temp) <= qc
                    m1_m2 = [m1_m2, [x, y]'];
                end

            end

        end

        q_each_n_temp = [];

        for kk = 1:length(m1_m2(1, :))
            m1 = m1_m2(1, kk);
            m2 = m1_m2(2, kk);
            G1 = [G1, m1 * a1 + m2 * a2];
            q_each_n_temp = [q_each_n_temp, k + m1 * a1 + m2 * a2 + n1 * b1 + n2 * b2];
            q = [q, k + m1 * a1 + m2 * a2 + n1 * b1 + n2 * b2];

            % Determine whether q is at the edge of the coupling cut-off in k space
            if (abs(n1 - nc) < ne) || (abs(n1 + nc) < ne) || (abs(n2 - nc) < ne) || (abs(n2 + nc) < ne)
                is_edge_q = [is_edge_q, true];
            else
                is_edge_q = [is_edge_q, false];
            end

        end

        q_each_n.q = q_each_n_temp;
        q_each_n.n = [n1, n2];
        q_each_n.k = k;
        q_n_k(count, ii) = q_each_n;
        count = count + 1;

    end

    is_edge_q = logical(is_edge_q);
    H = zeros(length(q), length(q));

    flag3 = [];
    flag4 = [];

    % construct the Hamiltonian H based on the central equation
    for x = 1:length(q)

        for y = x:length(q)

            if x == y
                H(x, y) = (B + A * S_fft(0, f, R) + C * S_fft(0, f, R)) * abs(q(x)) ^ 2;
            else
                k3 = min(abs(G1 - (q(x) - q(y))));

                if k3 < 1e-10
                    H(x, y) = A * S_fft(q(x) - q(y), f, R) * abs(q(x)) * abs(q(y));
                    H(y, x) = H(x, y);
                    flag3 = [flag3, true];
                else

                    k4 = min(abs(G2 - (q(x) - q(y))));

                    if k4 < 1e-10
                        H(x, y) = C * S_fft(q(x) - q(y), f, R) * abs(q(x)) * abs(q(y));
                        H(y, x) = H(x, y);
                        flag4 = [flag4, true];
                    end

                end

            end

        end

    end

    % Visualize the Hamiltonian at k=0 and all coupled wave vectors q
    if k == 0
        H0 = H;
        q0 = q;

        fig_hami = figure;
        imagesc(H0);
        xlabel('column')
        ylabel('row')
        colorbar;
        title('The constructed Hamiltonian at wave vector k=0')
        res = 600;
        filename = 'fig_hami';
        print(fig_hami, '-dpng', ['-r' num2str(res)], filename);

        fig_coupled_q = figure;
        plot(real(q0), imag(q0), 'o', 'linewidth', 1, 'color', 'red')
        xlabel('k_x')
        ylabel('k_y')
        title('All coupled wave vectors q at k=0')
        res = 600;
        filename = 'fig_coupled_q';
        print(fig_coupled_q, '-dpng', ['-r' num2str(res)], filename);
    end

    [V, D] = eig(H); % Diagonalize the Hamiltonian H

    eigen_val = diag(D);
    index_to_delete = find(eigen_val < 0);
    V(:, index_to_delete) = [];
    eigen_val(index_to_delete) = [];

    % Determine whether an eigenstate is a momentum edge state
    is_edge_state = [];

    for jj = 1:length(V(1, :))
        state = V(:, jj);

        if sum((abs(state(is_edge_q))) .^ 2) >= edge_state_crit
            is_edge_state = [is_edge_state, true];
        else
            is_edge_state = [is_edge_state, false];
        end

    end

    % Take >num_band< bands and eigenmodes
    V = V(:, 1:num_band);
    eigen_val = eigen_val(1:num_band);
    is_edge_state = is_edge_state(1:num_band);
    is_edge_state = logical(is_edge_state);
    not_edge_state = ~is_edge_state;

    freq_k.edge = sqrt(eigen_val(is_edge_state)) .* a ./ (2 * pi);
    freq_k.real = sqrt(eigen_val(not_edge_state)) .* a ./ (2 * pi);
    freq_k.k = k;
    freq(ii) = freq_k;
    eigenmode_k.edge = V(:, is_edge_state);
    eigenmode_k.real = V(:, not_edge_state);
    eigenmode_k.k = k;
    eigenmode(ii) = eigenmode_k;
    q_current_k.q = q;
    q_current_k.k = k;
    q_k(ii) = q_current_k;
    fprintf('calculation of k=%f+%fi is finished\n', real(k), imag(k));

end

close(wtbar);

%% plot the band structure
tmpx = 1:length(k0);
fig_2D_Moire_band = figure;
plot(tmpx(1), freq(1).real, '.', 'linewidth', 1, 'color', 'red')
hold on

for ii = 2:length(k0)
    plot(tmpx(ii), freq(ii).real, '.', 'linewidth', 1, 'color', 'red')
end

% plot the edge state band
% plot(tmpx(1), freq(1).edge, '.', 'linewidth', 1, 'color', 'black')
%
% for ii = 2:length(k0)
%     plot(tmpx(ii), freq(ii).edge, '.', 'linewidth', 1, 'color', 'black')
% end

hold off
title('Band structure of 2D Moire photonic crystal')
xlabel('wave vector')
ylabel('freq')
grid on

% save the results
clear('fig_k');
clear('fig_hami');
clear('fig_coupled_q');
clear('fig_2D_Moire_band');
save('phc_2D_moire_data.mat');
toc
diary off
