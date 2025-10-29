function S_fft = S_fft(G, f, R)
    %S_fft - Description
    %
    % Syntax: S_fft = S_fft(G, f, R)
    %
    % generate FFT of S analytically, which is used in moire photonic crystal calculations
    % see >phc_2D_moire_mynote.m< for more details

    % Author: Wenfei Guo, Huazhong University of Science and Technology
    % Date: 06/24/2024

    if abs(G) < 1e-5
        S_fft = f;
    else
        S_fft = 2 * f * besselj(1, abs(G) * R) ./ (abs(G) * R);
    end

end
