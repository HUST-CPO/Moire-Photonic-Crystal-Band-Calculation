function E_z = E_z(c_q, q, r)
    % Calculate the eigenmode E_z for TM mode in plane-wave expansion
    %
    % c_q: coeffcients at each coupled wave vector
    % q: all coupled wave vectors, denoted by complex numbers
    % r: position vector, denoted by a complex number

    % Author: Wenfei Guo, Huazhong University of Science and Technology
    % Date: 06/24/2024

    E_z = 0;

    for ii = 1:length(c_q)
        E_z = E_z + c_q(ii) * exp(1j * (real(q(ii)) * real(r) + imag(q(ii)) * imag(r)));
    end

end
