function el = build_electronic_energy_poly(P)
% build_electronic_energy_poly
% Construct electronic free-energy part using the SAME convention as your uploaded free_energy.m:
%
%   Fel(psi1) = (1/2)*a1*psi1^2 + (1/factorial(4))*b1*psi1^4
%   dFel/dpsi1 = a1*psi1 + (1/6)*b1*psi1^3
%
% Required fields in P:
%   P.a1, P.b1
%
% Output:
%   el.params : parameter struct
%   el.F(psi1): electronic energy
%   el.dF(psi1): derivative

    arguments
        P struct
    end

    % ---- defaults (if missing) ----
    if ~isfield(P,'a1'), P.a1 = 0.0; end
    if ~isfield(P,'b1'), P.b1 = 1.0; end

    a1 = double(P.a1);
    b1 = double(P.b1);

    el = struct();
    el.params = P;

    el.F  = @(psi1) 0.5*a1*(psi1.^2) + (1/factorial(4))*b1*(psi1.^4);

    % d/dpsi1: 0.5*a1*2*psi1 = a1*psi1
    %          (1/24)*b1*4*psi1^3 = (1/6)*b1*psi1^3
    el.dF = @(psi1) a1*psi1 + (1/6)*b1*(psi1.^3);
end