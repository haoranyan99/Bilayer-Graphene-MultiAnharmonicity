function lat = build_lattice_potential_simple(par)
% build_lattice_potential_simple
% Read lattice params from par.lat.
%
% Base:
%   V0(psi)= A1*(1 - cos(w*psi)) + 0.5*m*psi^2 + 0.25*k4*psi^4 + s*sqrt(psi^2+eps^2)
%
% NEW: soft cutoff with high-order penalty for |psi|>psi_cut:
%   V = V0 + K*((max(0,|psi|-psi_cut)/Delta)^p)
%   dV = dV0 + dV_pen
%
% par.lat fields:
%   A1, m, k4, s, eps, w
%   psi_cut (default 35)
%   pen_p   (default 16)   % high order
%   pen_K   (default 1e4)  % strength
%   pen_D   (default 1.0)  % wall width scale (Delta)

    if ~isstruct(par)
        error('build_lattice_potential_simple: input must be a struct (par).');
    end
    if ~isfield(par,'lat') || ~isstruct(par.lat)
        error('build_lattice_potential_simple: par.lat must exist and be a struct.');
    end

    P = par.lat;

    % ---------- defaults ----------
    if ~isfield(P,'A1'),      P.A1      = 2.0;    end
    if ~isfield(P,'m'),       P.m       = 0.1;    end
    if ~isfield(P,'k4'),      P.k4      = 0.002;  end
    if ~isfield(P,'s'),       P.s       = 0.20;   end
    if ~isfield(P,'eps'),     P.eps     = 0.02;   end
    if ~isfield(P,'w'),       P.w       = 2.0;    end

    if ~isfield(P,'psi_cut'), P.psi_cut = 30;     end
    if ~isfield(P,'pen_p'),   P.pen_p   = 16;     end
    if ~isfield(P,'pen_K'),   P.pen_K   = 1e4;    end
    if ~isfield(P,'pen_D'),   P.pen_D   = 1.0;    end

    % ---------- checks ----------
    if ~(isfinite(P.eps) && P.eps > 0)
        error('build_lattice_potential_simple: eps must be finite and > 0.');
    end
    if ~(isfinite(P.psi_cut) && P.psi_cut > 0)
        error('build_lattice_potential_simple: psi_cut must be finite and > 0.');
    end
    if ~(isfinite(P.pen_p) && P.pen_p >= 2 && mod(P.pen_p,1)==0)
        error('build_lattice_potential_simple: pen_p must be an integer >= 2.');
    end
    if ~(isfinite(P.pen_K) && P.pen_K > 0)
        error('build_lattice_potential_simple: pen_K must be finite and > 0.');
    end
    if ~(isfinite(P.pen_D) && P.pen_D > 0)
        error('build_lattice_potential_simple: pen_D must be finite and > 0.');
    end

    lat = struct();
    lat.params = P;

    lat.softabs = @(psi) sqrt(psi.^2 + P.eps^2);

    % ---- base (no cutoff) ----
    V0  = @(psi) ...
        P.A1*(1 - cos(P.w*psi)) ...
      + 0.5*P.m*(psi.^2) ...
      + 0.25*P.k4*(psi.^4) ...
      + P.s*sqrt(psi.^2 + P.eps^2);

    dV0 = @(psi) ...
        P.A1*P.w*sin(P.w*psi) ...
      + P.m*psi ...
      + P.k4*(psi.^3) ...
      + P.s*(psi ./ sqrt(psi.^2 + P.eps^2));

    % ---- penalty ----
    Vpen  = @(psi) penalty_V_(psi, P.psi_cut, P.pen_D, P.pen_p, P.pen_K);
    dVpen = @(psi) penalty_dV_(psi, P.psi_cut, P.pen_D, P.pen_p, P.pen_K);

    % ---- total ----
    lat.V  = @(psi) V0(psi)  + Vpen(psi);
    lat.dV = @(psi) dV0(psi) + dVpen(psi);
end

% ===================== penalty helpers =====================

function Vp = penalty_V_(psi, psi_cut, D, p, K)
    % t = max(0, |psi|-psi_cut)
    t = max(0, abs(psi) - psi_cut);
    Vp = K .* (t./D).^p;
end

function dVp = penalty_dV_(psi, psi_cut, D, p, K)
    t = max(0, abs(psi) - psi_cut);

    dVp = zeros(size(psi), 'like', psi);
    mask = t > 0;
    if any(mask(:))
        % derivative: K * p * (t/D)^(p-1) * (1/D) * sign(psi)
        dVp(mask) = K .* p .* (t(mask)./D).^(p-1) .* (1./D) .* sign(psi(mask));
    end
end