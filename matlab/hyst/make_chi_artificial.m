function chi = make_chi_artificial(par, Pchi)
% make_chi_artificial
% Artificial chi(T,dop):
%   low --(sigmoid up #1)-> plateau1 --(sigmoid up #2)-> plateau2
%       --(right decay onset)-> exponential decay (sigmoid-gated)
%
% Usage:
%   chi = make_chi_artificial(par);
%   chi = make_chi_artificial(par, Pchi);
%   y   = chi(T, dop);

    if nargin < 2 || isempty(Pchi)
        Pchi = struct();
    end

    % -------- defaults --------
    def.low       = 0.0;
    def.peak      = 1.0;

    def.dop_step1 = 1.72;    % first climb center
    def.w1_0      = 1e-3;    % first climb width
    def.h1        = 0.80;     % first plateau added height

    def.dop_step2 = 1.78;    % second climb center
    def.w2_0      = 2e-3;    % second climb width
    def.h2        = def.peak - def.h1 - def.low;     % second plateau added height (total top = low+h1+h2)

    def.dop_peak  = 1.86;    % optional Gaussian bump center
    def.sig0      = 5e-3;
    def.peakAmp   = 0;

    def.dop_decay = 1.9;    % decay onset
    def.wd_0      = 2e-3;    % decay "turn-on" width (sigmoid)
    def.decayLen0 = 0.10;    % exponential decay length scale
    def.decayAmp  = 1.0;     % 1 -> fully decay; 0.5 -> decay half; 0 disables
    % --------------------------

    Pchi = fill_defaults_(Pchi, def);

    if ~isfield(par,'el') || ~isfield(par.el,'T_eps')
        epsT = 1e-6;
    else
        epsT = par.el.T_eps;
    end

    chi = @(T,dop) chi_eval_(T, dop, Pchi, epsT);
end

% ================= INTERNAL =================

function y = chi_eval_(T, dop, P, epsT)

    Tuse = max(abs(T), epsT);

    % mild T broadening (keep small; tune factors if needed)
    w1      = P.w1_0  * (1 + 0.2*Tuse);
    w2      = P.w2_0  * (1 + 0.2*Tuse);
    wd      = P.wd_0  * (1 + 0.2*Tuse);
    sig     = P.sig0  * (1 + 0.5*Tuse);
    decayLn = P.decayLen0 * (1 + 0.1*Tuse);

    % two-step climb
    s1 = 1 ./ (1 + exp(-(dop - P.dop_step1)./w1));  % 0->1
    s2 = 1 ./ (1 + exp(-(dop - P.dop_step2)./w2));  % 0->1

    base = P.low + P.h1 .* s1 + P.h2 .* s2;

    % optional Gaussian bump (usually sits on the top plateau)
    peak = P.peakAmp .* exp(-0.5*((dop - P.dop_peak)./sig).^2);

    y0 = base + peak;

    % right decay: sigmoid-gated exponential drop
    if P.decayAmp > 0
        dOn = 1 ./ (1 + exp(-(dop - P.dop_decay)./wd)); % 0 before onset, ->1 after
        tail = exp(-max(dop - P.dop_decay, 0)./decayLn); % =1 before onset, decays after

        % Mix: before onset => factor ~1; after onset => factor -> (1 - decayAmp) + decayAmp*tail
        decayFactor = (1 - P.decayAmp).*dOn + (1 - dOn) + P.decayAmp.*dOn.*tail;
        y = y0 .* decayFactor;
    else
        y = y0;
    end

    y = max(y, 0);
end

function S = fill_defaults_(S, def)
    f = fieldnames(def);
    for i = 1:numel(f)
        if ~isfield(S,f{i}) || isempty(S.(f{i}))
            S.(f{i}) = def.(f{i});
        end
    end
end