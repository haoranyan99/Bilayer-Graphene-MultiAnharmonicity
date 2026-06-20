function out = plot_lattice_potential_simple_current_posonly()
% plot_lattice_potential_simple_current_posonly
% Non-UI script:
%   - x-domain: ONLY psi2 > 0 (and within par.range.psi2)
%   - show ONLY positive part of curves: V>0 and dV>0
%   - mark minima of V on V-curve
%   - mark maxima of dV on dV-curve
%
% Output:
%   out.fig : figure handle
%   out.png : saved png path ("" if not saved)

    out = struct('fig',[],'png',"");

    % --------- load params ----------
    par = make_params();   % expects par.lat and par.range.psi2
    assert(isfield(par,'lat'), 'make_params() must return par.lat struct.');
    assert(isfield(par,'range') && isfield(par.range,'psi2'), 'Need par.range.psi2 = [xmin xmax].');

    need = {'A1','m','k4','s','eps'};
    for k=1:numel(need)
        assert(isfield(par.lat,need{k}), 'Missing par.lat.%s in make_params()', need{k});
    end

    psi2_lim = par.range.psi2;   % exact range
    N = 4001;

    % --------- restrict x to psi>0 ----------
    xL = max(psi2_lim(1), 0);
    xR = psi2_lim(2);

    if ~(xR > xL)
        error('Positive x-domain is empty: par.range.psi2=[%g,%g].', psi2_lim(1), psi2_lim(2));
    end

    % --------- build functions ----------
    lat = build_lattice_potential_simple(par);

    psi = linspace(xL, xR, N);
    V   = lat.V(psi);
    dV  = lat.dV(psi);

    % --------- show only positive y parts ----------
    V_pos  = V;   V_pos(~(V_pos > 0))   = NaN;
    dV_pos = dV;  dV_pos(~(dV_pos > 0)) = NaN;

    % --------- find minima of V (within psi>0) ----------
    % minima when dV crosses 0 from (-) to (+)
    sgn = sign(dV);
    sgn(~isfinite(sgn)) = 0;
    idxMinSeed = find(sgn(1:end-1) < 0 & sgn(2:end) > 0) + 1;

    psiMin = [];
    VMin   = [];
    for ii=1:numel(idxMinSeed)
        j  = idxMinSeed(ii);
        jL = max(1, j-25);
        jR = min(numel(psi), j+25);
        a = psi(jL); b = psi(jR);

        try
            xm = fminbnd(@(x) lat.V(x), a, b);
            vm = lat.V(xm);
            % keep only those that will be displayed (psi>0 already; require V>0)
            if isfinite(vm) && (vm > 0) && (xm > 0)
                psiMin(end+1) = xm; %#ok<AGROW>
                VMin(end+1)   = vm; %#ok<AGROW>
            end
        catch
        end
    end

    % dedup close minima
    [psiMin, VMin] = dedupPairs_(psiMin, VMin, 5e-3);

    % --------- find maxima of dV (within psi>0) ----------
    % Use numerical second derivative ddV sign change (+ -> -) for peaks
    h = psi(2) - psi(1);
    ddV = gradient(dV, h);  % ~ d^2V/dpsi^2
    idxPkSeed = find(ddV(1:end-1) > 0 & ddV(2:end) < 0) + 1;

    psiPk = [];
    dVPk  = [];
    for ii=1:numel(idxPkSeed)
        j  = idxPkSeed(ii);
        jL = max(1, j-25);
        jR = min(numel(psi), j+25);
        a = psi(jL); b = psi(jR);

        try
            xp = fminbnd(@(x) -lat.dV(x), a, b);
            vp = lat.dV(xp);
            % keep only displayed part: dV>0 and psi>0
            if isfinite(vp) && (vp > 0) && (xp > 0)
                psiPk(end+1) = xp; %#ok<AGROW>
                dVPk(end+1)  = vp; %#ok<AGROW>
            end
        catch
        end
    end

    % dedup close peaks
    [psiPk, dVPk] = dedupPairs_(psiPk, dVPk, 5e-3);

    % --------- plot ----------
    fig = figure('Color','w','Position',[120 120 980 760]);
    out.fig = fig;

    ax1 = subplot(2,1,1);
    plot(ax1, psi, V_pos, 'LineWidth', 2);
    hold(ax1,'on');
    if ~isempty(psiMin)
        plot(ax1, psiMin, VMin, 'ro', 'MarkerFaceColor','r', 'LineWidth', 1.5);
        for i=1:numel(psiMin)
            text(ax1, psiMin(i), VMin(i), sprintf('  min @ %.4g', psiMin(i)), ...
                'FontSize', 10, 'Interpreter','none', 'VerticalAlignment','bottom');
        end
    end
    hold(ax1,'off');
    grid(ax1,'on'); box(ax1,'on');
    xlim(ax1, [xL xR]);
    ylabel(ax1,'$V(\psi_2)$','Interpreter','latex');
    title(ax1, sprintf('V(\\psi_2) (x>0, y>0 only) | A1=%g, m=%g, k4=%g, s=%g, eps=%g', ...
        par.lat.A1, par.lat.m, par.lat.k4, par.lat.s, par.lat.eps), ...
        'Interpreter','none','FontWeight','normal');

    ax2 = subplot(2,1,2);
    plot(ax2, psi, dV_pos, 'LineWidth', 2);
    hold(ax2,'on');
    yline(ax2, 0, '--', 'LineWidth', 1.0);
    if ~isempty(psiPk)
        plot(ax2, psiPk, dVPk, 'ko', 'MarkerFaceColor','k', 'LineWidth', 1.2);
        for i=1:numel(psiPk)
            text(ax2, psiPk(i), dVPk(i), sprintf('  max @ %.4g', psiPk(i)), ...
                'FontSize', 10, 'Interpreter','none', 'VerticalAlignment','bottom');
        end
    end
    hold(ax2,'off');
    grid(ax2,'on'); box(ax2,'on');
    xlim(ax2, [xL xR]);
    xlabel(ax2,'$\psi_2$','Interpreter','latex');
    ylabel(ax2,'$dV/d\psi_2$','Interpreter','latex');
    title(ax2, 'dV/d\psi_2 (x>0, y>0 only) + maxima markers', ...
        'Interpreter','none','FontWeight','normal');

    % --------- optional save ----------
    % out.png = "lattice_V_dV_posonly_marked.png";
    % exportgraphics(fig, out.png, 'Resolution', 200);

end

% ================= helpers =================
function [x2,y2] = dedupPairs_(x, y, tol)
    if isempty(x)
        x2 = []; y2 = [];
        return;
    end
    [x, ord] = sort(x);
    y = y(ord);

    keep = true(size(x));
    for i=2:numel(x)
        if abs(x(i)-x(i-1)) < tol
            keep(i) = false;
        end
    end
    x2 = x(keep);
    y2 = y(keep);
end