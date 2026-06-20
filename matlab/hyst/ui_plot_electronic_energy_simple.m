function ui_plot_electronic_energy_simple()
% UI plot for electronic energy built by build_electronic_energy_poly(P)
% Requirements:
%   - Plot Fel(psi1) AND dFel(psi1)
%   - X-range is EXACTLY par.range.psi1 (no L slider)
%   - T fixed from par.T (read from make_params)
%   - Slider: dop in par.range.dop
%   - Drag: update curves fast (skip minima); Release: compute minima

    % ---------------- initial params ----------------
    par = make_params();   % expects par.el + par.T + par.range.psi1/par.range.dop
    assert(isfield(par,'el'), 'make_params() must return par.el struct.');
    assert(isfield(par,'range') && isfield(par.range,'psi1'), 'Need par.range.psi1 = [xmin xmax].');
    assert(isfield(par,'range') && isfield(par.range,'dop'),  'Need par.range.dop  = [dmin dmax].');
    assert(isfield(par,'T'), 'Need par.T (fixed T used in UI).');

    needEl = {'a1_invV','a1_scaler','b1','T_eps'};
    for k=1:numel(needEl)
        assert(isfield(par.el,needEl{k}), 'Missing par.el.%s in make_params()', needEl{k});
    end

    psi1_lim = par.range.psi1;
    dop_lim  = par.range.dop;
    N0 = 3001;

    % --------- state (single source of truth) ---------
    state = struct();
    state.par = par;
    state.dop = mean(dop_lim);   % current doping point

    % -------------------- UI -------------------------
    fig = uifigure('Name','Electronic F(\psi_1) + dF UI', 'Position',[80 80 1250 720]);

    axF  = uiaxes(fig, 'Position',[360 380 850 300]);
    axdF = uiaxes(fig, 'Position',[360  80 850 260]);

    for ax = [axF axdF]
        ax.Box = 'on';
        grid(ax,'on');
        ax.TickLabelInterpreter = 'latex';
        xlim(ax, psi1_lim);
    end
    xlabel(axdF,'$\psi_1$','Interpreter','latex');
    ylabel(axF,'$F_{\mathrm{el}}(\psi_1)$','Interpreter','latex');
    ylabel(axdF,'$dF_{\mathrm{el}}/d\psi_1$','Interpreter','latex');

    title(axF,'Electronic energy (poly)','Interpreter','latex','FontWeight','normal');
    title(axdF,'Derivative','Interpreter','latex','FontWeight','normal');

    info = uilabel(fig, 'Position',[360 690 850 22], ...
                   'FontName','Consolas','FontSize',12);

    txt = uitextarea(fig, 'Position',[360 20 850 50], ...
                     'FontName','Consolas', 'Editable','off');

    % ---- slider row geometry ----
    y0=640;
    xName=20; xSld=120; wSld=190;
    xVal=320;

    % --- DOP slider only ---
    uilabel(fig,'Text','dop', 'Position',[xName y0 90 22], 'FontName','Consolas');

    sldDop = uislider(fig, 'Position',[xSld y0+10 wSld 3], ...
                      'Limits',dop_lim, 'Value',state.dop);
    sldDop.MajorTicks = linspace(dop_lim(1), dop_lim(2), 5);
    sldDop.MinorTicks = [];

    valDop = uilabel(fig,'Text',fmt_(state.dop), ...
                     'Position',[xVal y0 160 22], ...
                     'FontName','Consolas','HorizontalAlignment','left');

    % fixed T label
    uilabel(fig,'Text',sprintf('T (fixed) = %s', fmt_(par.T)), ...
        'Position',[20 y0-28 320 22], 'FontName','Consolas');

    uibutton(fig,'Text','Reset', 'Position',[20 20 90 34], ...
        'ButtonPushedFcn', @(~,~)resetAll_());

    % plot handles
    hF   = plot(axF,  nan, nan, 'LineWidth',2); hold(axF,'on');
    hMin = plot(axF,  nan, nan, 'ro', 'MarkerSize',7,'LineWidth',1.5,'MarkerFaceColor','r');
    hold(axF,'off');

    hdF  = plot(axdF, nan, nan, 'LineWidth',2); hold(axdF,'on');
    yline(axdF, 0, '--', 'LineWidth',1.0);
    hold(axdF,'off');

    % callbacks
    sldDop.ValueChangingFcn = @(~,evt) onChangingDop_(evt.Value);
    sldDop.ValueChangedFcn  = @(src,~)  onChangedDop_(src.Value);

    % initial draw: include minima
    redraw_(true);

    % =================== nested helpers ===================

    function onChangingDop_(v)
        state.dop = v;
        valDop.Text = fmt_(v);
        redraw_(false); % fast
    end

    function onChangedDop_(v)
        state.dop = v;
        valDop.Text = fmt_(v);
        redraw_(true); % with minima
    end

    function resetAll_()
        state.par = par;
        state.dop = mean(dop_lim);
        sldDop.Value = state.dop;
        valDop.Text  = fmt_(state.dop);
        redraw_(true);
    end

    function redraw_(doMinima)
        P = state.par;
        T = P.T;
        dop = state.dop;

        % ---------- chi(T,dop) ----------
        try
            chiFun = make_chi_artificial(P, struct()); % use defaults unless you pass Pchi elsewhere
            chiTD  = chiFun(T, dop);
        catch ME
            hF.XData = nan; hF.YData = nan;
            hdF.XData = nan; hdF.YData = nan;
            hMin.XData = nan; hMin.YData = nan;
            info.Text = "ERROR: make_chi_artificial not found / failed.";
            txt.Value = {ME.message};
            drawnow;
            return;
        end

        % ---------- build electronic poly energy ----------
        a1 = P.el.a1_scaler * (P.el.a1_invV - abs(chiTD));
        b1 = P.el.b1;

        Pel = struct('a1', a1, 'b1', b1);

        try
            el = build_electronic_energy_poly(Pel);
        catch ME
            hF.XData = nan; hF.YData = nan;
            hdF.XData = nan; hdF.YData = nan;
            hMin.XData = nan; hMin.YData = nan;
            info.Text = "ERROR: build_electronic_energy_poly failed.";
            txt.Value = {ME.message};
            drawnow;
            return;
        end

        % ---------- curves ----------
        psi = linspace(psi1_lim(1), psi1_lim(2), N0);
        F   = el.F(psi);
        dF  = el.dF(psi);

        hF.XData  = psi; hF.YData  = F;
        hdF.XData = psi; hdF.YData = dF;

        info.Text = sprintf('T=%s   dop=%s   chi=%s   a1=%s   b1=%s   X=[%s,%s]   minima=%s', ...
            fmt_(T), fmt_(dop), fmt_(chiTD), fmt_(a1), fmt_(b1), ...
            fmt_(psi1_lim(1)), fmt_(psi1_lim(2)), tern_(doMinima,'ON','OFF'));

        if ~doMinima
            hMin.XData = nan; hMin.YData = nan;
            txt.Value = {'(dragging) minima calculation skipped'};
            drawnow limitrate;
            return;
        end

        % ---- minima from dF sign-crossing (- -> +) ----
        sgn = sign(dF);
        sgn(~isfinite(sgn)) = 0;
        idx = find(sgn(1:end-1)<0 & sgn(2:end)>0) + 1;

        psiMin = [];
        FMin   = [];

        for ii=1:numel(idx)
            j  = idx(ii);
            jL = max(1, j-25);
            jR = min(numel(psi), j+25);
            a = psi(jL); b = psi(jR);
            try
                xm = fminbnd(@(x) el.F(x), a, b);
                psiMin(end+1) = xm; %#ok<AGROW>
                FMin(end+1)   = el.F(xm); %#ok<AGROW>
            catch
            end
        end

        % merge duplicates
        if ~isempty(psiMin)
            [psiMin, ord] = sort(psiMin);
            FMin = FMin(ord);
            keep = true(size(psiMin));
            tol = 5e-6 * (psi1_lim(2)-psi1_lim(1));
            for i2=2:numel(psiMin)
                if abs(psiMin(i2)-psiMin(i2-1)) < tol
                    keep(i2)=false;
                end
            end
            psiMin = psiMin(keep);
            FMin   = FMin(keep);
        end

        if isempty(psiMin)
            hMin.XData = nan; hMin.YData = nan;
            txt.Value = {'(no minima found in range)'};
        else
            hMin.XData = psiMin;
            hMin.YData = FMin;

            lines = cell(1, min(14, numel(psiMin)));
            for i3=1:numel(lines)
                lines{i3} = sprintf('%2d: psi=% .6f   Fmin=% .6e', i3, psiMin(i3), FMin(i3));
            end
            txt.Value = lines;
        end

        drawnow limitrate;
    end
end

% ---------------- small utils ----------------
function s = fmt_(x)
    if x == 0
        s = "0";
    elseif abs(x) < 1e-3 || abs(x) > 1e3
        s = string(sprintf('%.3e', x));
    else
        s = string(sprintf('%.4g', x));
    end
end

function out = tern_(cond, a, b)
    if cond, out = a; else, out = b; end
end