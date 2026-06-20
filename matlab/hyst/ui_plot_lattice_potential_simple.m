function ui_plot_lattice_potential_simple()
% UI plot for lattice potential built by build_lattice_potential_simple(P_lat)
% Requirements:
%   - Plot V(psi2) AND dV(psi2)
%   - X-range is EXACTLY par.range.psi2 (no L slider)
%   - Sliders: A1, m, k4, s, eps
%   - Drag: update curves fast (skip minima); Release: compute minima

    % ---------------- initial params ----------------
    par = make_params();   % expects par.lat.A1/m/k4/s/eps and par.range.psi2
    assert(isfield(par,'lat'), 'make_params() must return par.lat struct.');
    assert(isfield(par,'range') && isfield(par.range,'psi2'), 'Need par.range.psi2 = [xmin xmax].');

    need = {'A1','m','k4','s','eps'};
    for k=1:numel(need)
        assert(isfield(par.lat,need{k}), 'Missing par.lat.%s in make_params()', need{k});
    end

    psi2_lim = par.range.psi2;
    N0 = 3001;   % curve resolution (increase if you want smoother)

    % --------- state (single source of truth) ---------
    state = struct();
    state.par = par;

    % -------------------- UI -------------------------
    fig = uifigure('Name','Lattice V(\psi_2) + dV UI', 'Position',[80 80 1250 720]);

    % two axes: V (top), dV (bottom)
    axV  = uiaxes(fig, 'Position',[360 380 850 300]);
    axdV = uiaxes(fig, 'Position',[360  80 850 260]);

    for ax = [axV axdV]
        ax.Box = 'on';
        grid(ax,'on');
        ax.TickLabelInterpreter = 'latex';
        xlim(ax, psi2_lim);
    end
    xlabel(axdV,'$\psi_2$','Interpreter','latex');
    ylabel(axV,'$V(\psi_2)$','Interpreter','latex');
    ylabel(axdV,'$dV/d\psi_2$','Interpreter','latex');

    title(axV,'Lattice potential','Interpreter','latex','FontWeight','normal');
    title(axdV,'Derivative','Interpreter','latex','FontWeight','normal');

    info = uilabel(fig, 'Position',[360 690 850 22], ...
                   'FontName','Consolas','FontSize',12);

    txt = uitextarea(fig, 'Position',[360 20 850 50], ...
                     'FontName','Consolas', 'Editable','off');

    % ---- sliders panel geometry ----
    y0=640; dy=62;
    xName=20; xSld=120; wSld=190;
    xVal=320;

    rows = struct();
    rows.A1  = addRow_('A1',  par.lat.A1,  0,    10,   1);
    rows.m   = addRow_('m',   par.lat.m,   0,    10,   2);
    rows.k4  = addRow_('k4',  par.lat.k4,  0,    0.01, 3);
    rows.s   = addRow_('s',   par.lat.s,   0,    5.0,  4);
    rows.eps = addRow_('eps', par.lat.eps, 1e-4, 0.2,  5);

    uibutton(fig,'Text','Reset', 'Position',[20 20 90 34], ...
        'ButtonPushedFcn', @(~,~)resetAll_());

    % plot handles
    hV   = plot(axV,  nan, nan, 'LineWidth',2); hold(axV,'on');
    hMin = plot(axV,  nan, nan, 'ro', 'MarkerSize',7,'LineWidth',1.5,'MarkerFaceColor','r');
    hold(axV,'off');

    hdV  = plot(axdV, nan, nan, 'LineWidth',2); hold(axdV,'on');
    hz   = yline(axdV, 0, '--', 'LineWidth',1.0); %#ok<NASGU>
    hold(axdV,'off');

    attachCallbacks_();

    % initial draw: include minima
    redraw_(true);

    % =================== nested helpers ===================

    function row = addRow_(name, init, vmin, vmax, rowIdx)
        y = y0 - (rowIdx-1)*dy;

        uilabel(fig,'Text',name, 'Position',[xName y 90 22], 'FontName','Consolas');

        sld = uislider(fig, 'Position',[xSld y+10 wSld 3], ...
                       'Limits',[vmin vmax], 'Value',init);
        sld.MajorTicks = linspace(vmin, vmax, 5);
        sld.MinorTicks = [];

        val = uilabel(fig,'Text',fmt_(init), ...
                      'Position',[xVal y 140 22], ...
                      'FontName','Consolas', 'HorizontalAlignment','left');

        row = struct('sld',sld,'val',val,'name',name);
    end

    function attachCallbacks_()
        flds = fieldnames(rows);
        for ii=1:numel(flds)
            field = flds{ii};
            rows.(field).sld.ValueChangingFcn = @(~,evt) onChanging_(field, evt.Value);
            rows.(field).sld.ValueChangedFcn  = @(src,~)  onChanged_(field, src.Value);
        end
    end

    function onChanging_(field, v)
        rows.(field).val.Text = fmt_(v);
        state.par.lat.(field) = v;   % ✅ correct path
        redraw_(false);              % fast (skip minima)
    end

    function onChanged_(field, v)
        rows.(field).val.Text = fmt_(v);
        state.par.lat.(field) = v;
        redraw_(true);               % full (with minima)
    end

    function resetAll_()
        % reset slider values + labels
        for k=1:numel(need)
            f = need{k};
            rows.(f).sld.Value = par.lat.(f);
            rows.(f).val.Text  = fmt_(par.lat.(f));
        end
        % reset state
        state.par = par;
        redraw_(true);
    end

    function redraw_(doMinima)
        P = state.par;

        % builder takes .lat struct
        lat = build_lattice_potential_simple(P);

        psi = linspace(psi2_lim(1), psi2_lim(2), N0);
        V   = lat.V(psi);
        dV  = lat.dV(psi);

        % update curves
        hV.XData  = psi;  hV.YData  = V;
        hdV.XData = psi;  hdV.YData = dV;

        info.Text = sprintf('A1=%s   m=%s   k4=%s   s=%s   eps=%s   X=[%s,%s]   minima=%s', ...
            fmt_(P.lat.A1), fmt_(P.lat.m), fmt_(P.lat.k4), fmt_(P.lat.s), fmt_(P.lat.eps), ...
            fmt_(psi2_lim(1)), fmt_(psi2_lim(2)), tern_(doMinima,'ON','OFF'));

        if ~doMinima
            % fast path: hide minima marker
            hMin.XData = nan; hMin.YData = nan;
            txt.Value = {'(dragging) minima calculation skipped'};
            drawnow limitrate;
            return;
        end

        % ---- compute minima from dV sign-crossing (- -> +) ----
        sgn = sign(dV);
        sgn(~isfinite(sgn)) = 0;
        idx = find(sgn(1:end-1)<0 & sgn(2:end)>0) + 1;

        psiMin = [];
        VMin   = [];

        for ii=1:numel(idx)
            j  = idx(ii);
            jL = max(1, j-25);
            jR = min(numel(psi), j+25);
            a = psi(jL); b = psi(jR);
            try
                xm = fminbnd(@(x) lat.V(x), a, b);
                psiMin(end+1) = xm; %#ok<AGROW>
                VMin(end+1)   = lat.V(xm); %#ok<AGROW>
            catch
            end
        end

        % merge duplicates
        if ~isempty(psiMin)
            [psiMin, ord] = sort(psiMin);
            VMin = VMin(ord);
            keep = true(size(psiMin));
            tol = 5e-3;
            for i2=2:numel(psiMin)
                if abs(psiMin(i2)-psiMin(i2-1)) < tol
                    keep(i2)=false;
                end
            end
            psiMin = psiMin(keep);
            VMin   = VMin(keep);
        end

        % update minima markers
        if isempty(psiMin)
            hMin.XData = nan; hMin.YData = nan;
            txt.Value = {'(no minima found in range)'};
        else
            hMin.XData = psiMin;
            hMin.YData = VMin;

            % monotonic check by |psi|
            [absPsi, ordAbs] = sort(abs(psiMin));
            Vabs = VMin(ordAbs);
            bad = find(diff(Vabs) < -1e-6);

            lines = cell(1, min(14, numel(absPsi)));
            for i3=1:numel(lines)
                lines{i3} = sprintf('%2d: psi=%8.4f   |psi|=%8.4f   Vmin=% .6e', ...
                    i3, psiMin(ordAbs(i3)), absPsi(i3), Vabs(i3));
            end
            if ~isempty(bad)
                lines{end+1} = sprintf('NOT monotone at steps: %s', mat2str(bad(:)'));
            else
                lines{end+1} = 'OK: Vmin increases with |psi| (within tolerance).';
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