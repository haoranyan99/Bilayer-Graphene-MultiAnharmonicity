function freeEnergy_scape_dop_lambda_UI_new()
% UI scape with minimum (red dot)
% - main sliders: doping, lambda
% - control window: invV, scaler, T
% - b1 fixed = par.el.b1 (no slider)
% - lattice fixed here (no UI)
% - ranges from par.range.*

    % -------------------- load params --------------------
    par = make_params();


    % -------------------- ranges --------------------
    psi1_lim = par.range.psi1;
    psi2_lim = par.range.psi2;

    dop_min  = par.range.dop(1);
    dop_max  = par.range.dop(2);

    lam_min  = par.range.lambda(1);
    lam_max  = par.range.lambda(2);

    T_min    = par.range.T(1);
    T_max    = par.range.T(2);

    invV_min = par.range.inV(1);
    invV_max = par.range.inV(2);

    sc_min   = par.range.scaler(1);
    sc_max   = par.range.scaler(2);

    % -------------------- state --------------------
    coef = make_coeff(par);

    psi1_opt = 0;
    psi2_opt = 0;

    T_state = par.T; % controlled by 2nd window

    % grids
    Nsurf = 161;
    psi1_vals = linspace(psi1_lim(1), psi1_lim(2), Nsurf);
    psi2_vals = linspace(psi2_lim(1), psi2_lim(2), Nsurf);
    [P1,P2]   = meshgrid(psi1_vals, psi2_vals);

    % -------------------- main figure --------------------
    fig = uifigure('Name','Free Energy scape (with minimum)', ...
        'Position',[80 80 1020 610]);

    ax2d = uiaxes(fig,'Position',[50 200 540 400],'FontSize',12);
    ax2d.TickLabelInterpreter = 'latex';
    xlabel(ax2d,'$\psi_1$','Interpreter','latex');
    ylabel(ax2d,'$\psi_2$','Interpreter','latex');
    xlim(ax2d, psi1_lim); ylim(ax2d, psi2_lim);

    axPsi = uiaxes(fig,'Position',[620 410 360 190],'FontSize',12);
    axPsi.TickLabelInterpreter = 'latex';
    xlabel(axPsi,'$\psi_1$','Interpreter','latex');
    ylabel(axPsi,'$F_\psi(\psi_1)$','Interpreter','latex');
    title(axPsi,'Electronic part','Interpreter','latex','FontWeight','normal');

    axX = uiaxes(fig,'Position',[620 200 360 190],'FontSize',12);
    axX.TickLabelInterpreter = 'latex';
    xlabel(axX,'$\psi_2$','Interpreter','latex');
    ylabel(axX,'$F_X(\psi_2)$','Interpreter','latex');
    title(axX,'Lattice part (fixed)','Interpreter','latex','FontWeight','normal');

    info_box = uitextarea(fig,'Position',[620 40 360 140], ...
        'Editable','off','FontSize',12,'FontName','Consolas');

    % -------------------- sliders: doping & lambda --------------------
    y0 = 110; dy = 50;

    dop0 = 0.5*(dop_min + dop_max);
    dop_slider = uislider(fig,'Position',[150 y0+dy 340 3], ...
        'Limits',[dop_min, dop_max], 'Value', dop0);
    dop_slider.MajorTicks = linspace(dop_min, dop_max, 5);
    dop_slider.MinorTicks = [];

    lam0 = par.lambda;
    lam_slider = uislider(fig,'Position',[150 y0 340 3], ...
        'Limits',[lam_min, lam_max], 'Value', lam0);
    lam_slider.MajorTicks = linspace(lam_min, lam_max, 6);
    lam_slider.MinorTicks = [];

    uilabel(fig,'Position',[70 y0+dy-10 70 22],'Text','doping','FontSize',12);
    uilabel(fig,'Position',[70 y0-10    70 22],'Text','$\lambda$','FontSize',12,'Interpreter','latex');

    dop_text = uilabel(fig,'Position',[500 y0+dy-14 260 22],'Text','', 'FontSize',12,'FontName','Consolas');
    lam_text = uilabel(fig,'Position',[500 y0-14    260 22],'Text','', 'FontSize',12,'FontName','Consolas');

    dop_slider.ValueChangingFcn = @(~,evt) set(dop_text,'Text',sprintf('doping = %.6g',evt.Value));
    lam_slider.ValueChangingFcn = @(~,evt) set(lam_text,'Text',sprintf('lambda = %.6g',evt.Value));

    % IMPORTANT: release => recompute minimum
    dop_slider.ValueChangedFcn = @(~,~) update_();
    lam_slider.ValueChangedFcn = @(~,~) update_();

    % -------------------- control window --------------------
    ctrl = make_control_window_();

    % initial
    refresh_labels_();
    update_plot_only_();  % draw surface first (no jump)
    update_();            % then compute minimum once

    % ===================== nested =====================
    function refresh_labels_()
        dop_text.Text = sprintf('doping = %.6g', dop_slider.Value);
        lam_text.Text = sprintf('lambda = %.6g', lam_slider.Value);
        ctrl.T_text.Text = sprintf('T = %.6g', T_state);
    end

    function apply_ctrl_to_par_()
        par.el.a1_invV   = ctrl.invV.Value;
        par.el.a1_scaler = ctrl.scaler.Value;
        % b1 fixed:
        par.el.b1        = par.el.b1;
        par.el.T_eps     = 1e-6;
        % lattice fixed already
    end

    function update_()
        refresh_labels_();
        apply_ctrl_to_par_();
        coef = make_coeff(par);

        dop = dop_slider.Value;
        lam = lam_slider.Value;

        el  = coef.build_el(T_state, dop);
        lat = coef.build_lat();

        % continuation + multistart => robust red dot movement
        try
            [psi1_opt, psi2_opt, Fmin] = minimize_multistart_( ...
                el, lat, lam, psi1_opt, psi2_opt, psi1_lim, psi2_lim);
        catch
            psi1_opt = 0; psi2_opt = 0; Fmin = NaN;
        end

        draw_(dop, lam, el, lat, Fmin);
    end

    function update_plot_only_()
        refresh_labels_();
        apply_ctrl_to_par_();
        coef = make_coeff(par);

        dop = dop_slider.Value;
        lam = lam_slider.Value;

        el  = coef.build_el(T_state, dop);
        lat = coef.build_lat();

        draw_(dop, lam, el, lat, NaN);
    end

    function draw_(dop, lam, el, lat, Fmin)
        cla(ax2d); cla(axPsi); cla(axX);

        % 2D surface
        F2D = @(p1,p2) el.F(p1) + lat.V(p2) + lam*p1.*p2;
        Fv  = arrayfun(F2D, P1, P2);

        pcolor(ax2d, P1, P2, Fv);
        shading(ax2d,'interp');
        colormap(ax2d,'turbo');
        colorbar(ax2d,'Location','eastoutside','FontSize',12);
        hold(ax2d,'on');
        contour(ax2d, P1, P2, Fv, 22, 'LineColor','k','LineWidth',0.5);

        % mark minimum
        if isfinite(psi1_opt) && isfinite(psi2_opt)
            plot(ax2d, psi1_opt, psi2_opt, 'ro','MarkerFaceColor','r','MarkerSize',7);
        end
        hold(ax2d,'off');
        xlim(ax2d, psi1_lim); ylim(ax2d, psi2_lim);
        title(ax2d, sprintf('T=%.4g, dop=%.6g, \\lambda=%.6g', T_state, dop, lam), ...
            'FontWeight','normal');

        % 1D cuts
        x1 = linspace(psi1_lim(1), psi1_lim(2), 801);
        x2 = linspace(psi2_lim(1), psi2_lim(2), 801);

        plot(axPsi, x1, el.F(x1), 'LineWidth',2);
        grid(axPsi,'on'); box(axPsi,'on');
        hold(axPsi,'on');
        if isfinite(psi1_opt), plot(axPsi, psi1_opt, el.F(psi1_opt), 'ro','MarkerFaceColor','r','MarkerSize',6); end
        hold(axPsi,'off');

        plot(axX, x2, lat.V(x2), 'LineWidth',2);
        grid(axX,'on'); box(axX,'on');
        hold(axX,'on');
        if isfinite(psi2_opt), plot(axX, psi2_opt, lat.V(psi2_opt), 'ro','MarkerFaceColor','r','MarkerSize',6); end
        hold(axX,'off');

        % info
        chi_val = NaN;
        if isfield(coef,'chi')
            try, chi_val = coef.chi(T_state, dop); catch, end
        end

        lines = {
            sprintf('T = %.12g (broadening)', T_state)
            sprintf('doping = %.12g', dop)
            sprintf('lambda = %.12g', lam)
            sprintf('chi(T,dop) = %.12g', chi_val)
            sprintf('el.a1_invV   = %.12g', par.el.a1_invV)
            sprintf('el.a1_scaler = %.12g', par.el.a1_scaler)
            sprintf('el.b1        = %.12g (fixed)', par.el.b1)
            sprintf('a1 = %.12g', el.params.a1)
            sprintf('b1 = %.12g', el.params.b1)
            sprintf('psi1* = %.12g', psi1_opt)
            sprintf('psi2* = %.12g', psi2_opt)
        };
        if isfinite(Fmin)
            lines{end+1} = sprintf('Fmin = %.12g', Fmin);
        end
        info_box.Value = lines;
    end

    function C = make_control_window_()
        fig2 = uifigure('Name','Controls (electronic only): invV / scaler / T', ...
            'Position',[1120 120 520 260]);

        C = struct();

        y0c=190; dyc=60;
        xName=20; xSld=210; wSld=220; xVal=450;

        [C.invV,   ~] = addRow_(fig2,'invV',   par.el.a1_invV,   invV_min, invV_max, 1);
        [C.scaler, ~] = addRow_(fig2,'scaler', par.el.a1_scaler, sc_min,   sc_max,   2);

        % T slider
        uilabel(fig2,'Text','T (broadening)', 'Position',[xName y0c-2*dyc 180 22], 'FontName','Consolas');
        C.T = uislider(fig2,'Position',[xSld y0c-2*dyc+10 wSld 3], ...
            'Limits',[T_min, T_max], 'Value', T_state);
        C.T.MajorTicks = linspace(T_min, T_max, 6);
        C.T.MinorTicks = [];

        C.T_text = uilabel(fig2,'Text',sprintf('T = %.6g',T_state), ...
            'Position',[xVal y0c-2*dyc 70 22], 'FontName','Consolas');

        % Drag => just redraw surface; Release => minimize
        C.T.ValueChangingFcn = @(~,evt) onTChanging_(evt.Value);
        C.T.ValueChangedFcn  = @(~,~) onTChanged_();

        function onTChanging_(v)
            C.T_text.Text = sprintf('T = %.6g', v);
            T_state = v;
            update_plot_only_();
        end
        function onTChanged_()
            T_state = C.T.Value;
            update_();
        end

        function [sld, val] = addRow_(parent, name, init, vmin, vmax, rowIdx)
            y = y0c - (rowIdx-1)*dyc;
            uilabel(parent,'Text',name, 'Position',[xName y 180 22], 'FontName','Consolas');

            sld = uislider(parent, 'Position',[xSld y+10 wSld 3], ...
                           'Limits',[vmin vmax], 'Value',init);
            sld.MajorTicks = linspace(vmin, vmax, 5);
            sld.MinorTicks = [];

            val = uilabel(parent,'Text',fmt_(init), ...
                          'Position',[xVal y 70 22], ...
                          'FontName','Consolas','HorizontalAlignment','left');

            sld.ValueChangingFcn = @(~,evt) onChanging_(evt.Value);
            sld.ValueChangedFcn  = @(~,~) onChanged_();

            function onChanging_(v)
                val.Text = fmt_(v);
                update_plot_only_();
            end
            function onChanged_()
                val.Text = fmt_(sld.Value);
                update_();
            end
        end
    end
end

% ===================== multistart minimizer =====================
function [p1_best, p2_best, F_best] = minimize_multistart_(el, lat, lambda, p1_seed, p2_seed, psi1_lim, psi2_lim)
    n_start = 12;

    starts = zeros(n_start,2);
    starts(1,:) = [p1_seed, p2_seed];

    for k=2:4
        starts(k,1) = p1_seed + 0.4*randn();
        starts(k,2) = p2_seed + 0.8*randn();
    end

    for k=5:n_start
        starts(k,1) = psi1_lim(1) + (psi1_lim(2)-psi1_lim(1))*rand();
        starts(k,2) = psi2_lim(1) + (psi2_lim(2)-psi2_lim(1))*rand();
    end

    starts(:,1) = min(max(starts(:,1), psi1_lim(1)), psi1_lim(2));
    starts(:,2) = min(max(starts(:,2), psi2_lim(1)), psi2_lim(2));

    F_best = +Inf;
    p1_best = p1_seed; p2_best = p2_seed;

    for k=1:n_start
        [p1,p2,Fk] = minimize_local_(el, lat, lambda, starts(k,1), starts(k,2));
        if isfinite(Fk) && Fk < F_best
            F_best = Fk;
            p1_best = p1;
            p2_best = p2;
        end
    end
end

% ===================== local minimizer (CG + line search) =====================
function [psi1_opt, psi2_opt, Fmin] = minimize_local_(el, lat, lambda, psi1_seed, psi2_seed)
    tol = 1e-6;
    max_iter = 600;

    F = @(p1,p2) el.F(p1) + lat.V(p2) + lambda*p1.*p2;
    gradF = @(p1,p2) [el.dF(p1) + lambda*p2; lat.dV(p2) + lambda*p1];

    psi1 = psi1_seed; psi2 = psi2_seed;
    g = gradF(psi1,psi2);
    d = -g;

    it = 0;
    while norm(g) > tol && it < max_iter
        alpha = line_search_local_(gradF, psi1, psi2, d);

        psi1n = psi1 + alpha*d(1);
        psi2n = psi2 + alpha*d(2);
        gn = gradF(psi1n,psi2n);

        denom = g'*g;
        if denom <= 0, beta = 0; else, beta = max((gn'*(gn-g))/denom, 0); end
        dn = -gn + beta*d;
        if dot(dn,-gn) < 0, dn = -gn; end

        psi1 = psi1n; psi2 = psi2n;
        g = gn; d = dn;
        it = it + 1;
    end

    psi1_opt = psi1;
    psi2_opt = psi2;
    Fmin = F(psi1,psi2);
end

function alpha = line_search_local_(gradF, psi1, psi2, d)
    tol = 1e-6;
    max_iter = 80;
    alpha1 = 0;
    alpha2 = 0.1;
    step_size = 0.1;

    g1 = dirderiv_(gradF, psi1, psi2, d, alpha1);
    if g1 >= 0
        alpha = 0; return;
    end

    g2 = dirderiv_(gradF, psi1, psi2, d, alpha2);
    guard = 0;
    while g2 <= 0
        alpha1 = alpha2;
        alpha2 = alpha2 + step_size;
        g2 = dirderiv_(gradF, psi1, psi2, d, alpha2);
        guard = guard + 1;
        if guard > 2000
            alpha = alpha2; return;
        end
    end

    for i=1:max_iter
        am = 0.5*(alpha1+alpha2);
        gm = dirderiv_(gradF, psi1, psi2, d, am);
        if abs(gm) < tol || abs(alpha2-alpha1) < tol
            alpha = am; return;
        end
        if gm < 0, alpha1 = am; else, alpha2 = am; end
    end
    alpha = 0.5*(alpha1+alpha2);
end

function g = dirderiv_(gradF, psi1, psi2, d, alpha)
    p1 = psi1 + alpha*d(1);
    p2 = psi2 + alpha*d(2);
    gg = gradF(p1,p2);
    g = gg' * d;
end

function s = fmt_(x)
    if x == 0
        s = "0";
    elseif abs(x) < 1e-3 || abs(x) > 1e3
        s = string(sprintf('%.3e', x));
    else
        s = string(sprintf('%.4g', x));
    end
end