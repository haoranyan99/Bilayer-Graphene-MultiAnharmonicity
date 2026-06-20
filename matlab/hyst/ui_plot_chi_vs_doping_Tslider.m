function ui_plot_chi_vs_doping_Tslider()
% UI: plot chi(T,dop) vs doping with a Temperature slider.
% chi is provided by make_chi_artificial(par, Pchi).
%
% Controls:
%   - T slider updates curve in real time
%   - doping range comes from par.range.dop (fallback provided)

    % ---------------- params ----------------
    par = make_params();
    if ~isfield(par,'range') || ~isfield(par.range,'dop')
        par.range.dop = [1.6, 1.95];
    end
    if ~isfield(par,'range') || ~isfield(par.range,'T')
        if ~isfield(par,'T'), par.T = 6; end
        par.range.T = [0, 2*par.T];
    end
    if ~isfield(par,'el') || ~isfield(par.el,'T_eps')
        par.el.T_eps = 1e-6;
    end
    if ~isfield(par,'T')
        par.T = 6;
    end

    dop_lim = par.range.dop;
    T_lim   = par.range.T;
    if isempty(T_lim) || numel(T_lim)~=2
        T_lim = [0, 2*par.T];
    end

    N = 1500;
    dop = linspace(dop_lim(1), dop_lim(2), N);

    % ---------------- chi handle (single source of truth) ----------------
    % Optionally override chi-shape parameters here:
    Pchi = struct();
    % Pchi.dop_step  = 1.72;
    % Pchi.w_step0   = 1e-3;
    % Pchi.plateau   = 1.0;
    % Pchi.dop_peak  = 1.86;
    % Pchi.sig0      = 1e-3;
    % Pchi.peakAmp   = 0.9;
    % Pchi.dop_decay = 1.90;
    % Pchi.decayLen0 = 0.10;
    % Pchi.decayAmp  = 0.25;

    chi = make_chi_artificial(par, Pchi);

    % ---------------- UI ----------------
    fig = uifigure('Name','chi(T,doping) slider', 'Position',[100 100 1050 640]);

    ax = uiaxes(fig, 'Position',[80 160 900 440]);
    ax.Box = 'on';
    grid(ax,'on');
    ax.TickLabelInterpreter = 'latex';
    xlabel(ax,'doping','Interpreter','latex');
    ylabel(ax,'$\chi(T,\mathrm{dop})$','Interpreter','latex');
    title(ax,'$\chi$ vs doping (T slider)','Interpreter','latex','FontWeight','normal');
    xlim(ax, dop_lim);

    info = uilabel(fig, 'Position',[80 610 900 22], ...
        'FontName','Consolas','FontSize',12);

    % slider row
    uilabel(fig,'Text','T', 'Position',[80 95 30 22], 'FontName','Consolas');
    sld = uislider(fig, 'Position',[120 105 700 3], ...
        'Limits', T_lim, 'Value', clamp_(par.T, T_lim(1), T_lim(2)));
    sld.MajorTicks = linspace(T_lim(1), T_lim(2), 6);
    sld.MinorTicks = [];

    valT = uilabel(fig,'Text',fmt_(sld.Value), ...
        'Position',[840 95 140 22], 'FontName','Consolas');

    % plot handle
    h = plot(ax, nan, nan, 'LineWidth', 2);

    % callbacks
    sld.ValueChangingFcn = @(~,evt) onT_(evt.Value, true);
    sld.ValueChangedFcn  = @(src,~)  onT_(src.Value, false);

    % initial draw
    onT_(sld.Value, false);

    % ================= nested =================
    function onT_(T, isDragging)
        valT.Text = fmt_(T);

        y = chi(T, dop);
        h.XData = dop;
        h.YData = y;

        info.Text = sprintf('T=%s   dop=[%s,%s]   %s', ...
            fmt_(T), fmt_(dop_lim(1)), fmt_(dop_lim(2)), tern_(isDragging,'(dragging)',''));
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

function v = clamp_(x, a, b)
    v = min(max(x, a), b);
end