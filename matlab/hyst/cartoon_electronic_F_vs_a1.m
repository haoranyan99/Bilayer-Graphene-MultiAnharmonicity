function plot_electronic_free_energy_a1()

    % ---------------- parameters ----------------
    b1 = 1.0;
    a1_list = linspace(1, -1, 10);     % different a1 values
    phi = linspace(-6, 6, 1200);

    % ---------------- figure ----------------
    figure('Color','w','Position',[200 120 360 340]);
    ax = axes;
    hold(ax,'on')

    cmap = jet(256);
    colormap(cmap)

    a1_min = min(a1_list);
    a1_max = max(a1_list);

    % ---------------- plot curves ----------------
    for i = 1:numel(a1_list)

        a1 = a1_list(i);

        P.a1 = a1;
        P.b1 = b1;
        el = build_electronic_energy_poly(P);

        F = el.F(phi);

        % map a1 -> color
        t = (a1 - a1_min) / (a1_max - a1_min);
        idx = max(1, round(t*(size(cmap,1)-1))+1);
        color = cmap(idx,:);

        plot(phi, F, 'LineWidth', 1.2, 'Color', color)

    end

    % ---------------- axis style ----------------
    box on
    axis square

    set(gca,...
        'FontSize',16,...
        'LineWidth',1.0,...
        'XTick',[],...
        'YTick',[])

    xlabel('\varphi','FontSize',20,'Interpreter','tex')
    ylabel('Free energy $F$','FontSize',18,'Interpreter','latex')

    % ---------------- colorbar ----------------
    caxis([a1_min a1_max])
    cb = colorbar;
    cb.Label.String = 'a_1';
    cb.Label.FontSize = 16;

end