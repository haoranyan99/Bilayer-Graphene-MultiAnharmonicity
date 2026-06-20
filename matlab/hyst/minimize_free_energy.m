function [psi1_opt, psi2_opt, F_opt, out] = minimize_free_energy(par, T, dop, psi1_init, psi2_init, is_plot)
% minimize_free_energy (par-driven)
%
% Minimizes total free energy:
%   F(psi1,psi2) = Fel(psi1) + Vlat(psi2) + lambda*psi1*psi2
%
% where
%   Fel is built from coef.build_el(T,dop)   (4th order only)
%   Vlat is built from coef.build_lat()      (your simplified lattice form)
%
% Inputs:
%   par        : from make_params()
%   T, dop     : scalars
%   psi1_init  : initial psi1
%   psi2_init  : initial psi2
%   is_plot    : true/false
%
% Outputs:
%   psi1_opt, psi2_opt, F_opt
%   out: struct with diagnostics (iters, avg line search iters, runtime, trajectory)

    tic;

    % ---------------- settings ----------------
    tol      = 1e-8;
    max_iter = 1000;
    % -----------------------------------------

    % build coefficient interfaces
    coef = make_coeff(par);
    el   = coef.build_el(T, dop);
    lat  = coef.build_lat();

    if isfield(par,'lambda')
        lambda = par.lambda;
    else
        lambda = 0.0;
    end

    % total free energy and gradient
    F = @(psi1,psi2) el.F(psi1) + lat.V(psi2) + lambda*psi1.*psi2;

    grad_F = @(psi1,psi2) [ ...
        el.dF(psi1) + lambda*psi2; ...
        lat.dV(psi2) + lambda*psi1 ];

    % init CG
    psi1 = psi1_init;
    psi2 = psi2_init;

    grad = grad_F(psi1, psi2);
    d    = -grad;

    iter_points = [psi1, psi2];
    total_bisect_iter = 0;

    iter = 0;
    while norm(grad) > tol && iter < max_iter

        % line search (bisection on directional derivative)
        [alpha, bisect_iter] = line_search_(grad_F, psi1, psi2, d);
        total_bisect_iter = total_bisect_iter + bisect_iter;

        % update
        psi1_new = psi1 + alpha * d(1);
        psi2_new = psi2 + alpha * d(2);

        grad_new = grad_F(psi1_new, psi2_new);

        % Polak–Ribière+ (with restart safeguard)
        denom = (grad' * grad);
        if denom <= 0
            beta = 0;
        else
            beta = max((grad_new' * (grad_new - grad)) / denom, 0);
        end
        d_new = -grad_new + beta * d;

        % restart if not a descent direction
        if dot(d_new, -grad_new) < 0
            d_new = -grad_new;
        end

        psi1 = psi1_new;
        psi2 = psi2_new;
        grad = grad_new;
        d    = d_new;

        iter_points(end+1,:) = [psi1, psi2]; %#ok<AGROW>
        iter = iter + 1;
    end

    psi1_opt = psi1;
    psi2_opt = psi2;
    F_opt    = F(psi1_opt, psi2_opt);

    cg_iter_count = iter;
    if cg_iter_count > 0
        avg_bisect_iter = total_bisect_iter / cg_iter_count;
    else
        avg_bisect_iter = 0;
    end

    out = struct();
    out.cg_iters = cg_iter_count;
    out.avg_line_search_iters = avg_bisect_iter;
    out.runtime_sec = toc;
    out.traj = iter_points;
    out.final_grad_norm = norm(grad);

    if is_plot
        fprintf('T=%.6g, dop=%.6g\n', T, dop);
        fprintf('psi1 = %.8g, psi2 = %.8g, F = %.12g\n', psi1_opt, psi2_opt, F_opt);
        fprintf('CG iters = %d, Avg line-search(bisect) iters = %.4g, final |grad|=%.4g, time=%.3fs\n', ...
            cg_iter_count, avg_bisect_iter, out.final_grad_norm, out.runtime_sec);

        % plot 2D landscape around trajectory
        pad = 1.0;
        psi1_vals = linspace(min(iter_points(:,1))-pad, max(iter_points(:,1))+pad, 80);
        psi2_vals = linspace(min(iter_points(:,2))-pad, max(iter_points(:,2))+pad, 80);
        [Psi1, Psi2] = meshgrid(psi1_vals, psi2_vals);
        F_vals = F(Psi1, Psi2);

        figure('Color','w');
        pcolor(Psi1, Psi2, F_vals);
        shading interp;
        colorbar;
        colormap turbo;
        hold on;
        contour(Psi1, Psi2, F_vals, 25, 'LineColor', 'k', 'LineWidth', 0.5);

        plot(iter_points(:,1), iter_points(:,2), 'r-o', 'LineWidth', 0.6, 'MarkerSize', 3);

        xlabel('\psi_1');
        ylabel('\psi_2');
        title('Free energy (par-driven) + CG trajectory');
    end
end


% ============================================================
% Line search by bisection on directional derivative:
%   g(alpha) = grad_F(psi+alpha*d) · d
% Find root where g(alpha)=0, starting from alpha=0 with g(0)<0.
% ============================================================
function [alpha, bisect_iter] = line_search_(grad_F, psi1, psi2, d)

    tol = 1e-6;
    max_iter = 100;

    alpha1 = 0;
    alpha2 = 0.1;
    step_size = 0.1;

    bisect_iter = 0;

    g1 = dirderiv_(grad_F, psi1, psi2, d, alpha1);
    if g1 >= 0
        % not a descent direction (should be prevented by restart, but just in case)
        alpha = 0;
        return;
    end

    g2 = dirderiv_(grad_F, psi1, psi2, d, alpha2);
    guard = 0;
    while g2 <= 0
        alpha1 = alpha2;
        alpha2 = alpha2 + step_size;
        g1 = g2;
        g2 = dirderiv_(grad_F, psi1, psi2, d, alpha2);
        guard = guard + 1;
        if guard > 2000
            % failsafe: return something small
            alpha = alpha2;
            return;
        end
    end

    for i = 1:max_iter
        bisect_iter = bisect_iter + 1;
        alpha_mid = 0.5*(alpha1 + alpha2);
        g_mid = dirderiv_(grad_F, psi1, psi2, d, alpha_mid);

        if abs(g_mid) < tol || abs(alpha2 - alpha1) < tol
            alpha = alpha_mid;
            return;
        end

        if g_mid < 0
            alpha1 = alpha_mid;
        else
            alpha2 = alpha_mid;
        end
    end

    alpha = 0.5*(alpha1 + alpha2);
end


function g = dirderiv_(grad_F, psi1, psi2, d, alpha)
    p1 = psi1 + alpha * d(1);
    p2 = psi2 + alpha * d(2);
    grad_new = grad_F(p1, p2);
    g = grad_new' * d;
end