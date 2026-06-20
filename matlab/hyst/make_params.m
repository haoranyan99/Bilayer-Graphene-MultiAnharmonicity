function par = make_params()
% make_params (minimal, current layout)
% Only provides parameters needed to build:
%   - electronic energy builder (a1 from chi(T,dop))
%   - lattice potential builder (simple form)
% PLUS provides ranges for plotting scape UI:
%   par.range.psi1, par.range.psi2, par.range.dop, par.range.lambda, par.range.T
%
% NEW:
%   par.hyst.refine    : doping grid refinement used in main_hyst computation
%   par.hyst.highlight : plotting/analysis region (independent of refine)

    par = struct();

    par.lambda = 1;     % coupling
    par.T = 6;          % temperature

    % ===================== electronic (psi1) =====================
    par.el = struct();
    par.el.a1_invV   =  0.80;
    par.el.a1_scaler =  170;
    par.el.b1 = 1.2;
    par.el.T_eps = 1e-6;

    % ===================== lattice (psi2) =====================
    par.lat = struct();
    par.lat.A1  = 5;
    par.lat.kd  = 1e-3;
    par.lat.m   = 0.2;
    par.lat.k4  = 1e-4;
    par.lat.s   = 2.0;
    par.lat.w   = 1.5;
    par.lat.eps = 0.1;
    par.lat.psi_cut = 25;

    % ===================== ranges for scape/UI =====================
    par.range = struct();

    chi_peak = 1.0;
    delta_h = chi_peak - par.el.a1_invV;
    estimate_psi1_0 = sqrt(6 * delta_h * par.el.a1_scaler / par.el.b1);
    
    dop_start = 1.7;
    dop_end = 1.8;

    par.range.psi1   = [-1.2*estimate_psi1_0, 1.2*estimate_psi1_0];
    par.range.psi2   = [-25, 25];
    par.range.dop    = [dop_start, dop_end];
    par.range.lambda = [0, 2*par.lambda];
    par.range.T      = [0, 2*par.T];
    par.range.inV    = [0, 2*par.el.a1_invV];
    par.range.scaler = [0, 2*par.el.a1_scaler];

    % ===================== hysteresis (main_hyst) =====================
    par.hyst = struct();

    % base doping scan grid (coarse)
    par.hyst.N = 101;
    par.hyst.grid_mode = "N";    % "N" or "step"
    par.hyst.step = 0.02;        % used if grid_mode=="step"
    par.hyst.psi1_0 = 0;
    par.hyst.psi2_0 = 0;
    par.hyst.T_fixed = par.T;

    % -------- refinement for COMPUTATION (grid densify) --------
    par.hyst.refine = struct();
    par.hyst.refine.enable = true;
    par.hyst.refine.ranges = [1.77, dop_end];   % densify here for solving
    par.hyst.refine.mode   = "N";         % "step" or "N"
    par.hyst.refine.step   = 0.002;          % finer step inside refine range
    par.hyst.refine.N      = 101;            % used if mode=="N"

    % -------- highlight for PLOTTING/STAT (independent) --------
    par.hyst.highlight = struct();
    par.hyst.highlight.enable = true;
    par.hyst.highlight.ranges = [1.74, 1.75];  % <<<<<< you can set any region you want

    % ===================== plot setting =====================
    par.plot.ax1_xlim = [1.72 1.80];
    par.plot.ax1_ylim = [0 15];
    
    par.plot.ax2_xlim = [1.72 1.80];
    par.plot.ax2_ylim = [0 22];

    par.plot.ax1_pos = [0.08 0.22 0.39 0.63];
    par.plot.ax2_pos = [0.54 0.22 0.39 0.63];

    % par.io.save_file = "D:\xxx\hyst_turnpoint_summary.txt";  % Windows例子
    par.io.save_file = "/Users/haoranyan/Library/CloudStorage/OneDrive-Emory/inharmonicity/BLG_project/hyst_results/hyst_data1.txt"; % Mac例子
end