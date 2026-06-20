function H = main_hyst2()
%MAIN_HYST  Hysteresis using builder interface, plot ONLY |psi1|,|psi2|.
%
% Updated:
%   1) remove shading of highlight region
%   2) allow manual xlim / ylim for each subplot
%   3) set subplot size roughly with width:height = 2:1
%   4) font size = 18
%   5) mean labels show ONLY numeric value, font 14, orange, bold
%   6) xline label at doping end shows ONLY numeric value, font 14, gray, bold

% -------------------- load params --------------------
par  = make_params();
coef = make_coeff(par);

% optional IO settings
if ~isfield(par,'io') || ~isstruct(par.io), par.io = struct(); end
if ~isfield(par.io,'save_file'), par.io.save_file = ""; end  % if empty -> uiputfile

if ~isfield(par,'lambda') || isempty(par.lambda) || ~isfinite(par.lambda)
    if isfield(par,'range') && isfield(par.range,'lambda')
        par.lambda = 0.5*sum(par.range.lambda);
    else
        par.lambda = 0;
    end
end

% -------------------- hysteresis settings (defaults -> merge par.hyst) --------------------
hyst = struct();
hyst.N = 201;
hyst.grid_mode = "N";   % "N" or "step"
hyst.step = 0.02;
hyst.psi1_0 = 0;
hyst.psi2_0 = 0;
hyst.T_fixed = par.T;

hyst.refine = struct();
hyst.refine.enable = false;
hyst.refine.ranges = zeros(0,2);
hyst.refine.mode   = "step";
hyst.refine.step   = 0.005;
hyst.refine.N      = 301;

hyst.highlight = struct();
hyst.highlight.enable = false;
hyst.highlight.ranges = zeros(0,2);

if isfield(par,'hyst') && isstruct(par.hyst)
    fns = fieldnames(hyst);
    for k=1:numel(fns)
        if isfield(par.hyst, fns{k}), hyst.(fns{k}) = par.hyst.(fns{k}); end
    end
    if isfield(par.hyst,'refine') && isstruct(par.hyst.refine)
        rf = fieldnames(hyst.refine);
        for k=1:numel(rf)
            if isfield(par.hyst.refine, rf{k}), hyst.refine.(rf{k}) = par.hyst.refine.(rf{k}); end
        end
    end
    if isfield(par.hyst,'highlight') && isstruct(par.hyst.highlight)
        hf = fieldnames(hyst.highlight);
        for k=1:numel(hf)
            if isfield(par.hyst.highlight, hf{k}), hyst.highlight.(hf{k}) = par.hyst.highlight.(hf{k}); end
        end
    end
end

% -------------------- ranges --------------------
psi1_lim = [-100,100];
psi2_lim = [-100,100];
dop_lim  = [-100,100];
if isfield(par,'range')
    if isfield(par.range,'psi1'), psi1_lim = par.range.psi1; end
    if isfield(par.range,'psi2'), psi2_lim = par.range.psi2; end
    if isfield(par.range,'dop'),  dop_lim  = par.range.dop;  end
end

% -------------------- minima-search settings --------------------
scan = struct();
scan.Npsi1 = 101;
scan.Npsi2 = 101;
scan.use_gaussian_smooth = true;
scan.smooth_sigma = 0.7;
scan.keep_topK = 40;

opt = struct();
opt.seed_jitter  = 0.12;
opt.seed_repeats = 2;
opt.max_iter     = 450;
opt.tol_fun      = 1e-10;
opt.tol_x        = 1e-8;
opt.cluster_tol  = 3e-2;
opt.fd_h         = 2e-3;
opt.min_eig_eps  = 1e-8;
opt.force_origin_and_nearby = true;
opt.nearby_delta = 0.2;

if isfield(par,'min')
    f = fieldnames(scan);
    for k=1:numel(f)
        if isfield(par.min,f{k}), scan.(f{k}) = par.min.(f{k}); end
    end
    f = fieldnames(opt);
    for k=1:numel(f)
        if isfield(par.min,f{k}), opt.(f{k}) = par.min.(f{k}); end
    end
end

% -------------------- plot options --------------------
plotopt = struct();
plotopt.fontSize = 18;      % overall axes/title/legend font
plotopt.annoFont = 14;      % mean value + xline value label
plotopt.meanColor = [0.8500 0.3250 0.0980]; % orange
plotopt.turnColor = [0.45 0.45 0.45];       % gray

% manual axis limits (empty = auto)
plotopt.ax1_xlim = [];
plotopt.ax1_ylim = [];
plotopt.ax2_xlim = [];
plotopt.ax2_ylim = [];

% manual axes positions, roughly 2:1 width:height
plotopt.ax1_pos = [0.08 0.3 0.38 0.5];
plotopt.ax2_pos = [0.54 0.3 0.38 0.5];

if isfield(par,'plot') && isstruct(par.plot)
    pf = fieldnames(plotopt);
    for k = 1:numel(pf)
        if isfield(par.plot,pf{k})
            plotopt.(pf{k}) = par.plot.(pf{k});
        end
    end
end

fs  = plotopt.fontSize;
afs = plotopt.annoFont;

% -------------------- build doping grid (with refine) --------------------
dop_base = linspace(dop_lim(1), dop_lim(2), max(61, round(hyst.N))).';
[dop_fwd, dop_bwd] = build_hyst_doping_grid_(dop_base, hyst);
N = numel(dop_fwd);

dop_turn = dop_fwd(end);

% -------------------- build lattice once --------------------
lat = coef.build_lat();

% -------------------- forward --------------------
psi_f = nan(N,2);
prev = [hyst.psi1_0, hyst.psi2_0];

for i=1:N
    dop = dop_fwd(i);
    el  = coef.build_el(hyst.T_fixed, dop);
    F = @(p1,p2) el.F(p1) + lat.V(p2) + par.lambda*p1.*p2;

    mins = find_local_minima_2d_denseSeeds_(F, psi1_lim, psi2_lim, scan, opt);
    if isempty(mins)
        pick = fallback_fminsearch_from_prev_(F, psi1_lim, psi2_lim, prev, opt);
    else
        pick = pick_closest_min_(mins, prev);
    end
    psi_f(i,:) = pick;
    prev = pick;
end

% -------------------- backward --------------------
psi_b = nan(N,2);
prev = psi_f(end,:);

for i=1:N
    dop = dop_bwd(i);
    el  = coef.build_el(hyst.T_fixed, dop);
    F = @(p1,p2) el.F(p1) + lat.V(p2) + par.lambda*p1.*p2;

    mins = find_local_minima_2d_denseSeeds_(F, psi1_lim, psi2_lim, scan, opt);
    if isempty(mins)
        pick = fallback_fminsearch_from_prev_(F, psi1_lim, psi2_lim, prev, opt);
    else
        pick = pick_closest_min_(mins, prev);
    end
    psi_b(i,:) = pick;
    prev = pick;
end

% -------------------- data for plotting --------------------
psi_f_abs = abs(psi_f);
psi_b_abs = abs(psi_b);

mask_hi_b = highlight_mask_(dop_bwd, hyst, dop_lim);
mean_b_psi1 = mean(psi_b_abs(mask_hi_b,1), 'omitnan');
mean_b_psi2 = mean(psi_b_abs(mask_hi_b,2), 'omitnan');

psi_b_turn = psi_b(1,:);
psi_b_turn_abs = abs(psi_b_turn);

% -------------------- figure --------------------
fig = figure('Color','w','Units','pixels','Position',[120 120 1200 320], ...
    'Name','Hysteresis (psi only)');

ax1 = axes('Parent',fig,'Units','normalized','Position',plotopt.ax1_pos);
plot(ax1, dop_fwd, psi_f_abs(:,1), '-','LineWidth',2); hold(ax1,'on');
plot(ax1, dop_bwd, psi_b_abs(:,1), '--','LineWidth',2);

xl1 = xline(ax1, dop_turn, '--', 'LineWidth', 2.8, ...
    'HandleVisibility','off', 'Color', plotopt.turnColor);
xl1.Label = sprintf('%.6g', dop_turn);
xl1.LabelHorizontalAlignment = 'left';
xl1.LabelVerticalAlignment   = 'middle';
xl1.FontSize = afs;
xl1.FontWeight = 'bold';
xl1.LabelOrientation = 'horizontal';

hold(ax1,'off');
grid(ax1,'off'); box(ax1,'on');
xlabel(ax1,'$\mathrm{doping}$','Interpreter','latex');
ylabel(ax1,'$|\psi_1|$','Interpreter','latex');
title(ax1, "electronic hysteresis", 'Interpreter','latex','FontWeight','normal','FontSize',fs);
legend(ax1,{'forward','backward'},'Interpreter','latex','Location','best');
set(ax1,'FontSize',fs,'TickLabelInterpreter','latex','LineWidth',1,'TickDir','out');

if ~isempty(plotopt.ax1_xlim), xlim(ax1, plotopt.ax1_xlim); end
if ~isempty(plotopt.ax1_ylim), ylim(ax1, plotopt.ax1_ylim); end

if any(mask_hi_b(:)) && isfinite(mean_b_psi1)
    text(ax1, 0.50, 0.72, sprintf('%.6g', mean_b_psi1), ...
        'Units','normalized', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontSize', afs, ...
        'FontWeight','bold', ...
        'Color', plotopt.meanColor, ...
        'Interpreter','none');
end

ax2 = axes('Parent',fig,'Units','normalized','Position',plotopt.ax2_pos);
plot(ax2, dop_fwd, psi_f_abs(:,2), '-','LineWidth',2); hold(ax2,'on');
plot(ax2, dop_bwd, psi_b_abs(:,2), '--','LineWidth',2);

xl2 = xline(ax2, dop_turn, '--', 'LineWidth', 2.8, ...
    'HandleVisibility','off', 'Color', plotopt.turnColor);
xl2.Label = sprintf('%.6g', dop_turn);
xl2.LabelHorizontalAlignment = 'left';
xl2.LabelVerticalAlignment   = 'middle';
xl2.FontSize = afs;
xl2.FontWeight = 'bold';
xl2.LabelOrientation = 'horizontal';

hold(ax2,'off');
grid(ax2,'off'); box(ax2,'on');
xlabel(ax2,'$\mathrm{doping}$','Interpreter','latex');
ylabel(ax2,'$|\psi_2|$','Interpreter','latex');
title(ax2, "lattice hysteresis", 'Interpreter','latex','FontWeight','normal','FontSize',fs);
legend(ax2,{'forward','backward'},'Interpreter','latex','Location','best');
set(ax2,'FontSize',fs,'TickLabelInterpreter','latex','LineWidth',1,'TickDir','out');

if ~isempty(plotopt.ax2_xlim), xlim(ax2, plotopt.ax2_xlim); end
if ~isempty(plotopt.ax2_ylim), ylim(ax2, plotopt.ax2_ylim); end

if any(mask_hi_b(:)) && isfinite(mean_b_psi2)
    text(ax2, 0.50, 0.72, sprintf('%.6g', mean_b_psi2), ...
        'Units','normalized', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontSize', afs, ...
        'FontWeight','bold', ...
        'Color', plotopt.meanColor, ...
        'Interpreter','none');
end

linkaxes([ax1 ax2],'x');

% -------------------- output struct --------------------
H = struct();
H.par = par;
H.coef = coef;
H.hyst = hyst;

H.dop_fwd = dop_fwd;
H.dop_bwd = dop_bwd;
H.psi_fwd = psi_f;
H.psi_bwd = psi_b;

H.dop_turn = dop_turn;
H.psi_backward_at_turn = psi_b_turn;
H.abspsi_backward_at_turn = psi_b_turn_abs;

H.highlight_mask_backward = mask_hi_b;
H.mean_backward_abs_psi1_in_highlight = mean_b_psi1;
H.mean_backward_abs_psi2_in_highlight = mean_b_psi2;

% -------------------- ask to save --------------------
msg = sprintf(['Save results?\n\n', ...
    'turning dop = %.6g\n', ...
    'backward at turn: psi1=%.6g, psi2=%.6g (abs: %.6g, %.6g)\n'], ...
    dop_turn, psi_b_turn(1), psi_b_turn(2), psi_b_turn_abs(1), psi_b_turn_abs(2));

choice = questdlg(msg, 'Save hysteresis results', 'Yes', 'No', 'No');
if strcmp(choice,'Yes')
    save_hyst_turnpoint_result_(H, par);
end

end

% ================= helpers =================
function save_hyst_turnpoint_result_(H, par)
fname = "";
if isfield(par,'io') && isfield(par.io,'save_file')
    fname = string(par.io.save_file);
end

if strlength(fname) == 0
    [f,p] = uiputfile({'*.txt','Text file (*.txt)';'*.*','All files'}, ...
        'Save hysteresis results');
    if isequal(f,0), return; end
    fname = string(fullfile(p,f));
end

dop_turn = H.dop_turn;
psi_b    = H.psi_backward_at_turn;
ab       = H.abspsi_backward_at_turn;

m1 = H.mean_backward_abs_psi1_in_highlight;
m2 = H.mean_backward_abs_psi2_in_highlight;

has_header = file_has_header_(fname, "# main_hyst results");

fid = fopen(fname,'a');
if fid < 0
    warning('Cannot open file for writing: %s', fname);
    return;
end

if ~has_header
    fprintf(fid, '# main_hyst results\n');
    fprintf(fid, '# columns:\n');
    fprintf(fid, '# dop_turn  psi1_b_turn  psi2_b_turn  abspsi1_b_turn  abspsi2_b_turn  mean_abspsi1_b_in_highlight  mean_abspsi2_b_in_highlight\n');
end

fprintf(fid, '%.16g  %.16g  %.16g  %.16g  %.16g  %.16g  %.16g\n', ...
    dop_turn, psi_b(1), psi_b(2), ab(1), ab(2), m1, m2);

fclose(fid);
fprintf('[main_hyst] saved: %s\n', fname);
end

function tf = file_has_header_(fname, header_line)
tf = false;
if ~isfile(fname), return; end

fid = fopen(fname,'r');
if fid < 0, return; end
c = onCleanup(@() fclose(fid)); %#ok<NASGU>

while true
    t = fgetl(fid);
    if ~ischar(t), return; end
    t = strtrim(string(t));
    if strlength(t)==0, continue; end
    tf = contains(t, header_line);
    return;
end
end

function [dop_fwd, dop_bwd] = build_hyst_doping_grid_(dop0, hyst)
dmin = min(dop0); dmax = max(dop0);
mode = "N"; if isfield(hyst,"grid_mode"), mode = string(hyst.grid_mode); end

if mode == "step"
    step = hyst.step;
    dop_fwd = (dmin:step:dmax).';
    if abs(dop_fwd(end)-dmax) > 1e-12, dop_fwd(end+1,1)=dmax; end %#ok<AGROW>
else
    dop_fwd = linspace(dmin, dmax, round(hyst.N)).';
end

if isfield(hyst,'refine') && isstruct(hyst.refine) && isfield(hyst.refine,'enable') && hyst.refine.enable
    rr = hyst.refine.ranges;
    if ~isempty(rr)
        rr = double(rr);
        if size(rr,2) ~= 2
            error('hyst.refine.ranges must be Mx2: [a b; ...]');
        end

        rmode = "step"; if isfield(hyst.refine,'mode'), rmode = string(hyst.refine.mode); end
        extra = zeros(0,1);

        for k=1:size(rr,1)
            a = rr(k,1); b = rr(k,2);
            if ~isfinite(a) || ~isfinite(b), continue; end
            lo = min(a,b); hi = max(a,b);
            lo = max(lo, dmin); hi = min(hi, dmax);
            if hi <= lo, continue; end

            if rmode == "N"
                nloc = hyst.refine.N;
                nloc = max(3, round(nloc));
                x = linspace(lo, hi, nloc).';
            else
                st = hyst.refine.step;
                st = max(1e-12, double(st));
                x = (lo:st:hi).';
                if abs(x(end)-hi) > 1e-12, x(end+1,1)=hi; end %#ok<AGROW>
            end
            extra = [extra; x]; %#ok<AGROW>
        end

        if ~isempty(extra)
            dop_fwd = unique([dop_fwd; extra]);
        end
    end
end

dop_bwd = flipud(dop_fwd);
end

function mask = highlight_mask_(dop, hyst, dop_lim)
mask = false(size(dop));
if ~isfield(hyst,'highlight') || ~isstruct(hyst.highlight) || ~isfield(hyst.highlight,'enable') || ~hyst.highlight.enable
    return;
end
if ~isfield(hyst.highlight,'ranges') || isempty(hyst.highlight.ranges)
    return;
end
rr = double(hyst.highlight.ranges);
if size(rr,2) ~= 2, return; end

for k=1:size(rr,1)
    lo = min(rr(k,1), rr(k,2));
    hi = max(rr(k,1), rr(k,2));
    lo = max(lo, dop_lim(1));
    hi = min(hi, dop_lim(2));
    if hi <= lo, continue; end
    mask = mask | (dop >= lo & dop <= hi);
end
end

function pick = pick_closest_min_(mins, prev)
d = sqrt(sum((mins-prev).^2,2)); [~,ix]=min(d); pick=mins(ix,:);
end

function pick = fallback_fminsearch_from_prev_(F, lim1, lim2, prev, opt)
obj = @(x) safe_F_(F, x(1), x(2), lim1, lim2);
fopt = optimset('Display','off','MaxIter',opt.max_iter,'TolFun',opt.tol_fun,'TolX',opt.tol_x);
x0 = prev;
x0(1)=min(max(x0(1),lim1(1)),lim1(2));
x0(2)=min(max(x0(2),lim2(1)),lim2(2));
try
    xhat = fminsearch(obj, x0, fopt);
catch
    xhat = x0;
end
xhat(1)=min(max(xhat(1),lim1(1)),lim1(2));
xhat(2)=min(max(xhat(2),lim2(1)),lim2(2));
pick = xhat;
end

function mins = find_local_minima_2d_denseSeeds_(F, psi1_lim, psi2_lim, scan, opt)
psi1 = linspace(psi1_lim(1), psi1_lim(2), scan.Npsi1);
psi2 = linspace(psi2_lim(1), psi2_lim(2), scan.Npsi2);
[P1,P2] = meshgrid(psi1, psi2);
Fgrid = arrayfun(@(u,v) safe_F_(F,u,v,psi1_lim,psi2_lim), P1, P2);
if scan.use_gaussian_smooth, Fgrid = gaussian_smooth_2d_(Fgrid, scan.smooth_sigma); end

mask_min = discrete_local_min_mask_(Fgrid);
[rr,cc] = find(mask_min);

if isempty(rr)
    seeds = [0,0];
else
    seeds = [P1(sub2ind(size(P1), rr, cc)), P2(sub2ind(size(P2), rr, cc))];
    E = Fgrid(sub2ind(size(Fgrid), rr, cc));
    [~,ix] = sort(E,'ascend');
    ix = ix(1:min(scan.keep_topK, numel(ix)));
    seeds = seeds(ix,:);
end

if opt.force_origin_and_nearby
    d = opt.nearby_delta;
    seeds = [seeds; 0,0; d,0; -d,0; 0,d; 0,-d];
end

obj = @(x) safe_F_(F, x(1), x(2), psi1_lim, psi2_lim);
fopt = optimset('Display','off','MaxIter',opt.max_iter,'TolFun',opt.tol_fun,'TolX',opt.tol_x);

cand = zeros(0,2);
for k=1:size(seeds,1)
    x0 = seeds(k,:);
    for rep=1:opt.seed_repeats
        xstart = x0 + opt.seed_jitter*randn(1,2);
        xstart(1)=min(max(xstart(1),psi1_lim(1)),psi1_lim(2));
        xstart(2)=min(max(xstart(2),psi2_lim(1)),psi2_lim(2));
        try
            xhat = fminsearch(obj, xstart, fopt);
        catch
            continue;
        end
        xhat(1)=min(max(xhat(1),psi1_lim(1)),psi1_lim(2));
        xhat(2)=min(max(xhat(2),psi2_lim(1)),psi2_lim(2));
        cand(end+1,:) = xhat; %#ok<AGROW>
    end
end
if isempty(cand), mins=zeros(0,2); return; end

uniq = cluster_points_(cand, opt.cluster_tol);

keep = false(size(uniq,1),1);
for i=1:size(uniq,1)
    x=uniq(i,:);
    Hh = hessian_fd_(@(u,v)F(u,v), x(1), x(2), opt.fd_h);
    ev = eig((Hh+Hh')/2);
    if all(real(ev) > opt.min_eig_eps), keep(i)=true; end
end
mins = uniq(keep,:);
end

function val = safe_F_(F, psi1, psi2, lim1, lim2)
if psi1<lim1(1) || psi1>lim1(2) || psi2<lim2(1) || psi2>lim2(2)
    val=1e30; return;
end
val = F(psi1,psi2);
if ~isfinite(val), val=1e30; end
end

function uniq = cluster_points_(P, tol)
uniq = zeros(0,2);
for i=1:size(P,1)
    x=P(i,:);
    if isempty(uniq)
        uniq=x;
    else
        d=sqrt(sum((uniq-x).^2,2));
        if all(d>tol), uniq(end+1,:)=x; end %#ok<AGROW>
    end
end
end

function Hh = hessian_fd_(f, x, y, h)
f00=f(x,y);
fxx=(f(x+h,y)-2*f00+f(x-h,y))/h^2;
fyy=(f(x,y+h)-2*f00+f(x,y-h))/h^2;
fxy=(f(x+h,y+h)-f(x+h,y-h)-f(x-h,y+h)+f(x-h,y-h))/(4*h^2);
Hh=[fxx,fxy;fxy,fyy];
end

function mask = discrete_local_min_mask_(A)
[R,C]=size(A); mask=false(R,C);
if R<3 || C<3, return; end
center=A(2:R-1,2:C-1);
n1=A(1:R-2,2:C-1); n2=A(3:R,2:C-1);
n3=A(2:R-1,1:C-2); n4=A(2:R-1,3:C);
n5=A(1:R-2,1:C-2); n6=A(1:R-2,3:C);
n7=A(3:R,1:C-2);  n8=A(3:R,3:C);
m = center<=n1 & center<=n2 & center<=n3 & center<=n4 & ...
    center<=n5 & center<=n6 & center<=n7 & center<=n8;
mask(2:R-1,2:C-1)=m;
end

function B = gaussian_smooth_2d_(A, sigma)
if sigma<=0, B=A; return; end
rad=max(1,ceil(3*sigma));
x=(-rad:rad);
g=exp(-(x.^2)/(2*sigma^2)); g=g/sum(g);
B=conv2(A,g,'same'); B=conv2(B,g','same');
end