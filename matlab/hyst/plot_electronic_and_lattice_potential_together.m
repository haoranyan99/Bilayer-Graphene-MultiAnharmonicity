function out = plot_electronic_and_lattice_potential_together()

out = struct('fig',[],'png',"");

% ==========================================================
% USER SETTINGS
% ==========================================================
dop = 1.75;
N1  = 4001;
N2  = 4001;
% ==========================================================

par = make_params();

psi1_lim = [-15 15];
psi2_lim = [0 25];
T        = par.T; 

%% ==========================================================
% electronic energy
% ==========================================================

psi1 = linspace(psi1_lim(1), psi1_lim(2), N1);

chiFun = make_chi_artificial(par, struct());
chiTD  = chiFun(T, dop);

a1 = par.el.a1_scaler * (par.el.a1_invV - abs(chiTD));
b1 = par.el.b1;

Pel = struct('a1', a1, 'b1', b1);
el  = build_electronic_energy_poly(Pel);

Fel = el.F(psi1);

% ---- find minima ----
dFel = gradient(Fel, psi1);
sgn1 = sign(dFel);
sgn1(~isfinite(sgn1)) = 0;

idx1 = find(sgn1(1:end-1) < 0 & sgn1(2:end) > 0) + 1;

psi1Min = [];
FelMin  = [];

for ii = 1:numel(idx1)

    j  = idx1(ii);
    jL = max(1, j-25);
    jR = min(numel(psi1), j+25);

    try
        xm = fminbnd(@(x) el.F(x), psi1(jL), psi1(jR));
        psi1Min(end+1) = xm; %#ok<AGROW>
        FelMin(end+1)  = el.F(xm); %#ok<AGROW>
    catch
    end

end

[psi1Min, FelMin] = dedupPairs_(psi1Min, FelMin, 1e-4);

%% ==========================================================
% lattice potential
% ==========================================================

lat = build_lattice_potential_simple(par);

psi2 = linspace(psi2_lim(1), psi2_lim(2), N2);
Vlat = lat.V(psi2);
dV   = lat.dV(psi2);

Vlat_pos = Vlat;
Vlat_pos(~(Vlat_pos > 0)) = NaN;

% ---- minima ----
sgn2 = sign(dV);
sgn2(~isfinite(sgn2)) = 0;

idx2 = find(sgn2(1:end-1) < 0 & sgn2(2:end) > 0) + 1;

psi2Min = [];
VlatMin = [];

for ii = 1:numel(idx2)

    j  = idx2(ii);
    jL = max(1, j-25);
    jR = min(numel(psi2), j+25);

    try
        xm = fminbnd(@(x) lat.V(x), psi2(jL), psi2(jR));
        vm = lat.V(xm);

        if vm > 0
            psi2Min(end+1) = xm; %#ok<AGROW>
            VlatMin(end+1) = vm; %#ok<AGROW>
        end
    catch
    end

end

[psi2Min, VlatMin] = dedupPairs_(psi2Min, VlatMin, 5e-3);

%% ==========================================================
% plot
% ==========================================================

fig = figure('Color','w','Position',[120 120 980 760]);
out.fig = fig;

%% electronic subplot
ax1 = subplot(2,1,1);

plot(ax1, psi1, Fel, 'LineWidth',2)
hold(ax1,'on')

plot(ax1, psi1Min, FelMin, 'o', ...
    'LineWidth',1.5)

for i=1:numel(psi1Min)

    text(ax1, psi1Min(i), FelMin(i), ...
        sprintf(' %.4g',FelMin(i)), ...
        'FontSize',14, ...
        'Interpreter','tex')

end

hold(ax1,'off')

grid(ax1,'on')
ax1.GridLineStyle = '--';
box(ax1,'on')

xlim(ax1,[-15 15])

xlabel(ax1,'\psi_1','Interpreter','tex','FontSize',18)
ylabel(ax1,'F_{el}(\psi_1)','Interpreter','tex','FontSize',18)

set(ax1,'FontSize',18)

%% lattice subplot
ax2 = subplot(2,1,2);

plot(ax2, psi2, Vlat_pos, 'LineWidth',2)
hold(ax2,'on')

plot(ax2, psi2Min, VlatMin, 'o', ...
    'LineWidth',1.5)

for i=1:numel(psi2Min)

    text(ax2, psi2Min(i), VlatMin(i), ...
        sprintf(' %.4g',VlatMin(i)), ...
        'FontSize',14, ...
        'Interpreter','tex')

end

hold(ax2,'off')

grid(ax2,'on')
ax2.GridLineStyle = '--';
box(ax2,'on')

xlim(ax2,[0 25])

xlabel(ax2,'\psi_2','Interpreter','tex','FontSize',18)
ylabel(ax2,'V(\psi_2)','Interpreter','tex','FontSize',18)

set(ax2,'FontSize',18)

end


% ================= helper =================

function [x2,y2] = dedupPairs_(x,y,tol)

if isempty(x)
    x2=[]; y2=[];
    return
end

[x,ord]=sort(x);
y=y(ord);

keep=true(size(x));

for i=2:numel(x)
    if abs(x(i)-x(i-1))<tol
        keep(i)=false;
    end
end

x2=x(keep);
y2=y(keep);

end