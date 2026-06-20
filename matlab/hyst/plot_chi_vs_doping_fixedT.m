function out = plot_chi_vs_doping_fixedT()

out = struct('fig',[],'png',"",'T_list',[],'dop_lim',[]);

% ================= USER SETTINGS =================
T_list = [6];
N = 1500;

doSave   = false;
saveName = "chi_a1_vs_doping_fixedT.png";
saveDPI  = 200;
% =================================================

par = make_params();

if ~isfield(par,'range') || ~isfield(par.range,'dop') || numel(par.range.dop)~=2
    par.range.dop = [1.6, 1.95];
end

dop_lim = par.range.dop;
dop = linspace(dop_lim(1), dop_lim(2), N);

out.T_list = T_list;
out.dop_lim = dop_lim;

% ---------------- chi / a1 ----------------
chi = make_chi_artificial(par);

a1 = @(T,dop) par.el.a1_scaler .* ...
    (par.el.a1_invV - abs(chi(T,dop)));

% ---------------- figure ----------------
fig = figure('Color','w','Position',[150 150 900 520]);
out.fig = fig;

ax = axes(fig);
hold(ax,'on')

% ================= chi (left axis) =================
yyaxis left

for i = 1:numel(T_list)

    T = T_list(i);
    y = chi(T,dop);

    plot(ax,dop,y,'LineWidth',2);

end

ylabel('\chi(T,dop)','Interpreter','latex')

% ================= a1 (right axis) =================
yyaxis right

for i = 1:numel(T_list)

    T = T_list(i);
    y = a1(T,dop);

    plot(ax,dop,y,'--','LineWidth',2);

end

ylabel('a_1(T,dop)','Interpreter','latex')

% ================= formatting =================

xlabel(ax,'doping','Interpreter','latex')

xlim(ax,dop_lim)

grid(ax,'on')
box(ax,'on')

% 所有轴颜色设为黑色
ax.XColor = 'k';
ax.YColor = 'k';

yyaxis left
ax.YColor = 'k';

yyaxis right
ax.YColor = 'k';

title('\chi and a_1 vs doping','Interpreter','latex')

end