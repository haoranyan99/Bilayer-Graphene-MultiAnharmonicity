function plot_electronic_energy_poly_current()

% ==========================================================
% USER SETTINGS
% ==========================================================

dop = 1.7;     % <<<<<< 在这里改 doping

% ==========================================================
% load parameters
% ==========================================================

par = make_params();

psi1_lim = par.range.psi1;
T        = par.T;

N = 4001;
psi = linspace(psi1_lim(1), psi1_lim(2), N);

% ==========================================================
% compute chi(T,dop)
% ==========================================================

chiFun = make_chi_artificial(par, struct());
chiTD  = chiFun(T, dop);

% ==========================================================
% electronic energy coefficients
% ==========================================================

a1 = par.el.a1_scaler * (par.el.a1_invV - abs(chiTD));
b1 = par.el.b1;

Pel = struct('a1', a1, 'b1', b1);

el = build_electronic_energy_poly(Pel);

F = el.F(psi);

% ==========================================================
% find minima
% ==========================================================

dF = gradient(F, psi);

sgn = sign(dF);
sgn(~isfinite(sgn)) = 0;

idx = find(sgn(1:end-1)<0 & sgn(2:end)>0) + 1;

psiMin = [];
FMin   = [];

for ii = 1:numel(idx)

    j  = idx(ii);
    jL = max(1, j-25);
    jR = min(numel(psi), j+25);

    a = psi(jL);
    b = psi(jR);

    try
        xm = fminbnd(@(x) el.F(x), a, b);

        psiMin(end+1) = xm;
        FMin(end+1)   = el.F(xm);

    catch
    end

end

% remove duplicates
[psiMin,ord] = sort(psiMin);
FMin = FMin(ord);

tol = 5e-6*(psi1_lim(2)-psi1_lim(1));
keep = true(size(psiMin));

for i=2:numel(psiMin)
    if abs(psiMin(i)-psiMin(i-1)) < tol
        keep(i) = false;
    end
end

psiMin = psiMin(keep);
FMin   = FMin(keep);

% ==========================================================
% plot
% ==========================================================

figure('Color','w','Position',[150 150 900 520])

plot(psi, F, 'LineWidth',2)
hold on

plot(psiMin, FMin,'ro','MarkerFaceColor','r','LineWidth',1.5)

for i = 1:numel(psiMin)

    text(psiMin(i),FMin(i), ...
        sprintf('  min @ %.4g',psiMin(i)), ...
        'FontSize',10)

end

hold off

grid on
box on

xlim(psi1_lim)

xlabel('\psi_1','Interpreter','latex')
ylabel('$F_{el}(\psi_1)$','Interpreter','latex')

title(sprintf('Electronic Energy | T=%g  dop=%g  chi=%g  a1=%g  b1=%g', ...
    T, dop, chiTD, a1, b1),'Interpreter','none')

end