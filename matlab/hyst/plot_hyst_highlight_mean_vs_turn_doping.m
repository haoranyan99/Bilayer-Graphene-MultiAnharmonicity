function out = plot_hyst_highlight_mean_vs_turn_doping()
% plot_hyst_highlight_mean_vs_turn_doping
% Read saved hyst results and plot:
%   mean(|psi1|)_backward_in_highlight vs dop_turn
%   mean(|psi2|)_backward_in_highlight vs dop_turn
%
% Expected columns:
% dop_turn  psi1_b_turn  psi2_b_turn  abspsi1_b_turn  abspsi2_b_turn  mean_abspsi1_b_in_highlight  mean_abspsi2_b_in_highlight

out = struct();
out.file = "";
out.data = [];

% ---------------- USER: set file ----------------
fname = "hyst_results/hyst_data1.txt";   % <- change if needed
% fname = "/Users/xxx/hyst_results.txt";
% fname = "D:\xxx\hyst_results.txt";
% -----------------------------------------------

if ~isfile(fname)
    error('File not found: %s', fname);
end
out.file = fname;

% ---------- read numeric rows, skipping # ----------
fid = fopen(fname,'r');
if fid < 0, error('Cannot open file: %s', fname); end
c = onCleanup(@() fclose(fid));

rows = zeros(0,7);
while true
    line = fgetl(fid);
    if ~ischar(line), break; end
    s = strtrim(string(line));
    if strlength(s)==0, continue; end
    if startsWith(s,"#"), continue; end

    v = sscanf(char(s), '%f');
    if numel(v) < 7
        % if older files have fewer cols, pad with NaN
        vv = nan(7,1);
        vv(1:min(numel(v),7)) = v(1:min(numel(v),7));
        v = vv;
    else
        v = v(1:7);
    end
    rows(end+1,:) = v(:).'; %#ok<AGROW>
end

if isempty(rows)
    error('No numeric rows found in: %s', fname);
end

dop_turn = rows(:,1);
mpsi1    = rows(:,6);
mpsi2    = rows(:,7);

% keep finite
good = isfinite(dop_turn) & isfinite(mpsi1) & isfinite(mpsi2);
dop_turn = dop_turn(good);
mpsi1    = mpsi1(good);
mpsi2    = mpsi2(good);

% sort by dop_turn
[dop_turn, ix] = sort(dop_turn);
mpsi1 = mpsi1(ix);
mpsi2 = mpsi2(ix);

out.data = [dop_turn, mpsi1, mpsi2];

% ---------------- plot ----------------
fs = 13;

fig1 = figure('Color','w','Name','Backward mean(|psi1|) in highlight vs dop_turn');
ax1 = axes(fig1); %#ok<LAXES>
plot(ax1, dop_turn, mpsi1, '-o', 'LineWidth', 2, 'MarkerSize', 5);
box(ax1,'on'); grid(ax1,'off');
xlabel(ax1,'dop\_turn (scan turning point)','Interpreter','none');
ylabel(ax1,'mean\_abspsi1\_b\_in\_highlight','Interpreter','none');
title(ax1,'Electronic: backward mean(|\psi_1|) in highlight vs dop\_turn', ...
    'Interpreter','tex','FontWeight','normal');
set(ax1,'FontSize',fs,'LineWidth',1,'TickDir','out');

fig2 = figure('Color','w','Name','Backward mean(|psi2|) in highlight vs dop_turn');
ax2 = axes(fig2); %#ok<LAXES>
plot(ax2, dop_turn, mpsi2, '-o', 'LineWidth', 2, 'MarkerSize', 5);
box(ax2,'on'); grid(ax2,'off');
xlabel(ax2,'dop\_turn (scan turning point)','Interpreter','none');
ylabel(ax2,'mean\_abspsi2\_b\_in\_highlight','Interpreter','none');
title(ax2,'Lattice: backward mean(|\psi_2|) in highlight vs dop\_turn', ...
    'Interpreter','tex','FontWeight','normal');
set(ax2,'FontSize',fs,'LineWidth',1,'TickDir','out');

out.fig_meanpsi1 = fig1;
out.fig_meanpsi2 = fig2;

end