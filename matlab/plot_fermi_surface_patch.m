function plot_fermi_surface_patch(patch_file, varargin)
%PLOT_FERMI_SURFACE_PATCH Plot Fermi occupation and contour from a BLG patch .npz file.
%
% Usage:
%   plot_fermi_surface_patch
%   plot_fermi_surface_patch("data/patch_tb_gap0p012_K_cart_Nx41_Ny41_dx0p002_dy0p002_T0.npz")
%   plot_fermi_surface_patch(..., "mu_eV", 0.01)
%   plot_fermi_surface_patch(..., "band", 2, "xlim_Ainv", [-0.05 0.05], "save", true)
%
% Name-value options:
%   mu_eV          : target chemical potential in eV; [] opens a selector
%   band           : one-based band index to draw, default 2
%   source         : "energy" (default) for E_n(k)-mu=0 or "occ" for f_n(k)=0.5
%   grid_N         : interpolation grid size for contour, default 301
%   xlim_Ainv      : x limits after subtracting K, default []
%   ylim_Ainv      : y limits after subtracting K, default []
%   half_width_Ainv: used when xlim/ylim are empty, default 0.1
%   save           : false (default) or true
%   out_file       : output image path when save=true

if nargin < 1 || strlength(string(patch_file)) == 0
    patch_file = latest_patch_file();
end
patch_file = char(patch_file);

p = inputParser;
addParameter(p, "mu_eV", []);
addParameter(p, "band", 2);
addParameter(p, "source", "energy");
addParameter(p, "grid_N", 301);
addParameter(p, "xlim_Ainv", []);
addParameter(p, "ylim_Ainv", []);
addParameter(p, "half_width_Ainv", 0.1);
addParameter(p, "save", false);
addParameter(p, "out_file", "");
parse(p, varargin{:});

mu_eV = double(p.Results.mu_eV);
band = double(p.Results.band);
contour_source = lower(string(p.Results.source));
grid_N = double(p.Results.grid_N);
xlim_Ainv = double(p.Results.xlim_Ainv);
ylim_Ainv = double(p.Results.ylim_Ainv);
half_width_Ainv = double(p.Results.half_width_Ainv);
do_save = logical(p.Results.save);
out_file = string(p.Results.out_file);

tmp_dir = tempname;
mkdir(tmp_dir);
cleanup = onCleanup(@() cleanup_tmp(tmp_dir));
unzip(patch_file, tmp_dir);

kx = read_npy_numeric(fullfile(tmp_dir, "kx.npy"));
ky = read_npy_numeric(fullfile(tmp_dir, "ky.npy"));
evals = read_npy_numeric(fullfile(tmp_dir, "evals.npy"));
iq = read_npy_numeric(fullfile(tmp_dir, "iq.npy"));
jq = read_npy_numeric(fullfile(tmp_dir, "jq.npy"));
center_k = read_npy_numeric(fullfile(tmp_dir, "center_k.npy"));
mu_eV = choose_mu_eV(tmp_dir, mu_eV);
occ = read_occ_for_mu(tmp_dir, mu_eV);

if band < 1 || band > size(evals, 2)
    error("Band index %d out of range 1..%d.", band, size(evals, 2));
end

kx_rel = kx - center_k(1);
ky_rel = ky - center_k(2);
[iq_unique, jq_unique, kx_grid, ky_grid] = patch_grid(iq, jq, kx_rel, ky_rel);
energy = evals(:, band);
occ_band = occ(:, band);

[fig_occ, fig_contour] = plot_occ_and_contour( ...
    kx_rel, ky_rel, iq, jq, iq_unique, jq_unique, kx_grid, ky_grid, ...
    energy, occ_band, mu_eV, band, contour_source, grid_N, patch_file, ...
    xlim_Ainv, ylim_Ainv, half_width_Ainv);

if do_save
    [occ_file, contour_file] = output_file_pair(patch_file, out_file, mu_eV, band);
    figure(fig_occ);
    occ_file = save_current_figure(occ_file);
    figure(fig_contour);
    contour_file = save_current_figure(contour_file);
    fprintf("Wrote: %s\n", occ_file);
    fprintf("Wrote: %s\n", contour_file);
end
end

function [fig_occ, fig_contour] = plot_occ_and_contour( ...
    kx, ky, iq, jq, iq_unique, jq_unique, kx_grid, ky_grid, ...
    energy, occ_band, mu_eV, band, contour_source, grid_N, patch_file, ...
    xlim_Ainv, ylim_Ainv, half_width_Ainv)

fig_occ = figure("Color", "w");
scatter(kx, ky, 18, occ_band, "filled");
axis equal tight;
apply_k_limits(gca, xlim_Ainv, ylim_Ainv, half_width_Ainv);
box on;
grid on;
colormap(parula);
cb = colorbar;
cb.Label.String = sprintf("f_{%d}(k)", band);
xlabel("k_x-K_x (1/A)");
ylabel("k_y-K_y (1/A)");
title(sprintf("Occupation scatter, band %d, \\mu = %.6g eV", band, mu_eV));
set(gca, "FontSize", 11, "LineWidth", 1.1, "TickDir", "in");

fig_contour = figure("Color", "w");
hold on;
switch contour_source
    case "energy"
        z = energy - mu_eV;
        contour_level = 0.0;
        contour_label = sprintf("E_%d(k)-\\mu=0", band);
    case "occ"
        z = occ_band;
        contour_level = 0.5;
        contour_label = sprintf("f_%d(k)=0.5", band);
    otherwise
        error("Unknown source '%s'. Use energy or occ.", contour_source);
end

ok = isfinite(kx) & isfinite(ky) & isfinite(z);
if nnz(ok) < 3
    error("Too few finite k points for contour.");
end

xg = linspace(min(kx(ok)), max(kx(ok)), grid_N);
yg = linspace(min(ky(ok)), max(ky(ok)), grid_N);
[Xg, Yg] = meshgrid(xg, yg);
interp_z = scatteredInterpolant(kx(ok), ky(ok), z(ok), "linear", "none");
Zg = interp_z(Xg, Yg);

[C, h] = contour(Xg, Yg, Zg, [contour_level contour_level], ...
    "LineWidth", 2.0, "LineColor", [0.1 0.2 0.8]);
if isempty(C) || isempty(h)
    warning("No Fermi contour found for band %d at mu=%.6g eV. Try a mu inside this band's energy range.", band, mu_eV);
end

% Lightly show sampled patch points so an empty contour still has context.
plot(kx, ky, ".", "Color", [0.75 0.75 0.75], "MarkerSize", 4);
axis equal tight;
apply_k_limits(gca, xlim_Ainv, ylim_Ainv, half_width_Ainv);
box on;
grid on;
xlabel("k_x-K_x (1/A)");
ylabel("k_y-K_y (1/A)");
title(sprintf("Fermi contour, %s, %s", contour_label, strip_path(patch_file)), ...
    "Interpreter", "none");
set(gca, "FontSize", 11, "LineWidth", 1.1, "TickDir", "in");

% Keep direct-grid variables touched for debugging shape mismatches.
if isempty(values_on_grid(iq, jq, energy, iq_unique, jq_unique)) || isempty(kx_grid) || isempty(ky_grid)
    error("Internal grid construction failed.");
end
end

function occ = read_occ_for_mu(tmp_dir, mu_eV)
evals = read_npy_numeric(fullfile(tmp_dir, "evals.npy"));
occ_file = fullfile(tmp_dir, "occ.npy");
occ = read_npy_numeric(occ_file);

mu_values_file = fullfile(tmp_dir, "mu_values.npy");
occ_mu_file = fullfile(tmp_dir, "occ_mu.npy");
if exist(mu_values_file, "file") && exist(occ_mu_file, "file")
    mu_values = read_npy_numeric(mu_values_file);
    occ_mu = read_npy_numeric(occ_mu_file);
    [best, idx] = min(abs(mu_values(:) - mu_eV));
    if best < 1e-12
        occ = squeeze(occ_mu(idx, :, :));
        return;
    end
end

occ = double(evals < mu_eV);
occ(abs(evals - mu_eV) < 1e-12) = 0.5;
end

function mu_eV = choose_mu_eV(tmp_dir, mu_eV)
if ~isempty(mu_eV)
    mu_eV = double(mu_eV(1));
    return;
end

mu_values_file = fullfile(tmp_dir, "mu_values.npy");
if ~exist(mu_values_file, "file")
    answer = inputdlg({"mu_eV"}, "Select chemical potential", 1, {"0"});
    if isempty(answer)
        error("No mu_eV selected.");
    end
    mu_eV = str2double(answer{1});
    if ~isfinite(mu_eV)
        error("Invalid mu_eV.");
    end
    return;
end

mu_values = read_npy_numeric(mu_values_file);
labels = arrayfun(@(x) sprintf("%.8g eV", x), mu_values(:), "UniformOutput", false);
if usejava("desktop")
    [idx, ok] = listdlg( ...
        "PromptString", "Select mu_eV", ...
        "SelectionMode", "single", ...
        "ListString", labels, ...
        "InitialValue", max(1, ceil(numel(labels) / 2)), ...
        "ListSize", [220 360]);
    if ~ok
        error("No mu_eV selected.");
    end
    mu_eV = double(mu_values(idx));
else
    [~, idx] = min(abs(mu_values(:)));
    mu_eV = double(mu_values(idx));
    fprintf("No GUI selector available. Using nearest mu to 0: %.8g eV\n", mu_eV);
end
end

function apply_k_limits(ax, xlim_Ainv, ylim_Ainv, half_width_Ainv)
if isempty(xlim_Ainv)
    xlim_Ainv = half_width_Ainv * [-1 1];
end
if isempty(ylim_Ainv)
    ylim_Ainv = half_width_Ainv * [-1 1];
end
if numel(xlim_Ainv) ~= 2 || numel(ylim_Ainv) ~= 2
    error("xlim_Ainv and ylim_Ainv must each have two values.");
end
xlim(ax, xlim_Ainv);
ylim(ax, ylim_Ainv);
end

function patch_file = latest_patch_file()
root_dir = fileparts(fileparts(mfilename("fullpath")));
default_data_dir = fullfile(root_dir, "data");
files = dir(fullfile(default_data_dir, "patch_*.npz"));
if isempty(files)
    [file, folder] = uigetfile( ...
        {"*.npz", "Patch NPZ files (*.npz)"; "*.*", "All files (*.*)"}, ...
        "Select BLG patch file", ...
        fullfile(default_data_dir, "*.npz"));
    if isequal(file, 0)
        error("No patch file selected.");
    end
    patch_file = fullfile(folder, file);
    return;
end
[~, idx] = max([files.datenum]);
patch_file = fullfile(files(idx).folder, files(idx).name);
end

function cleanup_tmp(tmp_dir)
if exist(tmp_dir, "dir")
    rmdir(tmp_dir, "s");
end
end

function [iq_unique, jq_unique, kx_grid, ky_grid] = patch_grid(iq, jq, kx, ky)
iq_unique = unique(iq(:).');
jq_unique = unique(jq(:));
kx_grid = values_on_grid(iq, jq, kx, iq_unique, jq_unique);
ky_grid = values_on_grid(iq, jq, ky, iq_unique, jq_unique);
end

function grid = values_on_grid(iq, jq, values, iq_unique, jq_unique)
grid = nan(numel(jq_unique), numel(iq_unique));
for i = 1:numel(values)
    ix = find(iq_unique == iq(i), 1);
    iy = find(jq_unique == jq(i), 1);
    if isempty(ix) || isempty(iy)
        error("Failed to map patch point to grid.");
    end
    grid(iy, ix) = values(i);
end
end

function [occ_file, contour_file] = output_file_pair(patch_file, out_file, mu_eV, band)
if strlength(out_file) > 0
    out_file = convertStringsToChars(out_file);
    [folder, name, ext] = fileparts(out_file);
    if isempty(ext)
        ext = ".png";
    end
    occ_file = fullfile(folder, [name '_occ' ext]);
    contour_file = fullfile(folder, [name '_contour' ext]);
    return;
end

[folder, name, ~] = fileparts(patch_file);
mu_tag = strrep(strrep(sprintf("%.10g", mu_eV), "-", "m"), ".", "p");
occ_file = fullfile(folder, sprintf("%s_band%d_mu%s_occ.png", name, band, mu_tag));
contour_file = fullfile(folder, sprintf("%s_band%d_mu%s_contour.png", name, band, mu_tag));
end

function arr = read_npy_numeric(path_in)
fid = fopen(path_in, "rb");
if fid < 0
    error("Cannot open %s.", path_in);
end
cleanup = onCleanup(@() fclose(fid));

magic = fread(fid, 6, "uint8=>uint8").';
if numel(magic) ~= 6 || magic(1) ~= 147 || ~strcmp(char(magic(2:end)), "NUMPY")
    error("Not a NumPy .npy file: %s.", path_in);
end

major = fread(fid, 1, "uint8");
minor = fread(fid, 1, "uint8");
if major == 1
    header_len = fread(fid, 1, "uint16", 0, "ieee-le");
elseif major == 2 || major == 3
    header_len = fread(fid, 1, "uint32", 0, "ieee-le");
else
    error("Unsupported .npy version %d.%d.", major, minor);
end

header = char(fread(fid, double(header_len), "uint8=>uint8").');
descr = regexp(header, "'descr':\s*'([^']+)'", "tokens", "once");
fortran_order = regexp(header, "'fortran_order':\s*(False|True)", "tokens", "once");
shape_token = regexp(header, "'shape':\s*\(([^)]*)\)", "tokens", "once");
if isempty(descr) || isempty(fortran_order) || isempty(shape_token)
    error("Could not parse .npy header in %s.", path_in);
end
if strcmp(fortran_order{1}, "True")
    error("Fortran-order .npy arrays are not supported: %s.", path_in);
end

shape_txt = strtrim(shape_token{1});
shape_parts = regexp(shape_txt, ",", "split");
shape = [];
for i = 1:numel(shape_parts)
    val = strtrim(shape_parts{i});
    if strlength(string(val)) > 0
        shape(end + 1) = str2double(val); %#ok<AGROW>
    end
end
if isempty(shape)
    shape = 1;
end

dtype = descr{1};
switch dtype
    case {"<f8", "|f8"}
        precision = "double=>double";
    case {"<f4", "|f4"}
        precision = "single=>double";
    case {"<i8", "|i8"}
        precision = "int64=>double";
    case {"<i4", "|i4"}
        precision = "int32=>double";
    case {"<u8", "|u8"}
        precision = "uint64=>double";
    case {"<u4", "|u4"}
        precision = "uint32=>double";
    otherwise
        error("Unsupported dtype %s in %s.", dtype, path_in);
end

raw = fread(fid, prod(shape), precision, 0, "ieee-le");
if numel(raw) ~= prod(shape)
    error("Unexpected end of file in %s.", path_in);
end

if numel(shape) == 1
    arr = reshape(raw, shape(1), 1);
else
    arr = reshape(raw, fliplr(shape));
    arr = permute(arr, numel(shape):-1:1);
end
end

function s = strip_path(path_in)
[~, name, ext] = fileparts(path_in);
s = [name ext];
end

function out_file = save_current_figure(out_file)
out_file = convertStringsToChars(out_file);
if iscell(out_file)
    out_file = out_file{1};
end
out_file = char(out_file);
out_file = out_file(:).';
[folder, ~, ~] = fileparts(out_file);
if strlength(string(folder)) > 0 && ~exist(folder, "dir")
    mkdir(folder);
end
try
    exportgraphics(gcf, out_file, 'Resolution', 220);
catch
    saveas(gcf, out_file);
end
end
