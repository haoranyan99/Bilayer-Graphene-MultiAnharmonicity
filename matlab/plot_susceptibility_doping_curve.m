function plot_susceptibility_doping_curve(files, varargin)
%PLOT_SUSCEPTIBILITY_DOPING_CURVE Plot chi(q) versus doping.
%
% Usage:
%   plot_susceptibility_doping_curve
%   plot_susceptibility_doping_curve("data/chi_doping_curve_patch_...csv")
%   plot_susceptibility_doping_curve([], "q_index", [0 0])
%   plot_susceptibility_doping_curve(["data/chi_mu0.csv", "data/chi_mu0p01.csv"])
%   plot_susceptibility_doping_curve(..., "quantity", "abs", "save", true)
%
% Name-value options:
%   q_index  : integer q shift [iq jq], default [0 0]
%   q_Ainv   : physical q [qx qy] in 1/A; if set, overrides q_index
%   quantity : "real" (default), "imag", or "abs"
%   save     : false (default) or true
%   out_file : output image path when save=true

if nargin < 1 || isempty(files)
    files = default_or_selected_chi_files();
end
files = normalize_files(files);

p = inputParser;
addParameter(p, "q_index", [0 0]);
addParameter(p, "q_Ainv", []);
addParameter(p, "quantity", "real");
addParameter(p, "save", false);
addParameter(p, "out_file", "");
parse(p, varargin{:});

q_index = double(p.Results.q_index);
q_Ainv = double(p.Results.q_Ainv);
quantity = lower(string(p.Results.quantity));
do_save = logical(p.Results.save);
out_file = string(p.Results.out_file);

first_data = readmatrix(char(files(1)), "FileType", "text", "CommentStyle", "#");
if is_doping_curve_file(char(files(1)), first_data)
    [mu, chi_value, selected_q, y_label, suffix] = read_curve_file(char(files(1)), quantity);
else
    mu = nan(numel(files), 1);
    chi_value = nan(numel(files), 1);
    selected_q = nan(numel(files), 2);

    for i = 1:numel(files)
        csv_file = char(files(i));
        data = readmatrix(csv_file, "FileType", "text", "CommentStyle", "#");
        if size(data, 2) < 8
            error("Expected at least 8 columns in %s.", csv_file);
        end

        mu(i) = read_mu_from_header(csv_file);
        row_idx = select_q_row(data, q_index, q_Ainv);
        selected_q(i, :) = data(row_idx, 3:4);

        switch quantity
            case {"real", "re", "chi_re"}
                chi_value(i) = data(row_idx, 5);
                y_label = "Re \chi(q)";
                suffix = "re";
            case {"imag", "im", "chi_im"}
                chi_value(i) = data(row_idx, 6);
                y_label = "Im \chi(q)";
                suffix = "im";
            case {"abs", "magnitude", "mag"}
                chi_value(i) = hypot(data(row_idx, 5), data(row_idx, 6));
                y_label = "|\chi(q)|";
                suffix = "abs";
            otherwise
                error("Unknown quantity '%s'. Use real, imag, or abs.", quantity);
        end
    end
end

[mu, order] = sort(mu);
chi_value = chi_value(order);
selected_q = selected_q(order, :);

figure("Color", "w");
plot(mu, chi_value, "-o", "LineWidth", 1.5, "MarkerSize", 5);
box on;
grid on;
xlabel("\mu (eV)");
ylabel(y_label);
title(sprintf("%s vs doping at q=(%.6g, %.6g) 1/A", y_label, selected_q(1, 1), selected_q(1, 2)));

if do_save
    if strlength(out_file) == 0
        root_dir = fileparts(fileparts(mfilename("fullpath")));
        out_file = fullfile(root_dir, "data", char("chi_doping_curve_" + suffix + ".png"));
    end
    exportgraphics(gcf, out_file, "Resolution", 220);
    fprintf("Wrote: %s\n", out_file);
end
end

function files = default_or_selected_chi_files()
root_dir = fileparts(fileparts(mfilename("fullpath")));
default_data_dir = fullfile(root_dir, "data");
curve_files = dir(fullfile(default_data_dir, "chi_doping_curve_*.csv"));
if ~isempty(curve_files)
    [~, idx] = max([curve_files.datenum]);
    files = string(fullfile(curve_files(idx).folder, curve_files(idx).name));
    return;
end

found = dir(fullfile(default_data_dir, "chi_*.csv"));
if ~isempty(found)
    files = strings(numel(found), 1);
    for i = 1:numel(found)
        files(i) = fullfile(found(i).folder, found(i).name);
    end
    return;
end

[file, folder] = uigetfile( ...
    {"*.csv", "Susceptibility CSV files (*.csv)"; "*.*", "All files (*.*)"}, ...
    "Select susceptibility CSV files for doping curve", ...
    fullfile(default_data_dir, "*.csv"), ...
    "MultiSelect", "on");
if isequal(file, 0)
    error("No susceptibility CSV files selected.");
end
if iscell(file)
    files = strings(numel(file), 1);
    for i = 1:numel(file)
        files(i) = fullfile(folder, file{i});
    end
else
    files = string(fullfile(folder, file));
end
end

function tf = is_doping_curve_file(csv_file, data)
tf = size(data, 2) >= 9;
fid = fopen(csv_file, "r");
if fid < 0
    return;
end
cleanup = onCleanup(@() fclose(fid));
line = fgetl(fid);
if ischar(line)
    tf = tf && contains(line, "blg_susceptibility_doping_curve_csv_v1");
end
end

function [mu, chi_value, selected_q, y_label, suffix] = read_curve_file(csv_file, quantity)
data = readmatrix(csv_file, "FileType", "text", "CommentStyle", "#");
if size(data, 2) < 9
    error("Expected doping curve columns: mu_eV,iq,jq,qx_Ainv,qy_Ainv,chi_re,chi_im,nKpair,nK.");
end

mu = data(:, 1);
selected_q = data(:, 4:5);
switch quantity
    case {"real", "re", "chi_re"}
        chi_value = data(:, 6);
        y_label = "Re \chi(q)";
        suffix = "re";
    case {"imag", "im", "chi_im"}
        chi_value = data(:, 7);
        y_label = "Im \chi(q)";
        suffix = "im";
    case {"abs", "magnitude", "mag"}
        chi_value = hypot(data(:, 6), data(:, 7));
        y_label = "|\chi(q)|";
        suffix = "abs";
    otherwise
        error("Unknown quantity '%s'. Use real, imag, or abs.", quantity);
end
end

function files = normalize_files(files)
if ischar(files)
    files = string({files});
elseif iscell(files)
    files = string(files);
else
    files = string(files);
end
files = files(:);
end

function row_idx = select_q_row(data, q_index, q_Ainv)
if ~isempty(q_Ainv)
    dq = hypot(data(:, 3) - q_Ainv(1), data(:, 4) - q_Ainv(2));
    [best, row_idx] = min(dq);
    if best > 1e-9
        warning("Using nearest q point: requested (%.6g, %.6g), found (%.6g, %.6g).", ...
            q_Ainv(1), q_Ainv(2), data(row_idx, 3), data(row_idx, 4));
    end
    return;
end

mask = (round(data(:, 1)) == q_index(1)) & (round(data(:, 2)) == q_index(2));
row_idx = find(mask, 1);
if isempty(row_idx)
    error("q_index [%d %d] not found in susceptibility file.", q_index(1), q_index(2));
end
end

function mu = read_mu_from_header(csv_file)
mu = NaN;
fid = fopen(csv_file, "r");
if fid < 0
    error("Cannot open %s.", csv_file);
end
cleanup = onCleanup(@() fclose(fid));

while true
    line = fgetl(fid);
    if ~ischar(line)
        break;
    end
    if startsWith(line, "# patch_metadata = ")
        json_txt = extractAfter(string(line), "# patch_metadata = ");
        meta = jsondecode(json_txt);
        if isfield(meta, "mu_eV")
            mu = double(meta.mu_eV);
            return;
        end
        if isfield(meta, "reference_mu_eV")
            mu = double(meta.reference_mu_eV);
            return;
        end
    end
    if ~startsWith(line, "#")
        break;
    end
end

[~, name, ~] = fileparts(csv_file);
token = regexp(name, "mu([mp0-9]+)", "tokens", "once");
if ~isempty(token)
    mu = str2double(strrep(strrep(token{1}, "m", "-"), "p", "."));
end
if isnan(mu)
    error("Cannot read mu_eV from CSV metadata or filename: %s.", csv_file);
end
end
