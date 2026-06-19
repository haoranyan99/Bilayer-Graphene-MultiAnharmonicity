function plot_susceptibility_qdep(csv_file, varargin)
%PLOT_SUSCEPTIBILITY_QDEP Plot q-dependent BLG susceptibility.
%
% Usage:
%   plot_susceptibility_qdep
%   plot_susceptibility_qdep("data/chi_patch_tb_gap0p012_K_cart_Nh20_dk0p002_mu0_T0_eta0p0001_qmesh_Nh4.csv")
%   plot_susceptibility_qdep(..., "quantity", "imag")
%   plot_susceptibility_qdep(..., "quantity", "abs", "save", true)
%
% Name-value options:
%   quantity : "real" (default), "imag", or "abs"
%   save     : false (default) or true
%   out_file : output image path when save=true

if nargin < 1 || strlength(string(csv_file)) == 0
    csv_file = latest_chi_file();
end
csv_file = char(csv_file);

p = inputParser;
addParameter(p, "quantity", "real");
addParameter(p, "save", false);
addParameter(p, "out_file", "");
parse(p, varargin{:});

quantity = lower(string(p.Results.quantity));
do_save = logical(p.Results.save);
out_file = string(p.Results.out_file);

data = readmatrix(csv_file, "FileType", "text", "CommentStyle", "#");
if size(data, 2) < 8
    error("Expected at least 8 columns: iq,jq,qx_Ainv,qy_Ainv,chi_re,chi_im,nKpair,nK.");
end

qx = data(:, 3);
qy = data(:, 4);
switch quantity
    case {"real", "re", "chi_re"}
        z = data(:, 5);
        zlabel = "Re \chi(q)";
        suffix = "re";
    case {"imag", "im", "chi_im"}
        z = data(:, 6);
        zlabel = "Im \chi(q)";
        suffix = "im";
    case {"abs", "magnitude", "mag"}
        z = hypot(data(:, 5), data(:, 6));
        zlabel = "|\chi(q)|";
        suffix = "abs";
    otherwise
        error("Unknown quantity '%s'. Use real, imag, or abs.", quantity);
end

[qx_grid, qy_grid, z_grid] = regular_grid(qx, qy, z);

figure("Color", "w");
imagesc(qx_grid(1, :), qy_grid(:, 1), z_grid);
set(gca, "YDir", "normal");
axis image;
box on;
colormap(parula);
cb = colorbar;
cb.Label.String = zlabel;
xlabel("q_x (1/A)");
ylabel("q_y (1/A)");
title(sprintf("%s: %s", zlabel, strip_path(csv_file)), "Interpreter", "none");

if do_save
    if strlength(out_file) == 0
        [folder, name, ~] = fileparts(csv_file);
        out_file = fullfile(folder, char(name + "_qdep_" + suffix + ".png"));
    end
    exportgraphics(gcf, out_file, "Resolution", 220);
    fprintf("Wrote: %s\n", out_file);
end
end

function csv_file = latest_chi_file()
root_dir = fileparts(fileparts(mfilename("fullpath")));
default_data_dir = fullfile(root_dir, "data");
files = dir(fullfile(default_data_dir, "chi_*.csv"));
if isempty(files)
    [file, folder] = uigetfile( ...
        {"*.csv", "Susceptibility CSV files (*.csv)"; "*.*", "All files (*.*)"}, ...
        "Select susceptibility CSV file", ...
        fullfile(default_data_dir, "*.csv"));
    if isequal(file, 0)
        error("No susceptibility CSV file selected.");
    end
    csv_file = fullfile(folder, file);
    return;
end
[~, idx] = max([files.datenum]);
csv_file = fullfile(files(idx).folder, files(idx).name);
end

function [qx_grid, qy_grid, z_grid] = regular_grid(qx, qy, z)
qx_unique = unique(qx(:).');
qy_unique = unique(qy(:));
z_grid = nan(numel(qy_unique), numel(qx_unique));

for i = 1:numel(z)
    ix = find(abs(qx_unique - qx(i)) < 1e-12, 1);
    iy = find(abs(qy_unique - qy(i)) < 1e-12, 1);
    if isempty(ix) || isempty(iy)
        error("Failed to map q point to regular grid.");
    end
    z_grid(iy, ix) = z(i);
end

[qx_grid, qy_grid] = meshgrid(qx_unique, qy_unique);
end

function s = strip_path(path_in)
[~, name, ext] = fileparts(path_in);
s = [name ext];
end
