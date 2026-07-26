%% visualize_orbitals.m
% Loads orbitals.mat (produced by julia/hydrogen_orbitals.jl) and renders
% each hydrogen orbital as a 3D isosurface with lobes colored by the sign
% of the wavefunction (the classic textbook "plus/minus lobe" picture),
% then saves individual screenshots plus one composite gallery figure.
%
% Outputs (into ../outputs/):
%   orbital_<name>.png       -- one high-res screenshot per orbital
%   orbital_gallery.png      -- composite tiled figure, poster-style

clear; close all; clc;

crimson = [0.647 0.110 0.188];   % Harvard Crimson (positive lobe)
gold    = [0.780 0.647 0.016];   % positive/negative accent (negative lobe)
darkbg  = [0.06 0.06 0.08];

dataFile = fullfile('..','julia','orbitals.mat');
outDir   = fullfile('..','outputs');
if ~exist(outDir, 'dir'); mkdir(outDir); end

S = load(dataFile);

% Hardcoded to match julia/hydrogen_orbitals.jl's orbital list exactly --
% avoids relying on how MAT.jl serializes a Julia Vector{String} (it can
% come back as a cell array, a string array, or a padded char matrix
% depending on version, which is what caused the brace-indexing error).
names = {'1s', '2s', '2pz', '3dz2', '4fz3'};

nOrb = numel(names);

%% ---- Individual high-res screenshots ----
for i = 1:nOrb
    name = names{i};
    psi = S.(['psi_' name]);
    coords = S.(['coords_' name]);

    [X, Y, Z] = meshgrid(coords, coords, coords);
    psiPerm = permute(psi, [2 1 3]);

    fig = figure('Color', darkbg, 'Position', [100 100 800 700], 'Visible', 'off');
    ax = axes('Parent', fig, 'Color', darkbg);
    hold(ax, 'on'); axis(ax, 'equal'); axis(ax, 'vis3d'); axis(ax, 'off');
    view(ax, 35, 22);

    absMax = max(abs(psiPerm(:)));
    isoVal = 0.15 * absMax;   % render at 15% of peak amplitude

    fvPos = isosurface(X, Y, Z, psiPerm, isoVal);
    fvNeg = isosurface(X, Y, Z, psiPerm, -isoVal);

    if ~isempty(fvPos.vertices)
        patch(ax, fvPos, 'FaceColor', crimson, 'EdgeColor', 'none', ...
            'FaceAlpha', 0.9, 'FaceLighting', 'gouraud', 'SpecularStrength', 0.5);
    end
    if ~isempty(fvNeg.vertices)
        patch(ax, fvNeg, 'FaceColor', gold, 'EdgeColor', 'none', ...
            'FaceAlpha', 0.9, 'FaceLighting', 'gouraud', 'SpecularStrength', 0.5);
    end

    camlight(ax, 'headlight'); camlight(ax, 'left'); lighting(ax, 'gouraud');
    title(ax, sprintf('Hydrogen Orbital: %s', upper(name)), ...
        'Color', 'w', 'FontSize', 15, 'FontName', 'Georgia');

    fname = fullfile(outDir, sprintf('orbital_%s.png', name));
    exportgraphics(fig, fname, 'Resolution', 300, 'BackgroundColor', darkbg);
    fprintf('Saved screenshot: %s\n', fname);
    close(fig);
end

%% ---- Composite poster-style gallery ----
figGallery = figure('Color', darkbg, 'Position', [50 50 1600 500]);
t = tiledlayout(figGallery, 1, nOrb, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:nOrb
    name = names{i};
    psi = S.(['psi_' name]);
    coords = S.(['coords_' name]);
    [X, Y, Z] = meshgrid(coords, coords, coords);
    psiPerm = permute(psi, [2 1 3]);

    ax = nexttile(t);
    set(ax, 'Color', darkbg); hold(ax, 'on'); axis(ax, 'equal'); axis(ax, 'vis3d'); axis(ax, 'off');
    view(ax, 35, 22);

    absMax = max(abs(psiPerm(:)));
    isoVal = 0.15 * absMax;

    fvPos = isosurface(X, Y, Z, psiPerm, isoVal);
    fvNeg = isosurface(X, Y, Z, psiPerm, -isoVal);
    if ~isempty(fvPos.vertices)
        patch(ax, fvPos, 'FaceColor', crimson, 'EdgeColor', 'none', 'FaceAlpha', 0.9, 'FaceLighting', 'gouraud');
    end
    if ~isempty(fvNeg.vertices)
        patch(ax, fvNeg, 'FaceColor', gold, 'EdgeColor', 'none', 'FaceAlpha', 0.9, 'FaceLighting', 'gouraud');
    end
    camlight(ax, 'headlight'); lighting(ax, 'gouraud');
    title(ax, upper(name), 'Color', 'w', 'FontName', 'Georgia', 'FontSize', 13);
end

title(t, 'Hydrogen Atom Orbital Gallery', 'Color', 'w', 'FontSize', 20, 'FontName', 'Georgia', 'FontWeight', 'bold');

fname = fullfile(outDir, 'orbital_gallery.png');
exportgraphics(figGallery, fname, 'Resolution', 300, 'BackgroundColor', darkbg);
fprintf('Saved gallery: %s\n', fname);