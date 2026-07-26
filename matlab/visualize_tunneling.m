%% visualize_tunneling.m
% Loads tunneling_frames.mat (produced by julia/tunneling_solver.jl) and
% renders the 3D probability density |psi|^2 as it hits a potential
% barrier and partially tunnels through.
%
% Outputs (into ../outputs/):
%   tunneling_frame_XX.png   -- individual screenshots (key moments)
%   tunneling_animation.mp4  -- full animation
%
% Run this from the matlab/ folder (or adjust the paths below).

clear; close all; clc;

%% ---- Config ----
crimson = [0.647 0.110 0.188];   % Harvard Crimson
gold    = [0.780 0.647 0.016];
darkbg  = [0.06 0.06 0.08];

dataFile   = fullfile('..','julia','tunneling_frames.mat');
outDir     = fullfile('..','outputs');
if ~exist(outDir, 'dir'); mkdir(outDir); end

%% ---- Load data ----
S = load(dataFile);
psi2 = S.psi2;              % Nx x Ny x Nz x Nframes
x = S.x; y = S.y; z = S.z;
t = S.t;
barrier_x = S.barrier_x;
nFrames = size(psi2, 4);

[X, Y, Z] = meshgrid(y, x, z);   % note: meshgrid swaps first two dims vs ndgrid

%% ---- Figure setup ----
fig = figure('Color', darkbg, 'Position', [100 100 1000 800]);
ax = axes('Parent', fig, 'Color', darkbg);
hold(ax, 'on'); axis(ax, 'equal'); axis(ax, 'vis3d');
xlim(ax, [min(x) max(x)]); ylim(ax, [min(y) max(y)]); zlim(ax, [min(z) max(z)]);
xlabel(ax, 'x (a_0)', 'Color', 'w'); ylabel(ax, 'y (a_0)', 'Color', 'w'); zlabel(ax, 'z (a_0)', 'Color', 'w');
set(ax, 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w', 'GridColor', [0.4 0.4 0.4]);
view(ax, 35, 22);
camlight(ax, 'headlight'); lighting(ax, 'gouraud');
title(ax, 'Quantum Tunneling Through a Potential Barrier', ...
    'Color', 'w', 'FontSize', 16, 'FontName', 'Georgia');

% Draw the potential barrier as a translucent slab
x1 = barrier_x(1); x2 = barrier_x(2);
yb = [min(y) max(y)]; zb = [min(z) max(z)];
barrierVerts = [x1 yb(1) zb(1); x2 yb(1) zb(1); x2 yb(2) zb(1); x1 yb(2) zb(1); ...
                x1 yb(1) zb(2); x2 yb(1) zb(2); x2 yb(2) zb(2); x1 yb(2) zb(2)];
barrierFaces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
patch(ax, 'Vertices', barrierVerts, 'Faces', barrierFaces, ...
    'FaceColor', gold, 'FaceAlpha', 0.12, 'EdgeColor', 'none');

%% ---- Video writer ----
v = VideoWriter(fullfile(outDir, 'tunneling_animation.mp4'), 'MPEG-4');
v.FrameRate = 12;
v.Quality = 95;
open(v);

isoSurfHandle = [];

% Key frames to save as standalone screenshots: before, at, after barrier
screenshotFrames = unique(round([1, nFrames*0.35, nFrames*0.55, nFrames]));
shotLabels = {'before_barrier', 'at_barrier', 'tunneling', 'after_barrier'};

for f = 1:nFrames
    if ~isempty(isoSurfHandle) && isvalid(isoSurfHandle)
        delete(isoSurfHandle);
    end

    P = squeeze(psi2(:,:,:,f));
    Pmax = max(P(:));
    isoVal = 0.10 * Pmax;   % isosurface at 10% of peak density

    Pperm = permute(P, [2 1 3]);  % match meshgrid's (y,x,z) ordering
    fv = isosurface(X, Y, Z, Pperm, isoVal);

    if ~isempty(fv.vertices)
        isoSurfHandle = patch(ax, fv, 'FaceColor', crimson, 'EdgeColor', 'none', ...
            'FaceAlpha', 0.85, 'FaceLighting', 'gouraud', 'AmbientStrength', 0.3, ...
            'SpecularStrength', 0.4);
    else
        isoSurfHandle = patch(ax, 'Vertices', [], 'Faces', []);
    end

    subtitle(ax, sprintf('t = %.2f  (transmission so far grows as packet crosses barrier)', t(f)), ...
        'Color', [0.85 0.85 0.85], 'FontName', 'Georgia');

    drawnow;
    frame = getframe(fig);
    writeVideo(v, frame);

    idx = find(screenshotFrames == f, 1);
    if ~isempty(idx)
        fname = fullfile(outDir, sprintf('tunneling_%s.png', shotLabels{idx}));
        exportgraphics(fig, fname, 'Resolution', 300);
        fprintf('Saved screenshot: %s\n', fname);
    end
end

close(v);
fprintf('Saved animation: %s\n', fullfile(outDir, 'tunneling_animation.mp4'));
fprintf('Transmission probability (from Julia sim): %.4f\n', S.transmission_prob);
