clear; close all; clc;

%% ============================================================
%  1. Initialization
% ============================================================

init_eidors();

%% ============================================================
%  2. Load real measurement files ONLY to get NSkip and settings
% ============================================================
% Здесь мы используем реальные файлы не для реконструкции,
% а чтобы взять тот же NSkip, numOfChannels и amplitude.

imported_data_1 = load('C:\Users\Acer\Documents\GitHub\sciopy\examples\Data\square_ref_measurements\square_ref_skip07.mat'); 
imported_data_2 = load('C:\Users\Acer\Documents\GitHub\sciopy\examples\Data\square_anom_measurements\square_anom_skip07.mat');

amplitude     = imported_data_1.amplitude;
n_frames      = imported_data_1.n_frames;
numOfChannels = double(imported_data_1.numOfChannels);
NSkip         = double(imported_data_2.NSkip);

fprintf('numOfChannels = %d\n', numOfChannels);
fprintf('NSkip = %d\n', NSkip);
fprintf('amplitude from file = %g\n', amplitude);

if numOfChannels ~= 32
    error('Expected 32 electrodes, but got %d.', numOfChannels);
end

%% ============================================================
%  3. Create same 3D EIDORS model as in real reconstruction
% ============================================================

imdl3D = mk_common_model('b3cr', [16,2]);
fmdl = imdl3D.fwd_model;

nelec = length(fmdl.electrode);

if nelec ~= 32
    error('Expected 32 electrodes, but model has %d electrodes.', nelec);
end

%% ============================================================
%  4. Compute original electrode centers
% ============================================================

elec_centers_old = zeros(nelec, 3);

for i = 1:nelec
    nodes_i = fmdl.electrode(i).nodes;
    coords_i = fmdl.nodes(nodes_i, :);
    elec_centers_old(i,:) = mean(coords_i, 1);
end

x = elec_centers_old(:,1);
y = elec_centers_old(:,2);
z = elec_centers_old(:,3);

angles_deg = mod(rad2deg(atan2(y, x)), 360);

%% ============================================================
%  5. Split electrodes into upper and lower rings
% ============================================================

z_median = median(z);

lower_old_idx = find(z < z_median);
upper_old_idx = find(z > z_median);

if length(lower_old_idx) ~= 16 || length(upper_old_idx) ~= 16
    error('Could not split electrodes into upper/lower rings of 16 electrodes.');
end

%% ============================================================
%  6. Sort electrodes inside each ring by angle
% ============================================================

direction = -1;
start_angle_deg = 0;

upper_angles = angles_deg(upper_old_idx);
lower_angles = angles_deg(lower_old_idx);

if direction == +1
    [~, upper_order] = sort(mod(upper_angles - start_angle_deg, 360), 'ascend');
    [~, lower_order] = sort(mod(lower_angles - start_angle_deg, 360), 'ascend');
else
    [~, upper_order] = sort(mod(start_angle_deg - upper_angles, 360), 'ascend');
    [~, lower_order] = sort(mod(start_angle_deg - lower_angles, 360), 'ascend');
end

upper_sorted_old_idx = upper_old_idx(upper_order);
lower_sorted_old_idx = lower_old_idx(lower_order);

%% ============================================================
%  7. Apply same custom electrode numbering as in real code
% ============================================================

upper_new_numbers = [1 2 5 6 9 10 13 14 17 18 21 22 25 26 29 30];

% ВАЖНО:
% Я оставляю именно твою версию из оригинального кода.
% Там lower ring начинается с 32.
lower_new_numbers = [32 3 4 7 8 11 12 15 16 19 20 23 24 27 28 31];

if length(upper_new_numbers) ~= 16 || length(lower_new_numbers) ~= 16
    error('Upper and lower numbering vectors must both have length 16.');
end

if ~isequal(sort([upper_new_numbers lower_new_numbers]), 1:32)
    error('The new numbering must contain every number from 1 to 32 exactly once.');
end

remap = zeros(1, nelec);

for k = 1:16
    remap(upper_new_numbers(k)) = upper_sorted_old_idx(k);
    remap(lower_new_numbers(k)) = lower_sorted_old_idx(k);
end

if any(remap == 0)
    error('Remap vector contains zeros. Something went wrong.');
end

if ~isequal(sort(remap), 1:nelec)
    error('Invalid remap vector. It must contain all old electrodes exactly once.');
end

fprintf('\nRemap vector:\n');
disp(remap);

fmdl.electrode = fmdl.electrode(remap);

%% ============================================================
%  8. Apply same stimulation pattern as real measurement
% ============================================================

% В оригинальном коде:
% amplitude = 0.01;
% amplitude*1000 = 10
%
% Я оставляю это поведение, чтобы методика совпадала.

amplitude = 0.01;

[fmdl.stimulation, fmdl.meas_select] = ...
    mk_stim_patterns(numOfChannels, 1, ...
                     [1+NSkip, 0], ...
                     [0, 1+NSkip], ...
                     {'no_meas_current'}, ...
                     amplitude*1000);

imdl3D.fwd_model = fmdl;

%% ============================================================
%  9. Check simulated data size
% ============================================================

n_meas_expected = numOfChannels * (numOfChannels - 3);

fprintf('Expected number of measurements = %d\n', n_meas_expected);

%% ============================================================
%  10. Create ideal conductivity model
% ============================================================

% Conductivity values
% Salt water: approximate value.
% For maximum accuracy, replace this with measured conductivity.
sigma_water = 1.5;    

% Object conductivity.
% If the object is plastic/insulating, keep very low conductivity.
sigma_object = 1e-6;

img_hom = mk_image(fmdl, sigma_water);
img_hom.name = 'Ideal homogeneous salt water';

img_obj = mk_image(fmdl, sigma_water);
img_obj.name = 'Ideal salt water with cylindrical object';

%% ============================================================
%  11. Add cylindrical object
% ============================================================

% Initial manually estimated object center
obj_cx_original = -0.3;
obj_cy_original = -0.3625;

% Move object inward toward model center.
% 1.00 = original position
% 0.75 = 25% closer to center
% 0.65 = even closer to center
move_to_center_factor = 0.75;

obj_cx = obj_cx_original * move_to_center_factor;
obj_cy = obj_cy_original * move_to_center_factor;

% Object radius:
% real radius = 35 mm, tank radius = 100 mm
obj_radius = 0.35;

elem_centers = interp_mesh(fmdl);

ex = elem_centers(:,1);
ey = elem_centers(:,2);
ez = elem_centers(:,3);

elems = fmdl.elems;
nodes = fmdl.nodes;

% Method 1: element center inside object
inside_by_center = (ex - obj_cx).^2 + ...
                   (ey - obj_cy).^2 <= obj_radius^2;

% Method 2: at least one element node inside object
inside_by_nodes = false(size(elems,1),1);

for e = 1:size(elems,1)
    node_ids = elems(e,:);
    node_xy = nodes(node_ids,1:2);

    dist2_nodes = (node_xy(:,1) - obj_cx).^2 + ...
                  (node_xy(:,2) - obj_cy).^2;

    if any(dist2_nodes <= obj_radius^2)
        inside_by_nodes(e) = true;
    end
end

% Combined object mask
inside_object = inside_by_center | inside_by_nodes;

fprintf('Object center: x = %.4f, y = %.4f\n', obj_cx, obj_cy);
fprintf('Object radius: %.4f\n', obj_radius);
fprintf('Distance from center: %.4f\n', sqrt(obj_cx^2 + obj_cy^2));
fprintf('Approx. distance from object edge to tank wall: %.4f\n', ...
        1 - sqrt(obj_cx^2 + obj_cy^2) - obj_radius);
fprintf('Number of elements inside object = %d\n', sum(inside_object));

img_obj = mk_image(fmdl, sigma_water);
img_obj.name = 'Ideal salt water with cylindrical object';
img_obj.elem_data(inside_object) = sigma_object;

%% ============================================================
%  12. Forward solve: ideal measurement
% ============================================================

vh_sim = fwd_solve(img_hom);
vi_sim = fwd_solve(img_obj);

vh = vh_sim.meas;
vi = vi_sim.meas;

fprintf('Size of vh = %d x %d\n', size(vh,1), size(vh,2));
fprintf('Size of vi = %d x %d\n', size(vi,1), size(vi,2));

if size(vh,1) ~= n_meas_expected
    error('vh has wrong size. Expected %d measurements.', n_meas_expected);
end

if size(vi,1) ~= n_meas_expected
    error('vi has wrong size. Expected %d measurements.', n_meas_expected);
end

figure(1);
plot(abs(vi - vh), '.-');
grid on;
set(gca, 'yscale', 'log');
title('Ideal simulated voltage difference');
xlabel('Measurement index');
ylabel('|vi - vh| [V]');

%% ============================================================
%  13. Difference solver: same as real reconstruction
% ============================================================

imdl3D.solve = @inv_solve_diff_GN_one_step;
imdl3D.reconst_type = 'difference';
imdl3D.RtR_prior = @prior_noser;
imdl3D.hyperparameter.value = 0.17;

%% ============================================================
%  14. Reconstruction of ideal simulated data
% ============================================================

img_ideal_rec = inv_solve(imdl3D, vh, vi);

img_ideal_rec.calc_colours.ref_level = 0;
img_ideal_rec.calc_colours.greylev = -0.05;

%% ============================================================
%  15. Compute electrode centers after remapping for visualization
% ============================================================

elec_centers_new = zeros(nelec, 3);

for i = 1:nelec
    nodes_i = fmdl.electrode(i).nodes;
    coords_i = fmdl.nodes(nodes_i, :);
    elec_centers_new(i,:) = mean(coords_i, 1);
end

%% ============================================================
%  16. Visualization: true ideal object
% ============================================================

figure(2); clf;

show_fem(img_obj);
title('True ideal model: cylindrical object in salt water');
axis equal;
view(3);
hold on;

% Model z-limits
z_min = min(fmdl.nodes(:,3));
z_max = max(fmdl.nodes(:,3));

% Smooth cylinder surface
n_cyl = 150;
theta = linspace(0, 2*pi, n_cyl);
z_cyl = linspace(z_min, z_max, 50);

[Theta, Zc] = meshgrid(theta, z_cyl);

Xc = obj_cx + obj_radius*cos(Theta);
Yc = obj_cy + obj_radius*sin(Theta);

hs = surf(Xc, Yc, Zc);
set(hs, ...
    'FaceColor', [0 0.45 0.9], ...
    'FaceAlpha', 0.30, ...
    'EdgeColor', 'none');

% Top and bottom circles
plot3(obj_cx + obj_radius*cos(theta), ...
      obj_cy + obj_radius*sin(theta), ...
      z_min*ones(size(theta)), ...
      'k-', 'LineWidth', 1.5);

plot3(obj_cx + obj_radius*cos(theta), ...
      obj_cy + obj_radius*sin(theta), ...
      z_max*ones(size(theta)), ...
      'k-', 'LineWidth', 1.5);

% Object center line
plot3([obj_cx obj_cx], [obj_cy obj_cy], [z_min z_max], ...
      'k--', 'LineWidth', 1.5);

camlight headlight;
lighting gouraud;
axis vis3d;
hold off;

%% ============================================================
%  17. Visualization: reconstructed ideal object
% ============================================================

figure(3);
show_fem(img_ideal_rec);
title('Ideal reconstruction using same EIDORS method');
hold on;

scale_xy = 1.15;

for i = 1:nelec
    c = elec_centers_new(i,:);
    label_pos = [scale_xy*c(1), scale_xy*c(2), c(3)];

    text(label_pos(1), label_pos(2), label_pos(3), sprintf('%d', i), ...
        'Color', 'b', ...
        'FontSize', 12, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center');
end

% Show true object center line on reconstructed image
z_min = min(fmdl.nodes(:,3));
z_max = max(fmdl.nodes(:,3));

plot3([obj_cx obj_cx], [obj_cy obj_cy], [z_min z_max], ...
      'r--', 'LineWidth', 2);

hold off;
axis equal;
view(3);

%% ============================================================
%  18. Slice visualization
% ============================================================

%% ============================================================
%  Correct slice visualization using actual model coordinates
% ============================================================

x_min = min(fmdl.nodes(:,1));
x_max = max(fmdl.nodes(:,1));

y_min = min(fmdl.nodes(:,2));
y_max = max(fmdl.nodes(:,2));

z_min = min(fmdl.nodes(:,3));
z_max = max(fmdl.nodes(:,3));

x0 = mean([x_min, x_max]);
y0 = mean([y_min, y_max]);
z0 = mean([z_min, z_max]);

fprintf('Model limits:\n');
fprintf('x: %.3f to %.3f, center %.3f\n', x_min, x_max, x0);
fprintf('y: %.3f to %.3f, center %.3f\n', y_min, y_max, y0);
fprintf('z: %.3f to %.3f, center %.3f\n', z_min, z_max, z0);

figure;
show_3d_slices(img_ideal_rec, obj_cx, obj_cy, z0);
eidors_colourbar(img_ideal_rec);
title('Ideal reconstruction: slices through object center');
axis equal;
view(3);

%% ============================================================
%  18. Slice visualization with show_slices
% ============================================================

z_min = min(fmdl.nodes(:,3));
z_max = max(fmdl.nodes(:,3));

z1 = z_min + 0.25*(z_max - z_min);
z2 = z_min + 0.50*(z_max - z_min);
z3 = z_min + 0.75*(z_max - z_min);

horizontal_slices = [
    inf, inf, z1;
    inf, inf, z2;
    inf, inf, z3
];

vertical_slices = [
    obj_cx, inf, inf;
    inf, obj_cy, inf
];

figure;
show_slices(img_ideal_rec, horizontal_slices);
eidors_colourbar(img_ideal_rec);
title('Ideal reconstruction: horizontal slices');

figure;
show_slices(img_ideal_rec, vertical_slices);
eidors_colourbar(img_ideal_rec);
title('Ideal reconstruction: vertical slices through object center');


%% ============================================================
%  19. Save ideal voltages if needed
% ============================================================
rec_elem_data_ideal = img_ideal_rec.elem_data;
rec_elem_centers = interp_mesh(fmdl);

save('ideal_reconstruction_arrays.mat', ...
     'rec_elem_data_ideal', ...
     'rec_elem_centers');

