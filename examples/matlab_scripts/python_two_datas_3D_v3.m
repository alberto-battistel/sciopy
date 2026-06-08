clear; close all; clc;

%% 1.Initialization
% Start-up EIDORS
init_eidors()

%Upload files from the ScioPy
%%Import refernce.mat
imported_data_1 = load('C:\Users\Acer\Documents\GitHub\sciopy\examples\Data\square_ref_measurements\square_ref_skip03.mat'); 
%%Import data with anomaly.mat
imported_data_2 = load('C:\Users\Acer\Documents\GitHub\sciopy\examples\Data\square_anom_measurements\square_anom_skip03.mat');

%% 2.Data processing
Raw_data_Ref = permute(imported_data_1.data, [2 3 1]);
amplitude = imported_data_1.amplitude;
n_frames = imported_data_1.n_frames;
numOfChannels = double(imported_data_1.numOfChannels); %electrodes
NSkip = double(imported_data_2.NSkip); %inj_skip
measure_mode = 1;  % single-ended mode
Raw_data_Anom = permute(imported_data_2.data, [2 3 1]);
numOfRefFrames = n_frames;
numOfImagingFrames = n_frames;

%Devide data into  and data frames
VoltageRef  = Raw_data_Ref;
VoltageAnoMoving = Raw_data_Anom;

%% 3.Convert Sciospec data to EIT data
rmv_indx = func_rmv_skip(numOfChannels, NSkip);
VeitRef = nan(numOfChannels*(numOfChannels-3), numOfRefFrames);
for k=1:numOfRefFrames
    if measure_mode == 1 % single-ended mode
        VeitRef(:,k) = func_ConvertSciospecToEIT(VoltageRef(:,:,k).',numOfChannels,NSkip,false);
    elseif measure_mode == 2 % differntial mode
        tmp = VoltageRef(:,:,k);
        tmp(rmv_indx) = [];
        VeitRef(:,k) = tmp(:);
    end
    
end

VeitAnoMoving = nan(numOfChannels*(numOfChannels-3), numOfImagingFrames);
for k=1:numOfImagingFrames
    if measure_mode == 1 % single-ended mode
        VeitAnoMoving(:,k) = func_ConvertSciospecToEIT(VoltageAnoMoving(:,:,k).',numOfChannels,NSkip,false);
    elseif measure_mode == 2 % differntial mode
        tmp = VoltageAnoMoving(:,:,k);
        tmp(rmv_indx) = [];
        VeitAnoMoving(:,k) = tmp(:);
    end
end

v_all = real(VeitAnoMoving);
v_ref = real(VeitRef(:,1)); % Set first frame as the reference frame for time-difference imaging
figure(1);clf;plot(abs(v_all),'.-'); title('U-shape Plot'); grid on; ylabel('abs(voltage) [V]'); set(gca,'yscale','log') % display U-shape plot for quality checking

%% 4. 3D model setup

%% 4.1 Create 3D model using previous working EIDORS model

% ============================================================
% Real physical tank parameters
% These parameters are used for documentation, scaling and saving.
% The actual FEM model remains the previous normalized b3cr model.
% ============================================================

tank_radius_m      = 0.100;   % 100 mm
tank_full_height_m = 0.152;   % 152 mm
water_height_m     = 0.120;   % 120 mm conducting saline height

z_lower_ring_m = 0.055;       % lower electrode row, 55 mm
z_upper_ring_m = 0.105;       % upper electrode row, 105 mm

object_radius_m = 0.035;      % 35 mm
object_height_m = 0.120;      % 120 mm

% Approximate saline conductivity
sigma_water = 1.5;

% ============================================================
% Previous working EIDORS model
% ============================================================

imdl3D = mk_common_model('b3cr', [16,2]);
fmdl = imdl3D.fwd_model;

nelec = length(fmdl.electrode);

if nelec ~= 32
    error('Expected 32 electrodes, but model has %d electrodes.', nelec);
end

% Reconstruction background conductivity
imdl3D.jacobian_bkgnd.value = sigma_water;

%% 4.2 Compute original electrode centers

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

%% 4.3 Split electrodes into upper and lower rings

z_median = median(z);
lower_old_idx = find(z < z_median);
upper_old_idx = find(z > z_median);
if length(lower_old_idx) ~= 16 || length(upper_old_idx) ~= 16
    error('Could not split electrodes into upper/lower rings of 16 electrodes.');
end


%% 4.4 Sort electrodes inside each ring by angle

% direction = +1 means increasing angle counter-clockwise
% direction = -1 means clockwise
%
% start_angle_deg controls where the first electrode in each ring starts.
% 0 deg   = +x direction
% 90 deg  = +y direction
% 180 deg = -x direction
% 270 deg = -y direction

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


%% 4.5 Define your required new electrode numbers
%% triangle
upper_new_numbers = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31];

lower_new_numbers = [32 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30];
% %% Square
% upper_new_numbers = [1 2 5 6 9 10 13 14 17 18 21 22 25 26 29 30];
% 
% lower_new_numbers = [32 3 4 7 8 11 12 15 16 19 20 23 24 27 28 31];

if length(upper_new_numbers) ~= 16 || length(lower_new_numbers) ~= 16
    error('Upper and lower numbering vectors must both have length 16.');
end

if ~isequal(sort([upper_new_numbers lower_new_numbers]), 1:32)
    error('The new numbering must contain every number from 1 to 32 exactly once.');
end


%% 4.6 Build remap vector
% Meaning:
% remap(new_electrode_number) = old_electrode_number
%
% Example:
% remap(1) = old electrode physically located at the first upper position
% remap(3) = old electrode physically located at the first lower position

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


%% 4.7 Apply remapping to the model

fmdl.electrode = fmdl.electrode(remap);


%% 4.8 Apply continuous current injection pattern after remapping


numOfChannels = 32;
amplitude = 0.01;   % amplitude*1000 = 10

[fmdl.stimulation, fmdl.meas_select] = ...
    mk_stim_patterns(numOfChannels, 1, ...
                     [1+NSkip, 0], ...
                     [0, 1+NSkip], ...
                     {'no_meas_current'}, ...
                     amplitude*1000);

imdl3D.fwd_model = fmdl;


%% 4.9 Compute new electrode centers after remapping

elec_centers_new = zeros(nelec, 3);

for i = 1:nelec
    nodes_i = fmdl.electrode(i).nodes;
    coords_i = fmdl.nodes(nodes_i, :);
    elec_centers_new(i,:) = mean(coords_i, 1);
end

angles_new_deg = mod(rad2deg(atan2(elec_centers_new(:,2), ...
                                    elec_centers_new(:,1))), 360);


%% 5. Difference solver
imdl3D.solve = @inv_solve_diff_GN_one_step;
imdl3D.reconst_type = 'difference';
imdl3D.RtR_prior = @prior_noser;      % как в tutorial
imdl3D.hyperparameter.value = 0.17;   % стартовое значение из tutorial

%% 6. Data
% Static measurement: average all frames before reconstruction

vh = mean(real(VeitRef), 2);          % averaged reference, 928x1
vi = mean(real(VeitAnoMoving), 2);    % averaged anomaly, 928x1

fprintf('Size of averaged vh: %d x %d\n', size(vh,1), size(vh,2));
fprintf('Size of averaged vi: %d x %d\n', size(vi,1), size(vi,2));

%% 7. Reconstruction
img = inv_solve(imdl3D, vh, vi);

%% ============================================================
%  8. Visualization: 3D FEM reconstruction
% ============================================================

img.calc_colours.ref_level = 0;
img.calc_colours.greylev = -0.05;

figure(2); clf;
show_fem(img);
title('3D EIDORS reconstruction with custom electrode renumbering');
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

hold off;
axis equal;
view(3);
eidors_colourbar(img);

%% ============================================================
%  9. show_slices visualization
% ============================================================

% Model limits in normalized EIDORS coordinates
x_min = min(fmdl.nodes(:,1));
x_max = max(fmdl.nodes(:,1));

y_min = min(fmdl.nodes(:,2));
y_max = max(fmdl.nodes(:,2));

z_min = min(fmdl.nodes(:,3));
z_max = max(fmdl.nodes(:,3));

x0 = 0.5 * (x_min + x_max);
y0 = 0.5 * (y_min + y_max);
z0 = 0.5 * (z_min + z_max);

fprintf('\nNormalized model limits:\n');
fprintf('x: %.4f to %.4f\n', x_min, x_max);
fprintf('y: %.4f to %.4f\n', y_min, y_max);
fprintf('z: %.4f to %.4f\n', z_min, z_max);

% Horizontal slices at 25%, 50%, 75% of model height
z1 = z_min + 0.25*(z_max - z_min);
z2 = z_min + 0.50*(z_max - z_min);
z3 = z_min + 0.75*(z_max - z_min);

horizontal_slices = [
    inf, inf, z1;
    inf, inf, z2;
    inf, inf, z3
];

figure(3); clf;
d123 = show_slices(img, horizontal_slices);
eidors_colourbar(img);
title('Real reconstruction: horizontal show\_slices');

% Central vertical slices
vertical_slices_center = [
    x0,  inf, inf;
    inf, y0,  inf
];

figure(4); clf;
show_slices(img, vertical_slices_center);
eidors_colourbar(img);
title('Real reconstruction: central vertical show\_slices');

%% ============================================================
%  10. show_slices through expected object position
% ============================================================

% Original manually estimated normalized object center
obj_cx_norm_original = -0.3;
obj_cy_norm_original = -0.3625;

% Move object estimate inward if needed
% 1.00 = original manually estimated position
% 0.75 = 25% closer to center
% 0.65 = 35% closer to center
move_to_center_factor = 0.75;

obj_cx_norm = obj_cx_norm_original * move_to_center_factor;
obj_cy_norm = obj_cy_norm_original * move_to_center_factor;

obj_radius_norm = object_radius_m / tank_radius_m;

fprintf('\nExpected object position in normalized coordinates:\n');
fprintf('obj_cx_norm = %.4f\n', obj_cx_norm);
fprintf('obj_cy_norm = %.4f\n', obj_cy_norm);
fprintf('obj_radius_norm = %.4f\n', obj_radius_norm);

vertical_slices_object = [
    obj_cx_norm, inf, inf;
    inf, obj_cy_norm, inf
];

figure(5); clf;
show_slices(img, vertical_slices_object);
eidors_colourbar(img);
title('Real reconstruction: vertical show\_slices through expected object center');

%% ============================================================
%  11. Optional old-style 3D slices
% ============================================================

figure(6); clf;
show_3d_slices(img, obj_cx_norm, obj_cy_norm, z0);
eidors_colourbar(img);
title('Real reconstruction: show\_3d\_slices through expected object center');
axis equal;
view(3);

%% ============================================================
%  12. Save reconstruction data for comparison
% ============================================================

% Main reconstruction values on FEM elements
rec_elem_data_real = img.elem_data;

% FEM element centers in normalized coordinates
rec_elem_centers_norm = interp_mesh(fmdl);

% Convert normalized coordinates to approximate physical coordinates
% x,y: normalized radius 1 corresponds to tank_radius_m
% z: map model z-range to water height 0...water_height_m

rec_elem_centers_m = rec_elem_centers_norm;

rec_elem_centers_m(:,1) = rec_elem_centers_norm(:,1) * tank_radius_m;
rec_elem_centers_m(:,2) = rec_elem_centers_norm(:,2) * tank_radius_m;

rec_elem_centers_m(:,3) = ...
    (rec_elem_centers_norm(:,3) - z_min) / (z_max - z_min) * water_height_m;

% FEM mesh data
rec_nodes_norm = fmdl.nodes;
rec_elems = fmdl.elems;

rec_nodes_m = rec_nodes_norm;
rec_nodes_m(:,1) = rec_nodes_norm(:,1) * tank_radius_m;
rec_nodes_m(:,2) = rec_nodes_norm(:,2) * tank_radius_m;
rec_nodes_m(:,3) = ...
    (rec_nodes_norm(:,3) - z_min) / (z_max - z_min) * water_height_m;

% Voltage-level data
vh_real = vh;
vi_real = vi;
dV_real = vi - vh;

%% ============================================================
%  Save only reconstruction data for comparison with ideal case
% ============================================================

rec_elem_data_real = img.elem_data;
rec_elem_centers = interp_mesh(fmdl);

save('real_reconstruction_arrays.mat', ...
     'rec_elem_centers', ...
     'rec_elem_data_real');