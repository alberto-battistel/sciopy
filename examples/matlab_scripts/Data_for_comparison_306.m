clear; close all; clc;

%% ============================================================
%  1. Initialization
% ============================================================

init_eidors();

%% ============================================================
%  2. User settings
% ============================================================

base_dir = 'C:\Users\Acer\Documents\GitHub\sciopy\examples\Data';

output_dir = fullfile(base_dir, 'reconstruction_results_34');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Static measurements: average all frames before reconstruction
use_average_frames = true;

% Measurement mode
measure_mode = 1;  % 1 = single-ended, 2 = differential

% Physical parameters, only for documentation / later interpretation
tank_radius_m      = 0.100;   % 100 mm
tank_full_height_m = 0.152;   % 152 mm
water_height_m     = 0.120;   % 120 mm

z_lower_ring_m = 0.055;       % 55 mm
z_upper_ring_m = 0.105;       % 105 mm

object_radius_m = 0.035;      % 35 mm
object_height_m = 0.120;      % 120 mm

sigma_water = 1.5;

%% ============================================================
%  3. Define measurement patterns
% ============================================================

patterns = struct([]);

% Pattern 1: triangle
patterns(1).name = 'triangle';
patterns(1).ref_folder  = fullfile(base_dir, 'triangle_ref_measurements');
patterns(1).anom_folder = fullfile(base_dir, 'triangle_anom_measurements');
patterns(1).ref_prefix  = 'triangle_ref_skip';
patterns(1).anom_prefix = 'triangle_anom_skip';

patterns(1).upper_new_numbers = ...
    [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31];

patterns(1).lower_new_numbers = ...
    [32 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30];

% Pattern 2: square
patterns(2).name = 'square';
patterns(2).ref_folder  = fullfile(base_dir, 'square_ref_measurements');
patterns(2).anom_folder = fullfile(base_dir, 'square_anom_measurements');
patterns(2).ref_prefix  = 'square_ref_skip';
patterns(2).anom_prefix = 'square_anom_skip';

patterns(2).upper_new_numbers = ...
    [1 2 5 6 9 10 13 14 17 18 21 22 25 26 29 30];

patterns(2).lower_new_numbers = ...
    [32 3 4 7 8 11 12 15 16 19 20 23 24 27 28 31];

%% ============================================================
%  4. Main loop: 2 patterns x 17 skips = 34 files
% ============================================================

numOfChannels_expected = 32;
skip_values = [0:14 16];

result_counter = 0;

for p = 1:length(patterns)

    pattern_name = patterns(p).name;

    upper_new_numbers = patterns(p).upper_new_numbers;
    lower_new_numbers = patterns(p).lower_new_numbers;

    fprintf('\n============================================================\n');
    fprintf('Processing pattern: %s\n', pattern_name);
    fprintf('============================================================\n');

    for skip_id = skip_values

        fprintf('\n--- Processing %s, skip%02d ---\n', pattern_name, skip_id);

        %% ====================================================
        %  4.1 Build file paths
        % ====================================================

        skip_str = sprintf('%02d', skip_id);

        ref_file = fullfile( ...
            patterns(p).ref_folder, ...
            [patterns(p).ref_prefix, skip_str, '.mat'] ...
        );

        anom_file = fullfile( ...
            patterns(p).anom_folder, ...
            [patterns(p).anom_prefix, skip_str, '.mat'] ...
        );

        if ~exist(ref_file, 'file')
            warning('Reference file not found: %s', ref_file);
            continue;
        end

        if ~exist(anom_file, 'file')
            warning('Anomaly file not found: %s', anom_file);
            continue;
        end

        %% ====================================================
        %  4.2 Load data
        % ====================================================

        imported_data_1 = load(ref_file);
        imported_data_2 = load(anom_file);

        Raw_data_Ref  = permute(imported_data_1.data, [2 3 1]);
        Raw_data_Anom = permute(imported_data_2.data, [2 3 1]);

        amplitude = imported_data_1.amplitude;
        n_frames_ref  = imported_data_1.n_frames;
        n_frames_anom = imported_data_2.n_frames;

        numOfChannels = double(imported_data_1.numOfChannels);
        NSkip = double(imported_data_2.NSkip);

        if numOfChannels ~= numOfChannels_expected
            error('Expected 32 electrodes, but got %d.', numOfChannels);
        end

        if NSkip ~= skip_id
            warning('Filename says skip%02d, but file contains NSkip = %d.', ...
                    skip_id, NSkip);
        end

        VoltageRef = Raw_data_Ref;
        VoltageAnoMoving = Raw_data_Anom;

        %% ====================================================
        %  4.3 Convert Sciospec data to EIT data
        % ====================================================

        rmv_indx = func_rmv_skip(numOfChannels, NSkip);

        VeitRef = nan(numOfChannels*(numOfChannels-3), n_frames_ref);

        for k = 1:n_frames_ref
            if measure_mode == 1
                VeitRef(:,k) = func_ConvertSciospecToEIT( ...
                    VoltageRef(:,:,k).', ...
                    numOfChannels, ...
                    NSkip, ...
                    false ...
                );
            elseif measure_mode == 2
                tmp = VoltageRef(:,:,k);
                tmp(rmv_indx) = [];
                VeitRef(:,k) = tmp(:);
            else
                error('Unknown measure_mode.');
            end
        end

        VeitAnoMoving = nan(numOfChannels*(numOfChannels-3), n_frames_anom);

        for k = 1:n_frames_anom
            if measure_mode == 1
                VeitAnoMoving(:,k) = func_ConvertSciospecToEIT( ...
                    VoltageAnoMoving(:,:,k).', ...
                    numOfChannels, ...
                    NSkip, ...
                    false ...
                );
            elseif measure_mode == 2
                tmp = VoltageAnoMoving(:,:,k);
                tmp(rmv_indx) = [];
                VeitAnoMoving(:,k) = tmp(:);
            else
                error('Unknown measure_mode.');
            end
        end

        %% ====================================================
        %  4.4 Create 3D EIDORS model
        % ====================================================

        imdl3D = mk_common_model('b3cr', [16,2]);
        fmdl = imdl3D.fwd_model;

        nelec = length(fmdl.electrode);

        if nelec ~= 32
            error('Expected 32 electrodes, but model has %d electrodes.', nelec);
        end

        imdl3D.jacobian_bkgnd.value = sigma_water;

        %% ====================================================
        %  4.5 Compute original electrode centers
        % ====================================================

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

        %% ====================================================
        %  4.6 Split electrodes into upper and lower rings
        % ====================================================

        z_median = median(z);

        lower_old_idx = find(z < z_median);
        upper_old_idx = find(z > z_median);

        if length(lower_old_idx) ~= 16 || length(upper_old_idx) ~= 16
            error('Could not split electrodes into upper/lower rings of 16 electrodes.');
        end

        %% ====================================================
        %  4.7 Sort electrodes inside each ring by angle
        % ====================================================

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

        %% ====================================================
        %  4.8 Build remap vector
        % ====================================================

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

        fmdl.electrode = fmdl.electrode(remap);

        %% ====================================================
        %  4.9 Apply stimulation pattern
        % ====================================================

        amplitude = 0.01;

        [fmdl.stimulation, fmdl.meas_select] = ...
            mk_stim_patterns(numOfChannels, 1, ...
                             [1+NSkip, 0], ...
                             [0, 1+NSkip], ...
                             {'no_meas_current'}, ...
                             amplitude*1000);

        imdl3D.fwd_model = fmdl;

        %% ====================================================
        %  4.10 Difference solver
        % ====================================================

        imdl3D.solve = @inv_solve_diff_GN_one_step;
        imdl3D.reconst_type = 'difference';
        imdl3D.RtR_prior = @prior_noser;
        imdl3D.hyperparameter.value = 0.17;

        %% ====================================================
        %  4.11 Static data averaging
        % ====================================================

        if use_average_frames
            vh = mean(real(VeitRef), 2);
            vi = mean(real(VeitAnoMoving), 2);
        else
            vh = real(VeitRef(:,1));
            vi = real(VeitAnoMoving);
        end

        fprintf('Size of vh: %d x %d\n', size(vh,1), size(vh,2));
        fprintf('Size of vi: %d x %d\n', size(vi,1), size(vi,2));

        %% ====================================================
        %  4.12 Reconstruction
        % ====================================================

        img = inv_solve(imdl3D, vh, vi);

        %% ====================================================
        %  4.13 Save only two arrays
        % ====================================================

        rec_elem_data_real = img.elem_data;
        rec_elem_centers = interp_mesh(fmdl);

        output_file = fullfile( ...
            output_dir, ...
            sprintf('%s_reconstruction_skip%02d.mat', pattern_name, skip_id) ...
        );

        save(output_file, ...
             'rec_elem_centers', ...
             'rec_elem_data_real');

        result_counter = result_counter + 1;

        fprintf('Saved: %s\n', output_file);

    end
end

fprintf('\n============================================================\n');
fprintf('Finished. Saved %d reconstruction files.\n', result_counter);
fprintf('Output folder:\n%s\n', output_dir);
fprintf('============================================================\n');