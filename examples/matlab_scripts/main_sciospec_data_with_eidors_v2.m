%% Start-up EIDORS
run C:\EIDORS\eidors-v3.12\eidors\startup.m
%% Select example dataset
numExp = 1;
switch numExp
    case 1 % Adjacent for current injection and voltage measurement pattern
        fpath ='NewScript\Data\20260504 11.34.20'; setupName = 'setup'; numOfRefFrames = 20; numOfImagingFrames = 80;
        NSkip = 0; % skip pattern for the current injection, 0 is the adjacent(or neighboring) injection pattern  
        channels = 1:32; % use the channels 1,2,.., 16 for the image reconstruction
    % case 2 % 2-skip for current injection and voltage measurement pattern
    %     fpath ='WaterTankData/20250212 13.21.48'; setupName = 'setup'; numOfRefFrames = 20; numOfImagingFrames = 270;      
    %     NSkip = 2;
    %     channels = 1:16;
end

%% Load Sciospec EIT data
numOfChannels = length(channels);

clear VoltageAnoMany_temp VoltageRef_temp

% Reference, the water tank without object
VoltageRef = nan(numOfChannels, numOfChannels, numOfRefFrames);
for i=1:numOfRefFrames
    fname = [setupName '_'  sprintf('%.5d', i) '.eit'];
    sciospecData = fnc_read_SciospecData(fullfile(fpath,setupName,fname));
    voltages = sciospecData.Voltages.voltage(:,channels);
    VoltageRef(:,:,i) = reshape(voltages(:),numOfChannels,numOfChannels);
end

% Anomaly, the water tank with an insulating object (a small cup of class).
VoltageAnoMoving = nan(numOfChannels, numOfChannels, numOfImagingFrames);
for i=1:numOfImagingFrames
    fname = [setupName '_'  sprintf('%.5d', i+numOfRefFrames) '.eit'];
    sciospecData = fnc_read_SciospecData(fullfile(fpath,setupName,fname));
    voltages = sciospecData.Voltages.voltage(:,channels);
    VoltageAnoMoving(:,:,i) = reshape(voltages(:),numOfChannels,numOfChannels);
end

amplitudeStr = sciospecData.Amplitude;
disp(['Injected Current amplitude : ' amplitudeStr])
if is_octave()
    tmp = strsplit(amplitudeStr);
else
    tmp = split(amplitudeStr);
end
amplitude = str2double(tmp{1});
%% Check Mesaure Mode
measure_mode = fnc_get_measure_mode(fullfile(fpath,setupName,[setupName '.setUp']));
if measure_mode ~= 1 &&  measure_mode ~= 2
    error(['Not supported measure mode, measure_mode: ' num2str(measure_mode)])
end
%% Convert Sciospec data to EIT data
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


%% Make foward and invese models using eidors
imdl2D = mk_common_model('f2c',32); %  2D circular model
figure(2);clf;show_fem(imdl2D.fwd_model);
title('FEM model')
[imdl2D.fwd_model.stimulation,imdl2D.fwd_model.meas_select] = mk_stim_patterns(numOfChannels,1,[1+NSkip,0],[0, 1+NSkip],{'no_meas_current'},amplitude*1000);
%% Image reconstruction
img2D = inv_solve_diff_GN_one_step(imdl2D, v_ref, v_all);

img2D.calc_colours.ref_level=0;
img2D.type='image';
img2D.show_slices.img_cols=0;
img2D.calc_colours.ref_level=0;

figure(3);clf;show_slices(img2D);
title('Difference Imaging')
%% Image reconstruction (Movie)
img2D = inv_solve_diff_GN_one_step(imdl2D, v_ref, v_all);
elem_data_all = img2D.elem_data;

figure(4);clf;
for i = 1:size(v_all,2)
    img2D.elem_data = elem_data_all(:,i);
    figure(4);
    show_slices(img2D);
    title({'Difference Imaging' ; [num2str(i) 'th frame']})
    pause(0.01)
end