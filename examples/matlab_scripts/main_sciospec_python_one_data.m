%% Start-up EIDORS
run C:\EIDORS\eidors-v3.12\eidors\startup.m

%Upload file from the ScioPy

imported_data = load('C:\Users\Acer\Documents\GitHub\sciopy\examples\experiment_target_2.mat');
Raw_data_Sciopy = permute(imported_data.data, [2 3 1]);
amplitude = imported_data.amplitude;
n_frames = imported_data.n_frames;
numOfChannels = double(imported_data.numOfChannels); %electrodes
NSkip = double(imported_data.NSkip); %inj_skip


measure_mode = 1;  % single-ended mode


numOfRefFrames = 20;
numOfImagingFrames = n_frames - numOfRefFrames;

%Devide data into  and data frames


VoltageRef  = Raw_data_Sciopy(:,:,1:numOfRefFrames);
VoltageAnoMoving = Raw_data_Sciopy(:,:,numOfRefFrames+1:n_frames);


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