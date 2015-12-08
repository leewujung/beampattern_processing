% 2015 11 24  Use real measurement az/el angles for model beampattern

clear
usrn = getenv('username');
addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\EllipseDirectFit']);
addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\beampattern_processing']);

% Bat data path
save_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing\20151208_model_bp_fit_ellipse'];
if ~exist(save_path,'dir')
    mkdir(save_path);
end
base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing'];
bat_proc_path = './proc_output_checked';
bat_proc_file_all = dir(fullfile(base_path,bat_proc_path,'rousettus_20150825_*'));
freq_wanted = 35e3;

% Set bp param
bp_info.c = 344;  % sound speed [m/s]
bp_info.freq = freq_wanted;  % [Hz]
bp_type = 'piston';
if strcmp(bp_type,'piston');
    bp_info.type = 'piston';
    bp_info.a = 4e-3;  % aperture diameter [m]
    bp_info.k = 2*pi*bp_info.freq/bp_info.c;  % wavenumber
elseif strcmp(bp_type,'gaussian')
    bp_info.type = 'gaussian';
    bp_info.mu = 0;  % mean
    bp_info.sigma = 1;  % variance
end


E_max_all = cell(length(bat_proc_file_all),1);
ar_all = cell(length(bat_proc_file_all),1);

for iF = 1:length(bat_proc_file_all)

    bat_proc_file = bat_proc_file_all(iF).name;
    fprintf('File: %s\n',bat_proc_file);
    
    data = load(fullfile(base_path,bat_proc_path,bat_proc_file));
    mic_num = 1:data.mic_data.num_ch_in_file;
    good_call_idx = find(data.proc.chk_good_call);
    
    for iC = good_call_idx'
 
        iC_save = find(iC==good_call_idx);

        % Get az/el from measurement
        [call_dB,az,el,ch_include_idx] = get_call_azel_dB_data(data,freq_wanted,iC);
        [~,mmidx] = max(call_dB);
        call_max_azel = [az(mmidx),el(mmidx)];
        
        % Model mic output
        mic_dB = model_beam(bp_info,call_max_azel,[az el]);
        
        % Interpolation
        [vq,vq_norm,azq,elq] = interp_bp(az(:),el(:),mic_dB,'rbf');
        
        % Convert from rad to degree
        az = az/pi*180;
        el = el/pi*180;
        azq = azq/pi*180;
        elq = elq/pi*180;
        
        % Set map projection
        mstruct = defaultm('ortho');
        mstruct = defaultm(mstruct);

        % Rotate measurements to use max position as origin
        [mm,mmidx] = max(vq_norm(:));
        origin_max = [elq(mmidx),azq(mmidx)];  % [Lat Lon]
        [el_max,az_max] = rotatem(el,az,origin_max,'forward','degrees');
        [elq_max,azq_max] = rotatem(elq,azq,origin_max,'forward','degrees');

        % Project lat-lon to map projection distance
        [xq,yq] = mfwdtran(mstruct,elq,azq);
        [x_max,y_max] = mfwdtran(mstruct,el_max,az_max);  % map projection
        [xq_max,yq_max] = mfwdtran(mstruct,elq_max,azq_max);  % map projection

        % Plot raw data
        fig_elp = figure;
        set(fig_elp,'position',[150 150 800 400]);
        subplot(121)
        contour(xq,yq,vq_norm,0:-3:-39,'fill','on');
        axis equal; axis([-1.1 1.1 -1.1 1.1]); grid on
        title('Raw data');
        
        % Plot measurements rotated to max point and fit ellipse
        E_max = plot_bp_fit_ellipse(subplot(122),xq_max,yq_max,vq_norm);
        title('Best-fitting ellipse');

        suptitle(sprintf('File #%02d call #%02d, ar=%2.2f',iF,iC,E_max.ar));
        save_fname = sprintf('file%d_call%d.png',iF,iC);
        saveSameSize(fig_elp,'file',fullfile(save_path,save_fname),'format','png','renderer','painters');
%         pause(0.5)
        close(fig_elp)
        
        E_max_all{iF}(iC_save) = E_max;
        ar_all{iF}(iC_save) = E_max.ar;
              
    end
    
end

% Merge all data
ar_all_merge = reshape(cell2mat(ar_all'),1,[]);
