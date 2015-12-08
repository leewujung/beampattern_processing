% 2015 11 14  Fit ellipse to all good calls
% 2015 11 23  - Implement map rotation with origin at peak mic location
%             - Fit ellipse using orthogonal projected distance
% 2015 12 08  Use merge_call code to do the fitting and process all data

clear
warning off

usrn = getenv('username');
base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing'];
save_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing\20151208_good_call_fit_ellipse'];
if ~exist(save_path,'dir')
    mkdir(save_path);
end
bat_proc_path = './proc_output_checked';
bat_proc_file_all = dir(fullfile(base_path,bat_proc_path,'rousettus_20150825_*'));
% bat_proc_file = 'rousettus_20150825_36134_05_mic_data_bp_proc';
plot_file_opt = 0;
freq_wanted = 35e3;

addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\EllipseDirectFit']);
addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);


E_max_all = cell(length(bat_proc_file_all),1);
ar_all = cell(length(bat_proc_file_all),1);
click_side_all = cell(length(bat_proc_file_all),1);

for iF = 1:length(bat_proc_file_all)

    bat_proc_file = bat_proc_file_all(iF).name;
    fprintf('File: %s\n',bat_proc_file);
    
    data = load(fullfile(base_path,bat_proc_path,bat_proc_file));
    mic_num = 1:data.mic_data.num_ch_in_file;
    good_call_idx = find(data.proc.chk_good_call);
    
    if plot_file_opt
        fig_all = figure;
    end
    numrow = ceil(length(good_call_idx)/2);
    
    for iC = good_call_idx'

        iC_save = find(iC==good_call_idx);
        
        % Get call info
        [call_dB,az,el,ch_include_idx] = get_call_azel_dB_data(data,freq_wanted,iC);
        [vq,vq_norm,azq,elq] = interp_bp(az(ch_include_idx),el(ch_include_idx),call_dB(ch_include_idx),'rbf');
        call_dB_norm = call_dB-max(call_dB);
        az = az/pi*180;  % convert to [deg]
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
        
        % Determine right/left click
        if azq(mmidx)<0
            click_side = 0;  % left: click_side = 0
            click_side_t = 'left';
        else
            click_side = 1;  % right: click_side = 1
            click_side_t = 'right';
        end

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

        suptitle(sprintf('File #%02d call #%02d, ar=%2.2f, side=%s',iF,iC,E_max.ar,click_side_t));
        save_fname = sprintf('file%d_call%d.png',iF,iC);
        saveSameSize(fig_elp,'file',fullfile(save_path,save_fname),'format','png','renderer','painters');
%         pause(0.5)
        close(fig_elp)
        
        E_max_all{iF}(iC_save) = E_max;
        ar_all{iF}(iC_save) = E_max.ar;
        click_side_all{iF}(iC_save) = click_side;
        
        
        % Plot rotated bp using ellipse center
        if plot_file_opt
            figure(fig_all)
            subplot(2,numrow,find(iC==good_call_idx));
            [C,~] = contour(xq_rot_ecen,yq_rot_ecen,vq_norm,0:-3:-39,'fill','on');
            title(sprintf('Call #%02d',iC));
        end
    end
    
    if plot_file_opt
        figure(fig_all)
        suptitle(bat_proc_file);
    end
    
end
warning on


% Merge all data
click_side_merge = reshape(cell2mat(click_side_all'),1,[]);
ar_all_merge = reshape(cell2mat(ar_all'),1,[]);

