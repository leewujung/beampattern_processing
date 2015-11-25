% 2015 11 24  Use real measurement az/el angles for model beampattern

clear
usrn = getenv('username');
addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\EllipseDirectFit']);
addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\beampattern_processing']);

% Bat data path
save_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_model'];
base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing'];
bat_proc_path = './proc_output';
bat_proc_file = 'rousettus_20150825_36134_02_mic_data_bp_proc';
freq_wanted = 35e3;
goodcall_angle_range = [-60 60];
mic_area_frac = 0.4;

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

% iF =1;
data = load(fullfile(base_path,bat_proc_path,bat_proc_file));
mic_num = 1:data.mic_data.num_ch_in_file;
good_call_idx = find(data.proc.chk_good_call);

for iC = good_call_idx'
    
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

    % Rotate measurements to max position
    [mm,mmidx] = max(vq_norm(:));
    origin = [elq(mmidx),azq(mmidx)];  % [Lat Lon]
    [el_rot,az_rot] = rotatem(el,az,origin,'forward','degrees');
    [elq_rot,azq_rot] = rotatem(elq,azq,origin,'forward','degrees');
    
    % Project lat-lon to map projection distance
    mstruct = defaultm('ortho');
    mstruct = defaultm(mstruct);
    [x,y] = mfwdtran(mstruct,el,az);
    [x_rot,y_rot] = mfwdtran(mstruct,el_rot,az_rot);
    [xq,yq] = mfwdtran(mstruct,elq,azq);
    [xq_rot,yq_rot] = mfwdtran(mstruct,elq_rot,azq_rot);
    
    % Find -3dB contour
    figure
    [C,h] = contour(xq_rot,yq_rot,vq_norm,0:-3:-9,'fill','on');
    Cout = parse_contour_output(C);
    c3db_xy = [];
    for iT=1:length(Cout)  % in case contour break into pieces
        if Cout(iT).Level == -3
            c3db_xy = [c3db_xy; Cout(iT).X',Cout(iT).Y'];
        end
    end

    % Fit ellipse (direct fit)
    A = EllipseDirectFit(c3db_xy);
    E = get_ellipse_param(A);
    
    % Plot ellipse (direct fit)
    xmin = min(c3db_xy(:,1))-range(c3db_xy(:,1)*0.5);
    xmax = max(c3db_xy(:,1))+range(c3db_xy(:,1)*0.5);
    ymin = min(c3db_xy(:,2))-range(c3db_xy(:,1)*0.5);
    ymax = max(c3db_xy(:,2))+range(c3db_xy(:,1)*0.5);
    hold on
    plot(x_rot,y_rot,'r+');
    fit_df = ezplot(E.eqt,[xmin,xmax,ymin,ymax]);
    set(fit_df,'linecolor','b','linewidth',2);
    hold off
    axis equal
    axis([-pi pi -pi/2 pi/2]);
    axis equal
    title(sprintf('Call#%02d',iC));

    AR_dirfit(iC) = E.ar;
    fprintf('Call #%02d, AR_dirfit=%2.2f\n',iC,E.ar);
end
AR_dirfit(AR_dirfit<1) = 1./AR_dirfit(AR_dirfit<1);


