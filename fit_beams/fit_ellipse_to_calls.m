% 2015 11 14  Fit ellipse to all good calls
% 2015 11 23  - Implement map rotation with origin at peak mic location
%             - Fit ellipse using orthogonal projected distance

clear
warning off

usrn = getenv('username');
base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing'];
bat_proc_path = './proc_output';
bat_proc_file = 'rousettus_20150825_36134_02_mic_data_bp_proc';

addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\fitellipse']);
addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\EllipseDirectFit']);
addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\beampattern_processing']);

data = load(fullfile(base_path,bat_proc_path,bat_proc_file));
freq_wanted = 35e3;
mic_num = 1:data.mic_data.num_ch_in_file;
good_call_idx = find(data.proc.chk_good_call);

AR_dirfit = zeros(length(data.proc.chk_good_call),1);
AR_nonlin = zeros(length(data.proc.chk_good_call),1);
for iC = good_call_idx'

    [call_dB,az,el,ch_include_idx] = get_call_azel_dB_data(data,freq_wanted,iC);
    [vq,vq_norm,azq,elq] = interp_bp(az(ch_include_idx),el(ch_include_idx),call_dB(ch_include_idx),'rbf');
    az = az/pi*180;
    el = el/pi*180;
    azq = azq/pi*180;
    elq = elq/pi*180;
    
    % Rotate measurements to max position
    [mm,mmidx] = max(vq_norm(:));
    origin = [elq(mmidx),azq(mmidx)];  % [Lat Lon]
    [elq_rot,azq_rot] = rotatem(elq,azq,origin,'forward','degrees');
    
    % Project lat-lon to map projection distance
    mstruct = defaultm('ortho');
    mstruct = defaultm(mstruct);
    [xq,yq] = mfwdtran(mstruct,elq,azq);
    [xq_rot,yq_rot] = mfwdtran(mstruct,elq_rot,azq_rot);
    
    % Find -3dB contour
    figure
    [C,h] = contour(xq_rot,yq_rot,vq_norm,0:-3:-39,'fill','on');
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
    fit_df = ezplot(E.eqt,[xmin,xmax,ymin,ymax]);
    set(fit_df,'linecolor','b','linewidth',2);
    hold off
    
    AR_dirfit(iC) = E.ar;
    fprintf('Call #%02d, AR_dirfit=%2.2f\n',iC,E.ar);
    
    % Nonlinear least square
    try
        [zg, ag, bg, alphag] = fitellipse(c3db_xy);
        AR_nonlin(iC) = ag/bg;
        hold on
        fit_nonlin = plotellipse(zg,ag,bg,alphag,'r');
        hold off
        fprintf('Call #%02d, AR_nonlin=%2.2f\n',iC,ag/bg);
    catch
        fprintf('Call #%02d cannot be fitted with ellipse\n',iC);
    end
    title(sprintf('Call#%02d',iC));
    
    axis([-pi pi -pi/2 pi/2]);
    axis equal
    
end
warning on

AR_dirfit(AR_dirfit<1) = 1./AR_dirfit(AR_dirfit<1);
AR_nonlin(AR_nonlin<1) = 1./AR_nonlin(AR_nonlin<1);

% % Plot using geoshow
% elqm = elq;
% elqm(isnan(vq_norm)) = NaN;
% azqm = azq;
% azqm(isnan(vq_norm)) = NaN;
% 
% figure;
% axesm eckert4;
% framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
% gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
% axis off
% geoshow(elqm,azqm,vq_norm,'displaytype','texturemap');
% contourm(elq,azq,vq_norm,-3,'w','linewidth',2);
% textm(el,az,num2str(mic_num(ch_include_idx)'),'horizontalalignment','center');


