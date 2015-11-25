function [vq,vq_norm,azq,elq] = interp_bp(az,el,mic_dB,method)
% az,el    azimuth and elevation location for each mic [radian]
% mic_dB   corresponding mic recording amplitude [dB]
% method   interpolation method, 'natural' or 'rbf'

% Get data
az = az(:);
el = el(:);
mic_dB = mic_dB(:);

% Interpolation
maxref = max(mic_dB);
[azq,elq] = meshgrid(min(az):pi/180:max(az),min(el):pi/180:max(el));
if strcmp(method,'natural')  % natural neighbor interpolation
    vq = griddata(az,el,mic_dB,azq,elq,'natural');
elseif strcmp(method,'rbf')  % radial basis function interpolation
    vq = rbfinterp([azq(:)';elq(:)'],rbfcreate([az';el'],mic_dB','RBFFunction','multiquadrics'));
    vq = reshape(vq,size(azq));
end
vq_norm = vq-maxref;

% Find boundary
k = boundary(az,el,0);  % outer boundary of all measured points
[in,on] = inpolygon(azq,elq,az(k),el(k));
in = in|on;

% Set values outside of az-el boundary to NaN
vq(~in) = NaN;
vq_norm(~in) = NaN;

