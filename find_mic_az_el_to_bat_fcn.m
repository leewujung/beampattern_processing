function [mic2bat_2d,mic2bat_x] = find_mic_az_el_to_bat_fcn(mic_to_bat_vec,aim_v,norm_v)
% Find azimuth and elevation angle of mics with respective to the bat heam aim
%
% INPUT
%   aim_v    bat head aim
%   norm_v   normal vec to the bat head
%   mic_to_bat_vec  direction vector of mics to bat
% OUTPUT
%   mic_to_bat_angle   [azimuth,elevation] angle of the mics to bat [radian]
%
% Wu-Jung Lee
% 2014/10/12  
% 2015/10/13  use matrix rotation to incorporate bat head aim and roll angle

% First rotation, from global to bat frame
R_prox = [1,0,0;...
          0,1,0;...
          0,0,1]';
R_dist = [cross(aim_v,norm_v);...
          aim_v;...
          norm_v]';
T = R_dist'*R_prox;
mic_to_bat_vec_rot = (T*mic_to_bat_vec')';

% Second rotation, from bat aiming at +Y to bat aiming at +X
% (for convenience of calculating azimuth angles)
R_dist = [0,1,0;...
          -1,0,0;...
          0,0,1]';
T = R_dist'*R_prox;
mic_to_bat_vec_rot = (T*mic_to_bat_vec_rot')';

% Convert to spherical coord: used for 2D plotting
[mic2bat_2d(:,1),mic2bat_2d(:,2),mic2bat_2d(:,3)] =...
    cart2sph(mic_to_bat_vec_rot(:,1),mic_to_bat_vec_rot(:,2),mic_to_bat_vec_rot(:,3));
mic2bat_2d(:,1) = -mic2bat_2d(:,1);

% Compensate for over 90deg elevation angle: used for cross plotting
mic2bat_x = mic2bat_2d;

adj_idx = mic2bat_x(:,1)>pi/2 | mic2bat_x(:,1)<-pi/2;  % if fall into the back half
mic2bat_x(adj_idx,2) = pi-mic2bat_x(adj_idx,2);  % compensate for the elevation
% mic_to_bat_angle(adj_idx,1) = 2*pi+mic_to_bat_angle(adj_idx,1);  % change from -180:0 to 180:360

adj_idx = mic2bat_x(:,2)>pi;
mic2bat_x(adj_idx,2) = -(2*pi-mic2bat_x(adj_idx,2));

mic2bat_x(:,1) = -mic2bat_x(:,1);  % flip left-right

%% FOUND OUT THE AZIMUTH ANGLE SHOULD BE ROTATED, NOT MIRRORED
% adj_idx = mic_to_bat_angle(:,1)<0;
% mic_to_bat_angle(adj_idx,2) = pi-mic_to_bat_angle(adj_idx,2);
% mic_to_bat_angle(adj_idx,1) = -mic_to_bat_angle(adj_idx,1);  % mirror image [0:-180] to [0:180]
% 
% adj_idx = mic_to_bat_angle(:,2)>pi;
% mic_to_bat_angle(adj_idx,2) = -(2*pi-mic_to_bat_angle(adj_idx,2));
% 
% mic_to_bat_angle(:,1) = -(mic_to_bat_angle(:,1)-pi/2);  % shift by 90 deg and flip left-right

%% THE ADJUSTMENT BELOW ARE EQUIVALENT TO THOSE USED ABOVE
% adj_idx = mic_to_bat_angle(:,1)<0 & mic_to_bat_angle(:,2)>0;  % for pts on the back and anove horizon
% mic_to_bat_angle(adj_idx,2) = pi-mic_to_bat_angle(adj_idx,2);
% 
% adj_idx = mic_to_bat_angle(:,1)<0 & mic_to_bat_angle(:,2)<0;  % for pts on the back and below horizon
% mic_to_bat_angle(adj_idx,2) = -(pi-(-mic_to_bat_angle(adj_idx,2)));
% 
% adj_idx = mic_to_bat_angle(:,1)<0;
% mic_to_bat_angle(adj_idx,1) = -mic_to_bat_angle(adj_idx,1);  % mirror image [0:-180] to [0:180]

