function good_call = isgoodcall(az,el,mic_dB,angle_range,mic_area_frac)
% Check if this call fits the selection criteria
% All angles in *degree*

[maxref,maxref_idx] = max(mic_dB);
bnd_idx = boundary(az,el,0);  % outer boundary of all measured points

% more than 10 mic within az/el range around peak mic
azel_in_range = az<az(maxref_idx)+angle_range(2) & az>az(maxref_idx)+angle_range(1) &...
                el<el(maxref_idx)+angle_range(2) & el>el(maxref_idx)+angle_range(1);
angle_flag = sum(azel_in_range)>=10;

% highest amp within az/el range
max_angle_flag = az(maxref_idx)<angle_range(2) & az(maxref_idx)>angle_range(1) &...
                 el(maxref_idx)<angle_range(2) & el(maxref_idx)>angle_range(1);

% highest amp not on edge
max_edge_flag = ~ismember(maxref_idx,bnd_idx);

% at least 3 mics within 3dB contour
larger_3db_idx = (mic_dB-maxref)>-3.5;
ch_num_3db_flag = sum(larger_3db_idx)>=3;

% sampled area cover large percentage of az-el -90~90 plane
bnd_k = boundary(az,el,0);
bnd = [az(bnd_k),el(bnd_k)];
azel_box_xy = [-90 90];
azel_box = [azel_box_xy([1 1 2 2 1])',azel_box_xy([1 2 2 1 1])'];
[azx,elx]=polyxpoly(az(bnd_k),el(bnd_k),azel_box(:,1),azel_box(:,2));
[in,on] = inpolygon(bnd(:,1),bnd(:,2),azel_box(:,1),azel_box(:,2));
in = in|on;
joint_poly = [[azx,elx];bnd(in,:)];
bnd_joint_k = boundary(joint_poly(:,1),joint_poly(:,2),0);
joint_poly = joint_poly(bnd_joint_k,:);
frac_area = polyarea(joint_poly(:,1),joint_poly(:,2))/(180*180);
mic_area_flag = frac_area>mic_area_frac;

% Good call?
good_call = angle_flag & max_angle_flag & max_edge_flag & ch_num_3db_flag & mic_area_flag;

