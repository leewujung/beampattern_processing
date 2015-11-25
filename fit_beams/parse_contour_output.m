function Cout = parse_contour_output(C)
% Parse the output obtained from fcn contour
% Code adapted from contourcs.m written by Kesh Ikuma on Matlab filexchange

K = 0;
n0 = 1;
while n0<=size(C,2)
   K = K + 1;
   n0 = n0 + C(2,n0) + 1;
end
% initialize output struct
el = cell(K,1);
Cout = struct('Level',el,'Length',el,'X',el,'Y',el);

% fill the output struct
n0 = 1;
for k = 1:K
   Cout(k).Level = C(1,n0);
   idx = (n0+1):(n0+C(2,n0));
   Cout(k).Length = C(2,n0);
   Cout(k).X = C(1,idx);
   Cout(k).Y = C(2,idx);
   n0 = idx(end) + 1; % next starting index
end