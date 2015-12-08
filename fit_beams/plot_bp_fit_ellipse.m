function E = plot_bp_fit_ellipse(haxes,xq,yq,vq)
% Plot beampattern and best-fitting ellipse

axes(haxes)
[C,~] = contour(xq,yq,vq,0:-3:-39,'fill','on');
Cout = parse_contour_output(C);
c3db_xy = [];
for iT=1:length(Cout)  % in case contour break into pieces
    if Cout(iT).Level == -3
        c3db_xy = [c3db_xy; Cout(iT).X',Cout(iT).Y'];
    end
end
A = EllipseDirectFit(c3db_xy);  % fit ellipse (direct fit)
E = get_ellipse_param(A);       % get ellipse parameters
xmin = min(c3db_xy(:,1))-range(c3db_xy(:,1)*0.5);
xmax = max(c3db_xy(:,1))+range(c3db_xy(:,1)*0.5);
ymin = min(c3db_xy(:,2))-range(c3db_xy(:,1)*0.5);
ymax = max(c3db_xy(:,2))+range(c3db_xy(:,1)*0.5);
hold on
fit_df = ezplot(E.eqt,[xmin,xmax,ymin,ymax]);  % plot ellipse
set(fit_df,'linecolor','b','linewidth',2);
plot(E.x0,E.y0,'r*');
text(E.x0,E.y0,sprintf('%2.3f, %2.3f',E.x0,E.y0))
title('')
hold off
axis equal; grid on
axis([-1.1 1.1 -1.1 1.1]); 