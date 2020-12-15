function PlotEllipticity(pos, r, theta)
%
pos = pos(:);
r = r(:);

rot = [[cos(theta), -sin(theta)]; [sin(theta), cos(theta)]];
e = repmat(pos, 1, 2) + rot * diag(r);

hold on;
plot([pos(1) e(1, 1)], [pos(2) e(2, 1)], 'r', 'LineWidth', 2);
plot([pos(1) e(1, 2)], [pos(2) e(2, 2)], 'g', 'LineWidth', 2);
