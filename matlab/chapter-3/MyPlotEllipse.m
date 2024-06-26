function [] = MyPlotEllipse(V, centered_at, rx, ry, set_xlim, set_ylim, fPath, fName)
    % rx and ry are length of major and minor axis.
    cx = centered_at(1,1);
    cy = centered_at(2,1);

    % Set aside the major eigenvector for computing the angle.
    e1 = V(:,1);
    % If all values in major eigenvector are < 1, make them all positive.
    % This is because [1;0] is used for computing angle to rotate by.
    if all(V(:,1) < 0)
        e1 = -e1;
    end
    theta = acos(e1'*[1 ; 0]);

    % If the major axis is in quadrant 4, assign negative angle.
    if  (V(1,1) > 0) && (V(2,1) < 0)
        theta = -theta;
    end
    
    % Create data between 0 and 2*pi.
    t = linspace(0, 2*pi, 200);
    % Parametric values for an ellipse.
    data = [rx*cos(t); ry*sin(t)];
    % The rotation by angle theta (in radians).
    rotation = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    ellipse_values = rotation *data + centered_at;

    hold on
        % Plot the ellipse.
        p = plot(ellipse_values(1,:), ellipse_values(2,:), ...
            '-', 'Color', 'blue');
        
        % Plot the vector for the major axis.
        quiver(cx, cy, V(1,1), V(2,1), rx, 'color', 'k');
        quiver(cx, cy, -V(1,1), -V(2,1), rx, 'color', 'k');
        
        % Plot the vector for the minor axis.
        quiver(cx, cy, V(1,2), V(2,2), ry, 'color', 'k');
        quiver(cx, cy, -V(1,2), -V(2,2), ry, 'color', 'k');
        
        % Plot red point for major axis.
        plot(cx+rx*V(1,1), cy+rx*V(2,1), 'o', 'Color', 'red');
        plot(cx-rx*V(1,1), cy-rx*V(2,1), 'o', 'Color', 'red');

        % Plot green point for minor axis.
        plot(cx+ry*V(1,2), cy+ry*V(2,2), 'o', 'Color', 'green');
        plot(cx-ry*V(1,2), cy-ry*V(2,2), 'o', 'Color', 'green');

        grid on
        xlim(set_xlim)
        ylim(set_ylim)
        axis square
    hold off
    if (fPath(end) ~= '\')
        fPath = append(fPath,'\');
    end
    saved_file = append(fPath, fName, '.png');
    saveas(p, saved_file, 'png')
end