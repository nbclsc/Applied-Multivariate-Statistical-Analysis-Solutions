function [] = my_plot_ellipse(V, centered_at, rx, ry, set_xlim, set_ylim)
    % Doesn't matter order of eigenvectors. Rotation of the points based on
    % the first eigenvector. The largest eigenvalue will be the major axis
    % of the ellipse. The smallest eigenvalue will be the minor axis of the
    % ellipse.
    % V: Eigenvectors.
    % centered_at: Where to center the ellipse. xbar vector.
    % rx: The length in the direction of first eigenvector.
    % ry: The length in the direction of second eigenvector.
    % set_xlim: Vector of plot limits in the x-direction.
    % set_ylim: Vector of plot limits in the y-direction.
    % fPath: Path to where to save the ellipse file.
    % fName: File name to save as.

    % rx and ry are length of axis (vectors) in the ellipse.
    cx = centered_at(1,1);
    cy = centered_at(2,1);

    % Set aside the first eigenvector for computing the angle.
    e1 = V(:,1);
    % If all values in first eigenvector are < 1, make them all positive.
    % This is because [1;0] is used for computing angle to rotate by.
    if all(V(:,1) < 0)
        e1 = -e1;
    end
    theta = acos(e1'*[1 ; 0]);

    % If the first eigenvector is in quadrant 4, assign negative angle.
    if  (V(1,1) > 0) && (V(2,1) < 0)
        theta = -theta;
    end
    
    % Create data between 0 and 2*pi.
    t = linspace(0, 2*pi, 200);
    % Parametric values for an ellipse.
    data = [rx*cos(t); ry*sin(t)];
    % The rotation by angle theta (in radians) to first eigenvector.
    rotation = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    fprintf('Rotate to first eigenvector by %5.3f degrees.\n', rad2deg(theta));
    ellipse_values = rotation *data + centered_at;

    hold on
        % Plot the rotated ellipse.
        p = plot(ellipse_values(1,:), ellipse_values(2,:), ...
            '-', 'Color', 'blue');
        
        % Plot the first vector.
        quiver(cx, cy, V(1,1), V(2,1), rx, 'color', 'k');
        quiver(cx, cy, -V(1,1), -V(2,1), rx, 'color', 'k');
        
        % Plot the second vector.
        quiver(cx, cy, V(1,2), V(2,2), ry, 'color', 'k');
        quiver(cx, cy, -V(1,2), -V(2,2), ry, 'color', 'k');
        
        % Plot red point for one of the axis.
        plot(cx+rx*V(1,1), cy+rx*V(2,1), 'o', 'Color', 'red');
        plot(cx-rx*V(1,1), cy-rx*V(2,1), 'o', 'Color', 'red');

        % Plot green point for the other axis.
        plot(cx+ry*V(1,2), cy+ry*V(2,2), 'o', 'Color', 'green');
        plot(cx-ry*V(1,2), cy-ry*V(2,2), 'o', 'Color', 'green');

        grid on
        xlim(set_xlim)
        ylim(set_ylim)
        axis square
    hold off
end