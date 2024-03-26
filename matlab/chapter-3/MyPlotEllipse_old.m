function [] = MyPlotEllipse(V, D, c, fPath, fName)
    % The length of major and minor axis.
    rx = c/sqrt(D(1,1));
    ry = c/sqrt(D(2,2));

    % Set aside the major eigenvector for computing the angle.
    e1 = V(:,1);
    % If all values in major eigenvector are < 1, make them all positive.
    % This is because [1;0] is used for computing angle to rotate by.
    if all(V(:,1) < 0)
        e1 = -e1;
    end
    theta = acos(e1'*[1 ; 0]);
    
    % Create data between 0 and 2*pi.
    t = linspace(0, 2*pi, 200);

    % Parametric values for an ellipse.
    data = [rx*cos(t); ry*sin(t)];

    % The rotation by angle theta (in radians).
    rotation = [cos(theta) -sin(theta); sin(theta) cos(theta)];

    % Where to center the ellipse.
    centered_at = [0; 0];
    ellipse_values = rotation *data + centered_at;

    hold on
        % Plot the ellipse.
        p = plot(ellipse_values(1,:), ellipse_values(2,:), ...
            '-', 'Color', 'blue');
        
        % Plot the vector for the major axis.
        quiver(0, 0, V(1,1), V(2,1), c, 'color', 'k');
        quiver(0, 0, -V(1,1), -V(2,1), c, 'color', 'k');
        
        % Plot the vector for the minor axis.
        quiver(0, 0, V(1,2), V(2,2), c, 'color', 'k');
        quiver(0, 0, -V(1,2), -V(2,2), c, 'color', 'k');
        
        % Plot red point for major axis.
        plot((c/sqrt(D(1,1)))*V(1,1), (c/sqrt(D(1,1)))*V(2,1), 'o', 'Color', 'red');
        plot(-(c/sqrt(D(1,1)))*V(1,1), -(c/sqrt(D(1,1)))*V(2,1), 'o', 'Color', 'red');

        % Plot green point for minor axis.
        plot((c/sqrt(D(2,2)))*V(1,2), (c/sqrt(D(2,2)))*V(2,2), 'o', 'Color', 'green');
        plot(-(c/sqrt(D(2,2)))*V(1,2), -(c/sqrt(D(2,2)))*V(2,2), 'o', 'Color', 'green');

        % % Plot red point for major axis.
        % plot(rx*V(1,1), rx*V(2,1), 'o', 'Color', 'red');
        % plot(-rx*V(1,1), -rx*V(2,1), 'o', 'Color', 'red');
        % 
        % % Plot green point for minor axis.
        % plot(ry*V(1,2), ry*V(2,2), 'o', 'Color', 'green');
        % plot(-ry*V(1,2), -ry*V(2,2), 'o', 'Color', 'green');
        title(append('c^2 = ', num2str(c)))
        grid on
        % pbaspect([1 1 1])
        axis square
        % legend('asd', 'asda')
    hold off
    if (fPath(end) ~= '\')
        fPath = append(fPath,'\');
    end
    saved_file = append(fPath, fName, '.png');
    saveas(p, saved_file, 'png')
end