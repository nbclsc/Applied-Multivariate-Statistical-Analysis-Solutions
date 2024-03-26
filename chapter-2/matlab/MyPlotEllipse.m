function [] = MyPlotEllipse(V, D, c, fName)
    % Compute points corresponding to axis-oriented ellipse
    xc = 0;
    yc = 0;
    b = c/sqrt(D(1,1));
    a = c/sqrt(D(2,2));
    theta = acos(-V(:,1)'*[1 ; 0]); % acos(1/sqrt(3));
    
    t = linspace(0, 2*pi, 200);
    xt = b * cos(t) + xc;
    yt = a * sin(t) + yc;
    
    % Apply rotation by angle theta (in radians).
    cot = cos(theta); sit = sin(theta);
    x = xt * cot - yt * sit;
    y = xt * sit + yt * cot;

    hold on
        % Plot the ellipse.
        p=plot(x, y, '-', 'Color', 'blue');
        
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
        title(append('c^2 = ', num2str(c)))
        grid on
        pbaspect([1 1 1])
        % legend('asd', 'asda')
    hold off
    saved_file = append('C:\Users\nbclsc\Desktop\applied-multivariate-statistics\solutions\chapter-2\', fName, '.png');
    saveas(p, saved_file, 'png')
end