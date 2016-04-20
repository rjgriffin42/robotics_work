function draw_foot_2d(p_foot, linetype)
  if nargin == 1
    linetype = '-';
  end

  % configure appearance 
  line_color = [0 0 0];
  corners = [-0.125, -0.0675, 0; ...
             0.125, -0.0675, 0; ...
             0.125, 0.0675, 0; ...
            -0.125, 0.0675, 0];

  % compute corner locations
  for i = 1:4
    corner = makehgtform('zrotate', p_foot(4))* ...
             makehgtform('yrotate', p_foot(5))* ...
             makehgtform('zrotate', p_foot(6))* ...
             [corners(i, :) 1]';
    corners(i, 1) = corner(1)/corner(4) + p_foot(1);
    corners(i, 2) = corner(2)/corner(4) + p_foot(2);
    corners(i, 3) = corner(3)/corner(4) + p_foot(3);
  end

  % plot foot
  hold on;
  corners(5, :) = corners(1, :);
  for i = 1:4
    plot([corners(i, 1) corners(i + 1, 1)], ...
         [corners(i, 2) corners(i + 1, 2)], ...
         'Color', line_color, 'LineWidth', 1, 'LineStyle', linetype);
  end
end
