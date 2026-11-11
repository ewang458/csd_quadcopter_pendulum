function R = rotation(angle, xyz)
    switch lower(xyz)
        case 'x'
            R = [1, 0, 0;
                0, cos(angle), -sin(angle); 
                0, sin(angle), cos(angle)];
        case 'y'
            R = [cos(angle), 0, sin(angle); 
                0, 1, 0;
                -sin(angle), 0, cos(angle)];
        case 'z'
            R = [cos(angle), -sin(angle), 0;
                sin(angle), cos(angle), 0;
                0, 0, 1];
        otherwise
            R = eye(3);
    end
end