function [allwp, e] = waypoints(states, count, p)
    allwp = zeros(12,1);
    wp = zeros(3,1);
    xi = states(1:3);
    if count<9
        wp = p.wp(:,uint8(count+1));
    else
        wp(1) = 0;
        wp(2) = 0;
        wp(3) = 0;
    end
    e = norm(xi-wp);
    allwp(1:3) = wp;
    allwp(4:12) = 0;
end