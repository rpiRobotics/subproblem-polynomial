function R = half_tan_sp1_R(p1, p2, k)
    alpha_sin = dot(p2, cross(k,p1));
    alpha_cos = -dot(p2, cross(k,cross(k,p1)));
    % alpha = dot(cross(k,p2), cross(k,p2));
    alpha = dot(cross(k,p1), cross(k,p1)); % Seems better
    s = alpha_sin/alpha;
    c = alpha_cos/alpha;
    % R = eye(3,3)+s*hat(k)+(1-c)*hat(k)*hat(k);
    R = k*k' + s*hat(k) - c * hat(k)*hat(k);
end