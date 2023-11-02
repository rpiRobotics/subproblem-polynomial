function kin = define_RRC_fixed_q6()
    kin = get_kin_partial();
end

function kin = get_kin()
    zv = [0;0;0];
    ex = [1;0;0];
    ey = [0;1;0];
    ez = [0;0;1];
    
    % Define in inches
    p01 = zv;
    p12 = 20*ex -4*ey;
    p23 = 4*ey;
    p34 = 21.5*ex +3.375*ey;
    p45 = -3.375*ey;
    p56 = 21.5*ex+3.325*ey;
    p67 = -3.325*ey;
    p7T = 7*ex;
    
    kin.P = [p01 p12 p23 p34 p45 p56 p67 p7T];
    kin.H = [ex ez ex ez ex ez ex];
    kin.joint_type=[0 0 0 0 0 0 0];
end

function [kin_partial, R_6T] = get_kin_partial()
    kin = get_kin();
    [kin_partial, R_6T] = fwdkin_partial(kin,sym(pi)/6, 6);

    % Move O_5 along h_5 and move O_6 along h_6 s.t. p_56 = 0
    % p_56 = alpha_1 * h_5 + alpha_2 * h_6
    zv = [0;0;0];
    alpha = pinv(kin_partial.H(:,5:6))*kin_partial.P(:,6);
    delta_p_45 = alpha(1)*kin_partial.H(:,5);
    delta_p_6T = alpha(2)*kin_partial.H(:,6);
    kin_partial.P(:,[5 6 7]) = [kin_partial.P(:,5)+delta_p_45 zv kin_partial.P(:,7)+delta_p_6T];
end