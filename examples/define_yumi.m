function kin = define_yumi()

% alpha_vec = 2*pi_sym/360*[-90 90 -90 -90 -90 90 0];
% a_vec     = [-30 30 40.5 40.5 27 -27 0];
% d_vec     = [166 0 251.5 0 265 0 36];
% 
% kin_full = dh_to_kin(alpha_vec, a_vec, d_vec);
% q_0 = [0 0 0 -90 180 0 0]

ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];
zv = [0;0;0];

kin.joint_type = zeros([1 7]);
kin.H = [ez ey ez ey ex ey ex];
kin.P = [zv 166*ez-30*ex 30*ex 251.5*ez+40.5*ex 40.5*ez 265*ex-27*ez 27*ez zv];
end