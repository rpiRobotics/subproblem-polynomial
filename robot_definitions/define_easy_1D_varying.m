function kin = define_easy_1D_varying()
% R_06 = eye(3);
% p_06 = [1;1;1]*2;

ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];
zv = [0;0;0];

kin.joint_type = [0 0 0 0 0 0];
kin.P = [zv 4*ex 3*ex 2*ex ex zv zv];
kin.H = [ez ey ez ey ez ey];
end