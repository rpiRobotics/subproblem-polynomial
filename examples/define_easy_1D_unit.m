function kin = define_easy_1D_unit()
% R_06 = eye(3);
% p_06 = [1;1;1]*1.25;

zv = [0;0;0];
ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];


kin.joint_type = [0 0 0 0 0 0];
kin.P = [zv ex ex ex ex zv zv];
kin.H = [ez ey ez ey ez ey];
end