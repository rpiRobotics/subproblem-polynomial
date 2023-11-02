function kin = define_RRC_fixed_q7()
zv = [0;0;0];
ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];

% Define in inches,
p01 = zv;
p12 = 20*ex -4*ey;
p23 = 4*ey;
p34 = sym(21.5)*ex + sym(3.375)*ey;
p45 = -sym(3.375)*ey;
p56 = sym(21.5)*ex+sym(3.325)*ey;
p67 = zv;

kin.joint_type=[0 0 0 0 0 0];
kin.H = [ex ez ex ez ex ez];
kin.P = [p01 p12 p23 p34 p45 p56 p67];

end