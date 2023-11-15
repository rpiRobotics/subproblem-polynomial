function kin = define_puma_762()

Pi = sym(pi);

alpha_vec = [Pi/2 0 Pi/2 Pi/2 Pi/2 0];
a_vec     = [0 65/100 0 0 0 0];
d_vec     = [0 0 -19/100 6/10 0 0];

kin = dh_to_kin(alpha_vec, a_vec, d_vec);

end