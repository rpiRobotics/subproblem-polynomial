function kin = define_manseur_doty()

Pi = sym(pi);

alpha_vec = [Pi/2 0 Pi/2 0 Pi/2 0];
a_vec     = [3/10 1 0 15/10 0 0];
d_vec     = [0 0 2/10 0 0 0];

kin = dh_to_kin(alpha_vec, a_vec, d_vec);

end