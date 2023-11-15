function kin = define_raghavan_roth()

Pi = sym(pi);

alpha_vec = [20*Pi/180 31*Pi/180 45*Pi/180 81*Pi/180 12*Pi/180 100*Pi/180];
a_vec     = [8/10 12/10 33/100 18/10 6/10 22/10];
d_vec     = [9/10 37/10 10/10 5/10 21/10 63/100];

kin = dh_to_kin(alpha_vec, a_vec, d_vec);

end