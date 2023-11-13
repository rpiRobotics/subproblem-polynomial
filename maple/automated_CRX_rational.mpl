# Find all candidates for x_4 for the CRX robot
# Given R_06, p_06


# Input: R_06 and p_06
#b := convert(225/100*degrees, radians);
#g := convert(312/100*degrees, radians);
xb := 13/662;
xg := 32/1175;
t1 := 3141/10000;
t2 := 5926/10000;
t3 := 5358/10000;

# Define the system of equations
P1 := (xb, xg, t1, t2, t3) -> 2*x1 + 4*x2*x4 + 4*x3*x4 - 4*x4*xb + 4*xb*xg + 2*x1*x2^2 + 2*x1*x3^2 - 2*x1*x4^2 - 2*x1*xb^2 + 2*x1*xg^2 - 4*x1^2*x2*x4 - 4*x1^2*x3*x4 - 4*x2*x3^2*x4 - 4*x2^2*x3*x4 - 4*x1^2*x4*xb - 4*x2*x4*xb^2 + 4*x2^2*x4*xb - 4*x3*x4*xb^2 + 4*x3^2*x4*xb + 4*x2*x4*xg^2 + 4*x3*x4*xg^2 - 4*x1^2*xb*xg + 4*x2^2*xb*xg + 4*x3^2*xb*xg + 4*x4*xb*xg^2 - 4*x4^2*xb*xg + 2*x1*x2^2*x3^2 - 2*x1*x2^2*x4^2 - 2*x1*x3^2*x4^2 - 2*x1*x2^2*xb^2 - 2*x1*x3^2*xb^2 + 2*x1*x4^2*xb^2 + 2*x1*x2^2*xg^2 + 2*x1*x3^2*xg^2 - 2*x1*x4^2*xg^2 - 2*x1*xb^2*xg^2 + 16*x2*x3*x4*xb - 2*x1*x2^2*x3^2*x4^2 - 2*x1*x2^2*x3^2*xb^2 + 2*x1*x2^2*x4^2*xb^2 + 2*x1*x3^2*x4^2*xb^2 + 2*x1*x2^2*x3^2*xg^2 - 2*x1*x2^2*x4^2*xg^2 - 2*x1*x3^2*x4^2*xg^2 - 2*x1*x2^2*xb^2*xg^2 - 2*x1*x3^2*xb^2*xg^2 + 2*x1*x4^2*xb^2*xg^2 + 4*x1^2*x2*x3^2*x4 + 4*x1^2*x2^2*x3*x4 + 4*x1^2*x2*x4*xb^2 + 4*x1^2*x2^2*x4*xb + 4*x1^2*x3*x4*xb^2 + 4*x1^2*x3^2*x4*xb + 4*x2*x3^2*x4*xb^2 + 4*x2^2*x3*x4*xb^2 - 4*x2^2*x3^2*x4*xb - 4*x1^2*x2*x4*xg^2 - 4*x1^2*x3*x4*xg^2 - 4*x2*x3^2*x4*xg^2 - 4*x2^2*x3*x4*xg^2 - 4*x1^2*x2^2*xb*xg - 4*x1^2*x3^2*xb*xg + 4*x1^2*x4*xb*xg^2 + 4*x1^2*x4^2*xb*xg + 4*x2^2*x3^2*xb*xg - 4*x2*x4*xb^2*xg^2 - 4*x2^2*x4*xb*xg^2 - 4*x2^2*x4^2*xb*xg - 4*x3*x4*xb^2*xg^2 - 4*x3^2*x4*xb*xg^2 - 4*x3^2*x4^2*xb*xg + 16*x1^2*x2*x3*x4*xb - 16*x2*x3*x4*xb*xg^2 + 2*x1*x2^2*x3^2*x4^2*xb^2 - 2*x1*x2^2*x3^2*x4^2*xg^2 - 2*x1*x2^2*x3^2*xb^2*xg^2 + 2*x1*x2^2*x4^2*xb^2*xg^2 + 2*x1*x3^2*x4^2*xb^2*xg^2 - 32*x1*x2*x4*xb*xg - 32*x1*x3*x4*xb*xg - 4*x1^2*x2*x3^2*x4*xb^2 - 4*x1^2*x2^2*x3*x4*xb^2 - 4*x1^2*x2^2*x3^2*x4*xb + 4*x1^2*x2*x3^2*x4*xg^2 + 4*x1^2*x2^2*x3*x4*xg^2 - 4*x1^2*x2^2*x3^2*xb*xg + 4*x1^2*x2*x4*xb^2*xg^2 - 4*x1^2*x2^2*x4*xb*xg^2 + 4*x1^2*x2^2*x4^2*xb*xg + 4*x1^2*x3*x4*xb^2*xg^2 - 4*x1^2*x3^2*x4*xb*xg^2 + 4*x1^2*x3^2*x4^2*xb*xg + 4*x2*x3^2*x4*xb^2*xg^2 + 4*x2^2*x3*x4*xb^2*xg^2 + 4*x2^2*x3^2*x4*xb*xg^2 - 4*x2^2*x3^2*x4^2*xb*xg + 2*x1*x2^2*x3^2*x4^2*xb^2*xg^2 + 32*x1*x2*x3^2*x4*xb*xg + 32*x1*x2^2*x3*x4*xb*xg - 4*x1^2*x2*x3^2*x4*xb^2*xg^2 - 4*x1^2*x2^2*x3*x4*xb^2*xg^2 + 4*x1^2*x2^2*x3^2*x4*xb*xg^2 + 4*x1^2*x2^2*x3^2*x4^2*xb*xg - 16*x1^2*x2*x3*x4*xb*xg^2:

P2 := (xb, xg, t1, t2, t3) -> 20*t1 + 3*x1^2*x4^2 + 40*t2*x1 - 20*t1*x1^2 + 20*t1*x4^2 - 3*x1^2 + 3*x4^2 + 40*t2*x1*x4^2 - 20*t1*x1^2*x4^2 - 3:

P3 := (xb, xg, t1, t2, t3) -> 100*t3 - 108*x2 - 108*x3 + 30*x4 + 71*x2^2*x3^2 + 71*x2^2*x4^2 - 71*x3^2*x4^2 + 100*t3*x2^2 + 100*t3*x3^2 + 100*t3*x4^2 + 108*x2*x3^2 + 108*x2^2*x3 - 108*x2*x4^2 - 30*x2^2*x4 - 108*x3*x4^2 - 30*x3^2*x4 + 71*x2^2 - 71*x3^2 - 71*x4^2 + 100*t3*x2^2*x3^2 + 100*t3*x2^2*x4^2 + 100*t3*x3^2*x4^2 + 108*x2*x3^2*x4^2 + 108*x2^2*x3*x4^2 + 30*x2^2*x3^2*x4 - 120*x2*x3*x4 + 71*x2^2*x3^2*x4^2 + 100*t3*x2^2*x3^2*x4^2 - 71:

P4 := (xb, xg, t1, t2, t3) -> 2130*x4 - 7668*x3 + 5000*t1^2*x3^2 + 5000*t1^2*x4^2 + 5000*t2^2*x3^2 + 5000*t2^2*x4^2 + 5000*t3^2*x3^2 + 5000*t3^2*x4^2 - 4091*x3^2*x4^2 - 7668*x3*x4^2 - 2130*x3^2*x4 + 5000*t1^2 + 5000*t2^2 + 5000*t3^2 - 4091*x3^2 - 4091*x4^2 + 5000*t1^2*x3^2*x4^2 + 5000*t2^2*x3^2*x4^2 + 5000*t3^2*x3^2*x4^2 - 4091:

# Substitute in EE pose
p1 := P1(xb, xg, t1, t2, t3):
p2 := P2(xb, xg, t1, t2, t3):
p3 := P3(xb, xg, t1, t2, t3):
p4 := P4(xb, xg, t1, t2, t3):

# Combine p1 and p2
p_12_1 := resultant(p1, p2, x1):
f_p_12_1 := factors(p_12_1):
r1 := select(p -> depends(p, x2) and depends(p, x3) and depends(p, x4), f_p_12_1[2])[1][1]:

# Combine in p3
r1_p3_2 := resultant(r1, p3, x2):
f_r1_p3_2 := factors(r1_p3_2):
r2 := select(p -> depends(p, x3) and depends(p, x4), f_r1_p3_2[2])[1][1]:

# Combine in p4
r2_p4_x3 := resultant(r2, p4, x3):
# f_r2_p4_x3 := factors(r2_p4_x3):
# r3 := f_r2_p4_x3[2][1][1]:

# Find zeros of resultant
fsolve(r2_p4_x3);
degree(r2_p4_x3);

f_o1 := factors(r2_p4_x3):
print("Degree of each factor:"):
map(v -> degree(v[1]), f_o1[2]);
print("Multiplicity of each factor:");
map(v -> v[2], f_o1[2]);
print("Roots of each factor:");
map(v -> [fsolve(v[1])], f_o1[2]);