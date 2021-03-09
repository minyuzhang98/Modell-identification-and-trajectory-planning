%%%%%%%%% objective function and penalty function caluculation
syms q1;
syms q2;
syms q3;
syms q4;

syms d_q1;
syms d_q2;
syms d_q3;
syms d_q4;

syms dd_q1;
syms dd_q2;
syms dd_q3;
syms dd_q4;
syms sigma tao t_begin;
% 
% syms D1;
% syms L2;
% syms L3;
% syms m2;
% syms m3;
% syms g;
% syms I1;
% syms I2;
% syms I3;

D1 = 0.086;
L2 = 0.146;
L3 = 0.187;
% m1 = 0.3;
m2 = 0.2;
m3 = 0.5;
g = 9.81;
I1 = 0.000277;
I2 = 0.000085;
I3 = 0.0000874;




syms a1;
syms a2;
syms a3;
syms a4;
syms a5;
syms a6;
syms a7;
syms a8;
syms a9;
syms a10;
syms a11;
syms a12;
syms a13;
syms a14;
syms a15;
syms a16;
syms a17;
syms a18;

syms t;
syms alpha;


syms X_obs Y_obs Z_obs

x1 = 0;
y1 = 0;
z1 = D1;

Xel = L2*cos(q1)*cos(q2);
Yel = L2*cos(q2)*sin(q1);
Zel = D1 + L2*sin(q2);

Xef = cos(q1)*(L3*cos(q2 + q3) + L2*cos(q2));
Yef = sin(q1)*(L3*cos(q2 + q3) + L2*cos(q2));
Zef = D1 + L3*sin(q2 + q3) + L2*sin(q2);


x2 = (x1 + Xel)/2;
y2 = (y1 + Yel)/2;
z2 = (z1 + Zel)/2;


x3 = (Xel + Xef)/2;
y3 = (Yel + Yef)/2;
z3 = (Zel + Zef)/2;




s_el = [Xel; Yel; Zel] - [x1; y1; z1];
e_s_el = simplify(s_el/sqrt(s_el(1)^2 + s_el(2)^2 + s_el(3)^2));

el_ef = [Xef; Yef; Zef] - [Xel; Yel; Zel];
e_el_ef = simplify(el_ef/sqrt(el_ef(1)^2 + el_ef(2)^2 + el_ef(3)^2));


s_obs = simplify([X_obs; Y_obs; Z_obs] - [x1; y1; z1]);
el_obs = simplify([X_obs; Y_obs; Z_obs] - [Xel; Yel; Zel]);

v_d1 = simplify(s_obs - sum(e_s_el.*s_obs) * e_s_el);
v_d2 = simplify(el_obs - sum(e_el_ef.* el_obs) * e_el_ef);


d1 = simplify(sqrt(sum(v_d1.*v_d1)));
d2 = simplify(sqrt(sum(v_d2.*v_d2)));


d_x2 = derivation(x2,q1,q2,q3,d_q1,d_q2,d_q3);
d_y2 = derivation(y2,q1,q2,q3,d_q1,d_q2,d_q3);
d_z2 = derivation(z2,q1,q2,q3,d_q1,d_q2,d_q3);

d_x3 = derivation(x3,q1,q2,q3,d_q1,d_q2,d_q3);
d_y3 = derivation(y3,q1,q2,q3,d_q1,d_q2,d_q3);
d_z3 = derivation(z3,q1,q2,q3,d_q1,d_q2,d_q3);
% % 
% % % 
v2 = simplify(sqrt(d_x2^2 + d_y2^2 + d_z2^2));
v3 = simplify(sqrt(d_x3^2 + d_y3^2 + d_z3^2));
% % % 
K = simplify(0.5*m2*v2^2+0.5*m3*v3^2+0.5*I1*d_q1^2+0.5*I2*d_q2^2+0.5*I3*d_q3^2);
% % % 
U = m2*g*(z2-D1)+m3*g*(z3-D1);
% % % 
L =  simplify(K - U);
% 
L_q1 = diff(L,q1);
L_q2 = diff(L,q2);
L_q3 = diff(L,q3);

% 
L_dq1 = diff(L,d_q1);
L_dq2 = diff(L,d_q2);
L_dq3 = diff(L,d_q3);


Q1 = a1*t^4 + a2*t^3 + a3*t^2 + a4*t +a5;
Q2 = a6*t^4 + a7*t^3 + a8*t^2 + a9*t +a10;
Q3 = a11*t^4 + a12*t^3 + a13*t^2 + a14*t +a15;



% % % 
d_Q1 = diff(Q1,t);
d_Q2 = diff(Q2,t);
d_Q3 = diff(Q3,t);

dd_Q1 = diff(d_Q1,t);
dd_Q2 = diff(d_Q2,t);
dd_Q3 = diff(d_Q3,t);

% 
L_q1 = subs(L_q1, [q1 q2 q3 d_q1 d_q2 d_q3], [Q1 Q2 Q3 d_Q1 d_Q2 d_Q3]);
L_q2 = subs(L_q2, [q1 q2 q3 d_q1 d_q2 d_q3], [Q1 Q2 Q3 d_Q1 d_Q2 d_Q3]);
L_q3 = subs(L_q3, [q1 q2 q3 d_q1 d_q2 d_q3], [Q1 Q2 Q3 d_Q1 d_Q2 d_Q3]);

% 
L_dq1 = subs(L_dq1, [q1 q2 q3 d_q1 d_q2 d_q3], [Q1 Q2 Q3 d_Q1 d_Q2 d_Q3]);
L_dq2 = subs(L_dq2, [q1 q2 q3 d_q1 d_q2 d_q3], [Q1 Q2 Q3 d_Q1 d_Q2 d_Q3]);
L_dq3 = subs(L_dq3, [q1 q2 q3 d_q1 d_q2 d_q3], [Q1 Q2 Q3 d_Q1 d_Q2 d_Q3]);


L_dq1_t = diff(L_dq1, t);
L_dq2_t = diff(L_dq2, t);
L_dq3_t = diff(L_dq3, t);


syms threshold p_tao



d1 = subs(d1, [q1 q2 q3 d_q1 d_q2 d_q3], [Q1 Q2 Q3 d_Q1 d_Q2 d_Q3]);
d2 = subs(d2, [q1 q2 q3 d_q1 d_q2 d_q3], [Q1 Q2 Q3 d_Q1 d_Q2 d_Q3]);
d1(p_tao, sigma, t_begin, X_obs, Y_obs, Z_obs) = subs(d1, [t], [sigma*p_tao+t_begin])
d2(p_tao, sigma, t_begin, X_obs, Y_obs, Z_obs) = subs(d2, [t], [sigma*p_tao+t_begin])


t1 = L_dq1_t - L_q1;
t2 = L_dq2_t - L_q2;
t3 = L_dq3_t - L_q3;


P1 = abs(t1 * d_Q1);

P2 = abs(t2 * d_Q2);

P3 = abs(t3 * d_Q3);

P = P1 + P2 + P3;


P(tao, sigma, t_begin) = subs(P, [t], [sigma*tao+t_begin])



function [deri] = derivation(y,q1,q2,q3,d_q1,d_q2,d_q3)

deri = simplify(diff(y,q1)*d_q1 + diff(y,q2)*d_q2 + diff(y,q3)*d_q3);

end
% 
