syms y1 y2 y3 y4 y5 y6 y7 y8 y9 t 
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


%%without noise
cof1 = [ 0.00415787  0.0040592   0.00220592 -0.00458814 -0.00250607  0.00071375  0.01260407  0.00314527 -0.00743689  0.00141626 -0.00416885]

cof2 = [ 0.00595673  0.00507937  0.01113552 -0.00807071  0.00650885 -0.00389041 -0.01952883 -0.00408277  0.83789304  0.4794819  -0.00685992]

cof3 = [ 0.00185267  0.00497993  0.00436784  0.00197585 -0.00177349  0.00328277  0.00583527  0.45633212  0.004329  ]


feature1 = y4 * [y7*cos(y3) y7*cos(2*y2+y3) y7*cos(2*y2+2*y3) y8*cos(y2+y3) y5*y6*cos(y3) y4^2*sin(2*y2+2*y3) y5^2*sin(y3) y6^2*sin(y3) y7*sin(y3) y4*y5*sin(y3) y5*y6*sin(2*y2+y3)]


feature2 = y5 * [y8*cos(y3) y8*cos(2*y2+2*y3) y9*cos(y3) y5*y6*cos(2*y2) y7*sin(2*y2+y3) y7*sin(2*y2+2*y3) y8*sin(y3) y5*y6*sin(y3) cos(y2) cos(y2+y3) y7]


feature3 = y6 * [y7*cos(2*y2+y3) y8*cos(y3) y8*cos(y2+y3) y5*y6*cos(y3) y5*y6*cos(2*y3) y4^2*sin(2*y2+2*y3) y5^2*sin(y3) cos(y2+y3) y9]


%%%with 5% measurement noise
% cof1 = [ 0.00982341  0.00112673  0.0105573  -0.00101766 -0.00880836  0.00299311 -0.00664093 -0.0022332  -0.00513588 -0.00404829 -0.00473655]
% 
% feature1 = y4 * [y7*cos(y3) y7*cos(2*y2) y7*cos(2*y2+y3) y8*cos(y2+y3) y4*y5*cos(y2+y3) y5^2*sin(y3) y7*sin(y3) y7*sin(2*y2+2*y3) y5*y6*sin(y3) y5*y6*sin(2*y2+y3) y5*y6*sin(2*y2+2*y3)]
% 
% cof2 = [ 0.00703253  0.00188951  0.00816304 -0.00901684  0.00991944 -0.00588536 -0.02400965 -0.00581517  0.83493659  0.48395803 -0.00859779]
% 
% feature2 = y5 * [y8*cos(y3) y9*cos(y2) y9*cos(y3) y5*y6*cos(2*y2) y7*sin(2*y2+y3) y7*sin(2*y2+2*y3) y8*sin(y3) y5*y6*sin(y3) cos(y2) cos(y2+y3) y7]
% 
% cof3 = [0.0089653  0.00148149 0.00308855 0.00281317 0.00543578 0.00303632 0.45307103 0.00261267 0.00320265]
% 
% feature3 = y6 * [y8*cos(y3) y9*cos(y2+y3) y4^2*sin(2*y2+y3) y4^2*sin(2*y2+2*y3) y4*y5*sin(y3) cos(y2) cos(y2+y3) y8 y9]


%%with 20% measurement noise
% cof1 = [ 0.00310892  0.00767966  0.00896143 -0.00125902 -0.01679052 -0.0118086   0.00113808 -0.01398205 -0.00550262 -0.00465585 -0.00581769  0.00345448]
% 
% feature1 = y4 * [y7*cos(y2) y7*cos(y3) y7*cos(2*y2+y3) y9*cos(y3) y4*y5*cos(y2) y4*y6*cos(2*y2+y3) y9*sin(2*y2+y3) y4*y5*sin(2*y2+y3) y5*y6*sin(y3) y5*y6*sin(2*y2+y3) y5*y6*sin(2*y2+2*y3) y7]
% 
% cof2 = [ 0.00653756  0.0084816  -0.00566202  0.00610501  0.00247268  0.00369222  0.84921949  0.47703245 -0.00281969  0.01783928  0.00778038]
% 
% feature2 = y5 * [y7*cos(2*y2+2*y3) y9*cos(2*y2) y5*y6*cos(2*y2) y6^2*sin(2*y2+y3) y7*sin(y3) y9*sin(2*y2+2*y3) cos(y2) cos(y2+y3) y4^2*sin(2*y2) y8 y9]
% 
% cof3 = [ 0.00848382  0.00130296  0.00239829  0.00295957 -0.00298258  0.00531639  0.00275109  0.45369728  0.00338323]
% 
% feature3 = y6 * [y8*cos(y3) y9*cos(y2+y3) y4^2*sin(2*y2+y3) y4^2*sin(2*y2+2*y3) y8*sin(y3) y4*y5*sin(y3) cos(y2) cos(y2+y3) y9]


q1 = a1*t^4 + a2*t^3 + a3*t^2 + a4*t +a5;
q2 = a6*t^4 + a7*t^3 + a8*t^2 + a9*t +a10;
q3 = a11*t^4 + a12*t^3 + a13*t^2 + a14*t +a15;

dq1 = diff(q1, t);
dq2 = diff(q2, t);
dq3 = diff(q3, t);

ddq1 = diff(dq1, t);
ddq2 = diff(dq2, t);
ddq3 = diff(dq3, t);

feature1 = subs(feature1 , [y1 y2 y3 y4 y5 y6 y7 y8 y9], [q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3]);
feature2 = subs(feature2 , [y1 y2 y3 y4 y5 y6 y7 y8 y9], [q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3]);
feature3 = subs(feature3 , [y1 y2 y3 y4 y5 y6 y7 y8 y9], [q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3]);

p1 = sum(cof1 .* feature1);
p2 = sum(cof2 .* feature2);
p3 = sum(cof3 .* feature3);

syms sigma tao t_begin;

P = abs(p1) + abs(p2) + abs(p3);
P(tao, sigma, t_begin) = subs(P, [t], [sigma*tao+t_begin])

% E_total = ( E(0.2)+E(0.4)+E(0.6)+E(0.8)+E(1.0)+E(1.2)+E(1.4)+E(1.6)+E(1.8)+E(2.0) ) * 0.2;




