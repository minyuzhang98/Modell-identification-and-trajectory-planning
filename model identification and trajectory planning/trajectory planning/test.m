%%%% from p1 t0 p2
theta_0 = [pi/2, pi/5, -pi/4];
d_theta_0 = [0, 0, 0];
theta_end = [-pi/2, 2*pi/5, -pi/3];


%%% obstacle
P = zeros(201,3);
P(:,1) = linspace(0.15,0.18,201);
P(:,2) = linspace(0.4,-0.35,201);
P(:,3) = linspace(0.22,0.28,201);
object = P(1,:);


dtg = DTG(theta_0, d_theta_0, theta_end);
dtg.Object = P;
dtg.Point_number = 200;
init = dtg.getResult();
dtg.dynamicTrajGen();
dtg.plot_distance();
