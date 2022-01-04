% For your reference

% t is an array of timestamps for the simulation outputs
% y(:,1) is the x position of the tractor
% y(:,2) is the y position of the tractor
% y(:,3) is the angle of the tractor (angle along the body)
% y(:,4) is the x position of the trailer
% y(:,5) is the y position of the trailer
% y(:,6) is the angle the trailer hitch makes with the tractor
% 
% y(:,7) is the distance between the tractor's center (xt,yt) and the circle it is trying to track (L-offset)
% y(:,8) is the difference between the tractor's angle and the ideal angle for the circle it is trying to track (theta offset)
% y(:,9) is the difference between the trailer hitches' angle and the ideal angle for the circle it is trying to track (phi offset)
% 
% y(:,10) is the model reference for y(:,7) which is tracking r(t)
% y(:,11) is the model reference for y(:,8) which is going to 0
% 
% y(:,12) is the same as y(:,7) if things are going correctly (in the code, it is calculated a slightly different way, so agreement is a debugging tool)
% y(:,13) is -sigma*V*sin(y(:,11)) if things are going correctly (again, a debugging tool of sorts)
% 
% y(:,14) is L1hat, i.e. the estimate of L1 achieved through adaptation
% 
% y(:,15) is the model reference for y(:,9) which is going to ri(t)
% 
% y(:,16) is kyhat i.e. the estimate of the ideal gain for one of the terms in the trailer's controller (The physical meaning of ky is confusing and not-intuitive)
% y(:,17) is krhat i.e. the estimate of the ideal gain for another of the terms in the trailer's controller (The physical meaning of kr is confusing)
% y(:,18) is kFhat (same as above, this is a weird value which appears in the trailer's controller)
% -- All of these terms are estimates of something which when accurate results in cancellation of nonlinearities and decent control
% -- The estimates adapt when ri(t) is a persistent excitation
% 
% y(:,19) is the same as y(:,9) if we assume that the tractor is perfectly following its intended trajectory
% 
% ALSO I used some global variables to extract the controller values and the values which kyhat, krhat, and kFhat are approximating (They have formulas, those formulas are ugly)
% These are stored in the global variables u_record, t_record, and adaptive_params_target_record. Because the "t" values spit out by ode45() are different than those used by the simulation, 
% t and t_record are not the same!!!!

% Model parameters
V     = 1;
%sigma = 1;
L1    = 3.1;
L2    = 1;
a     = 0.5;

% Disturbance parameters
betaf = 0;
Vlr   = 0;
Vsr   = 0;
Vsi   = 0;

% Trajectory path following
R1 = @(t) 1000.*(t<20) + 3.*(25>t && t>=20)+ 1000.*(t>=25);
sigma = @(t) 1; % +1 for CCW, -1 for CW

r  = @(t) 0.1*sin(t) + 0.1*sin(2*t) + 0.1*sin(3*t);
ri = @(t) 0.1*sin(t) + 0.1*sin(2*t) + 0.1*sin(3*t);

zeta = 1;
phid =@(t) -sigma(t)*(pi - atan(R1(t)/a)-acos( (L2^2+a^2)/(2*R1(t)*sqrt(a^2+R1(t)^2)) ));

% Initial Conditions
xt0 = 0;
yt0 = 0;
thetat0 = 0;
phi0 = 0;
xi0 = xt0 - a*cos(thetat0) - L2*cos(thetat0+phi0);
yi0 = yt0 - a*sin(thetat0) - L2*sin(thetat0+phi0);

Los0     = 0;
thetaos0 = 0;
phios0   = phi0-phid(0);

L1hat0 = 2.5;

ym0 = phios0;
kyhat0 = 0;
krhat0 = 0;
kFhat0 = 0;
y0 = phios0;

ic = [xt0; yt0; thetat0; xi0; yi0; phi0; Los0; thetaos0; phios0; Los0; thetaos0; Los0; -sigma(0)*V*sin(thetaos0); L1hat0; ym0; kyhat0; krhat0; kFhat0; y0];

% control laws

omegan = 1;
zetad = 0.7;

% Tuning law
gamma = 15;

%Trailer control parameters
gammay = 25;
gammar = 25;
gammaF = 25;

%Trailer MRAC
am = 1;
bm = 1;


%% Simulation
global u_record t_record adaptive_params_target_record
u_record = [0,0];
t_record = [0];
adaptive_params_target_record = [0, 0, 0];

tspan = [0,500]
opts = odeset('MaxStep',1e-1);
[t,y] = ode45(@(t,x) tractor_trailer(t,x, V, L1, L2, a, betaf, Vlr, Vsr, Vsi, R1, sigma, zeta, phid, r, omegan, zetad, gamma, gammay, gammar, gammaF, ri, am, bm), tspan, ic, opts);

close all


pt  = [y(:,1), y(:,2)];
ptf = [y(:,1)+L1*cos(y(:,3)), y(:,2)+L1*sin(y(:,3))];
ptr = [y(:,1)-a*cos(y(:,3)), y(:,2)-a*sin(y(:,3))];
pti  = [y(:,4),y(:,5)];

px = [ptf(:,1), pt(:,1), ptr(:,1), pti(:,1)];
py = [ptf(:,2), pt(:,2), ptr(:,2), pti(:,2)];

pxt = [ptf(:,1), pt(:,1), ptr(:,1)];
pyt = [ptf(:,2), pt(:,2), ptr(:,2)];

pxi = [ptr(:,1), pti(:,1)];
pyi = [ptr(:,2), pti(:,2)];

Los     = y(:,7);
thetaos = y(:,8);
phios   = y(:,9);

%% Graphing
close all
subplot(7,2,[1,3,5]);
sn = 10; sn2 = 1000;sn3 = 1500;
figure(1)
plot(y(:,1),y(:,2),y(:,4),y(:,5))
% xlim([-10,10]);
% ylim([-10,10]);
hold on
plot(pxt(sn,:),pyt(sn,:),'k.-','LineWidth',2,'MarkerSize',20)
plot(pxi(sn,:),pyi(sn,:),'b.-','LineWidth',2,'MarkerSize',20)

plot(pxt(sn2,:),pyt(sn2,:),'k.-','LineWidth',2,'MarkerSize',20)
plot(pxi(sn2,:),pyi(sn2,:),'b.-','LineWidth',2,'MarkerSize',20)

plot(pxt(sn3,:),pyt(sn3,:),'k.-','LineWidth',2,'MarkerSize',20)
plot(pxi(sn3,:),pyi(sn3,:),'b.-','LineWidth',2,'MarkerSize',20)


%plot(R1(0)*cos(0:0.1:2*pi),R1(0)*sin(0:0.1:2*pi)-R1(0),"r--") 

legend("Tractor", "Trailer", "Path")

subplot(7,2,2);
plot(t, Los)
legend("Los")
subplot(7,2,4);
plot(t, thetaos, t, phios)
hold on 
plot(t, r(t), "--")
legend("Thetaos", "\phi_{os}")

subplot(7,2,6)
plot(t,y(:,14),t, ones(size(t))*L1);
legend("L1hat", "L1");

subplot(7,2,7:10)
plot(t_record,adaptive_params_target_record(:,1), "r--"); hold on; grid on;
plot(t_record,adaptive_params_target_record(:,2), "b--")
plot(t_record,adaptive_params_target_record(:,3), "g--")
plot(t, y(:,16),"r");
plot(t, y(:,17),"b");
plot(t, y(:,18),"g");


legend("ky^*","kr^*","kF^*", "kyhat", "krhat", "kFhat")

subplot(7,2,11)
plot(t_record, u_record(:,1));
hold on;
plot(t_record, u_record(:,2));
legend("Front wheel angle", "Trailer Wheel angle")


subplot(7,2,12)
plot(t, y(:,15)); 
hold on; 
plot(t, y(:,19));
legend("ym", "y");

% figure(2)
% subplot(3,1,1)
% plot(t, -sigma(t).*V.*sin(y(:,8)),t, y(:,13))
% subplot(3,1,2)
% plot(t, y(:,7),t, y(:,12))
% subplot(3,1,3)
% plot(t,y(:,14))

% Graphing interesting values


