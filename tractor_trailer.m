function [dxdt] = tractor_trailer(t,x, V, L1, L2, a, betaf, Vlr, Vsr, Vsi, R1, sigma, zeta, phid, r, omegan, zetad, gamma, gammay, gammar, gammaF, ri, am, bm)

dxdt = zeros(19,1);

xt     = x(1);
yt     = x(2);
thetat = x(3);
xi     = x(4);
yi     = x(5);
phi    = x(6);

Los     = x(7);
thetaos = x(8);
phios   = x(9);

% reference model for the tractor
xm1 = x(10);
xm2 = x(11);

x1  = x(12);
xi2 = x(13);

% adaptive tuning laws
L1hat = x(14);

% Trailer MRAC
ym = x(15);
kyhat = x(16);
krhat = x(17);
kFhat = x(18);
y = x(19);

% Tractor control law
f = (V^2-xi2^2)/(R1(t)+x1);
g = -sigma(t)*V*sqrt(V^2-xi2^2);
e1 = xm1-x1;
e2 = xm2-xi2;
mu = 1;
z = e2 + mu*e1;
beta = e1-omegan^2*xm1-2*zetad*omegan*xm2+omegan^2*r(t) + mu*e2;

if abs(g) > 0.1 && abs(z) > 0.1
    u = L1hat/g*(e1-omegan^2*xm1-2*zetad*omegan*xm2+omegan^2*r(t)-f+mu*e2 +z+mu^2/z);%(1/1+exp(-10*(abs(g)-1)))*(1/1+exp(-10*(abs(z)-1))) 
elseif abs(g) > 0.1 
    u = L1hat/g*(e1-omegan^2*xm1-2*zetad*omegan*xm2+omegan^2*r(t)-f+mu*e2);
else
    u = 0;
end

%Trailer control law
W1 = [-V/L2*cos(phid(t))+sigma(t)*V*a/(R1(t)*L2)*sin(phid(t)); -V/L2*sin(phid(t))-sigma(t)*V*a/(R1(t)*L2)*cos(phid(t)); -sigma(t)*V/R1(t)];
W2 = [V/L2*sin(phid(t))+sigma(t)*V*a/(R1(t)*L2)*cos(phid(t)); -V/L2*cos(phid(t))+sigma(t)*V*a/(R1(t)*L2)*sin(phid(t))];
Phi1 = [sin(thetaos); cos(thetaos); 1];
Phi2 = [sin(thetaos); cos(thetaos)];

ui = kyhat*y+krhat*ri(t)+kFhat;
%ui = 0;

 
%Control inputs 
delta = atan(u); %delta  = 0;
deltai = atan(ui);
global u_record t_record adaptive_params_target_record
u_record = [u_record; delta deltai];
t_record = [t_record; t];
adaptive_params_target_record = [adaptive_params_target_record; (-am)/(W2'*Phi2), (bm)/(W2'*Phi2),(-W1'*Phi1)/(W2'*Phi2)];

%Stationary reference frame ODEs
thetai = thetat + phi;
dxdt(1) = (V-Vlr)*cos(thetat)-Vsr*sin(thetat);
dxdt(2) = (V-Vlr)*sin(thetat)+Vsr*cos(thetat);
Thetat  = (V-Vlr)/L1*tan(delta + betaf) + Vsr/L1;
dxdt(3) = Thetat;

N = L2*cos(deltai);
M1 = -(V-Vlr)*sin(deltai+phi);
M2 = Vsr*cos(deltai+phi);
M3 = -Thetat*(a*cos(deltai+phi)+L2*cos(deltai));
M4 = -Vsi;

M3p = -a*Thetat*cos(deltai+phi);

Thetai = 1/N*(M1+M2+M3p+M4);

dxdt(4) = dxdt(1) + Thetat*a*sin(thetat) + Thetai*L2*sin(thetai);
dxdt(5) = dxdt(2) - Thetat*a*cos(thetat) - Thetai*L2*cos(thetai);
dxdt(6) = 1/N*(M1+M2+M3+M4);

% Rotating reference frame ODEs
dxdt(7) = -sigma(t)*abs(V-Vlr)*sin(thetaos)-sigma(t)*zeta*Vsr*cos(thetaos);
dxdt(8) = (V-Vlr)/L1*tan(delta+betaf)+Vsr/L1-sigma(t)*abs(V-Vlr)*cos(thetaos)/(R1(t)+Los)+sigma(t)*zeta*Vsr*sin(thetaos)/(R1(t)+Los);
dxdt(9) = 1/(L2*cos(deltai))*( -(V-Vlr)*sin(deltai+phios+phid(t))+Vsr*cos(deltai+phios+phid(t))-((V-Vlr)/L1*tan(delta+betaf)+Vsr/L1)*(a*cos(deltai+phios+phid(t))+L2*cos(deltai))-Vsi );


% MRAC Tractor 
dxdt(10) = xm2;
dxdt(11) = -omegan^2*xm1-2*zetad*omegan*xm2+omegan^2*r(t);
dxdt(12) = xi2;
dxdt(13) = f+1/L1*g*u; %%%%% This is L1 because what this is simulating is -sigma(t)*V*sin(thetaos), and thetaos is measureable. 


if L1hat > 0.1 || (L1hat == 0.1 && z*beta+z^2+z^2*mu^2 > 0)
    dxdt(14) = gamma*z*beta+z^2+z^2*mu^2; 
else
    dxdt(14) = 0;
end

% MRAC trailer
dxdt(15) = -am*ym + bm*ri(t);

e = ym - y;

dxdt(16) = gammay*sign(W2'*Phi2)*y*e;
dxdt(17) = gammar*sign(W2'*Phi2)*ri(t)*e;
dxdt(18) = gammaF*sign(W2'*Phi2)*e;

dxdt(19) = W1'*Phi1+W2'*Phi2*ui;

end

