clear; clc; close all;
%% computing transformation matrix

syms tta1(t) tta2(t) tta3(t) tta4(t) tta5(t) tta6(t);

PI = sym(pi);
T1 = compute_DH_mod_tf(0,0,335,tta1(t));
T2 = compute_DH_mod_tf(75,-PI/2,0,tta2(t));
T3 = compute_DH_mod_tf(270,0,0,tta3(t)-PI);
T4 = compute_DH_mod_tf(90,-PI/2,295,tta4(t));
T5 = compute_DH_mod_tf(0,PI/2,0,tta5(t));
T6 = compute_DH_mod_tf(0,-PI/2,80,tta6(t));

T_ee = simplify( T1*T2*T3*T4*T5*T6);
R_ee = T_ee(1:3,1:3); 

% for i = 1:4
%     for j = 1:4
%         fprintf("T_ee(%d,%d)\n",i,j)
%         disp((T_ee(i,j)))
%     end
% end

%% evaluation
vpa(subs(T_ee,[tta1 tta2 tta3 tta4 tta5 tta6],[0 ,-PI/2,PI,0,0,0]),3);
%% Euler Angles

syms phi_E theta_E psi_E

% R_Euler = Rz(phi_E) * Ry(theta_E) * Rz(psi_E);

psi_E = atan2(R_ee(3,2),-R_ee(3,1));
phi_E = atan2(R_ee(2,3),R_ee(1,3));
syms e
e = R_ee(1,3)/cos(phi_E);
theta_E = atan2(e,R_ee(3,3));

phi_E = simplify(phi_E);
theta_E = simplify(theta_E);
psi_E = simplify(psi_E);

%% Roll Pitch Yaw Angles

syms phi_RPY theta_RPY psi_RPY

% R_RPY = Rz(phi_RPY) * Ry(theta_RPY) * Rx(psi_RPY);

psi_RPY = atan2(R_ee(3,2),R_ee(3,3));
phi_RPY = atan2(R_ee(2,1),R_ee(1,1));
syms rpy
rpy = R_ee(3,3)/cos(psi_RPY);
theta_RPY = atan2(-R_ee(3,1),rpy);

phi_RPY = simplify(phi_RPY);
theta_RPY = simplify(theta_RPY);
psi_RPY = simplify(psi_RPY);

%% Axis/Angle Representation

syms theta
k = sym('k',[3,1]);

thetak = acos((R_ee(1,1)+R_ee(2,2)+R_ee(3,3)-1)/2);
thetak = simplify(thetak);

k = (1/(2*sin(thetak))) * [R_ee(3,2)-R_ee(2,3); R_ee(1,3)-R_ee(3,1); R_ee(2,1)-R_ee(1,2)];
k = simplify(k);

%% Unit Quaternion Representation 

syms qr qx qy qz

qr = acos(((R_ee(1,1)+R_ee(2,2)+R_ee(3,3)+1)^(1/2))/2);
qr = simplify(qr);

qx = (R_ee(3,2)-R_ee(2,3))/(4*qr);
qy = (R_ee(1,3)-R_ee(3,1))/(4*qr);
qz = (R_ee(2,1)-R_ee(1,2))/(4*qr);

qx = simplify(qx);
qy = simplify(qy);
qz = simplify(qz);