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
% for i = 1:4
%     for j = 1:4
%         fprintf("T_ee(%d,%d)\n",i,j)
%         disp((T_ee(i,j)))
%     end
% end

%% evaluation
vpa(subs(T_ee,[tta1 tta2 tta3 tta4 tta5 tta6],[0 ,-PI/2,PI,0,0,0]),3);