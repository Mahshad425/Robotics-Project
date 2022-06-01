clear; clc; close all;
%% computing transformation matrix

syms tta1(t) tta2(t) tta3(t) tta4(t) tta5(t) tta6(t);

tta1dot = diff(tta1(t),t,1);
tta2dot = diff(tta2(t),t,1);
tta3dot = diff(tta3(t),t,1);
tta4dot = diff(tta4(t),t,1);
tta5dot = diff(tta5(t),t,1);
tta6dot = diff(tta6(t),t,1);

q = [tta1(t); tta2(t); tta3(t); tta4(t); tta5(t); tta6(t)];
qdot = [tta1dot; tta2dot; tta3dot; tta4dot; tta5dot; tta6dot];

PI = sym(pi);
T1 = compute_DH_mod_tf(0,0,335,tta1(t));
T2 = compute_DH_mod_tf(75,-PI/2,0,tta2(t));
T3 = compute_DH_mod_tf(270,0,0,tta3(t)-PI);
T4 = compute_DH_mod_tf(90,-PI/2,295,tta4(t));
T5 = compute_DH_mod_tf(0,PI/2,0,tta5(t));
T6 = compute_DH_mod_tf(0,-PI/2,80,tta6(t));

T = cat(3,T1,T2,T3,T4,T5,T6);

T_ee = simplify( T1*T2*T3*T4*T5*T6);
% for i = 1:4
%     for j = 1:4
%         fprintf("T_ee(%d,%d)\n",i,j)
%         disp((T_ee(i,j)))
%     end
% end

%% evaluation
vpa(subs(T_ee,[tta1 tta2 tta3 tta4 tta5 tta6],[0 ,-PI/2,PI,0,0,0]),3);
%% Some of Equations

T_eedot = diff(T_ee,t);
P_eedot = T_eedot(1:3,4);

R_ee = T_ee(1:3,1:3);    
R_eedot = T_eedot(1:3,1:3);
S = R_eedot*transpose(R_ee);

T0 = sym('T0',[4,4,6]);
T0(:,:,1) = T(:,:,1);
for n=2:6
    T0(:,:,n) = T0(:,:,n-1) * T(:,:,n);
end

R = sym('R',[3,3,6]);
p = sym('p',[3,1,6]);
for n=1:6
    R(:,:,n) = T(1:3,1:3,n);
    p(:,:,n) = T(1:3,4,n);
end

%% Computing Jacobian with Derivation method

Jv1 = sym('Jv1',[3,6]);
for n=1:3
    for i=1:6
        Jv1(n,i) = diff(P_eedot(n),qdot(i));
    end
end
simplify(Jv1);  

w1 = sym('w1',[3,1]); 
w1(1) = -S(2,3);
w1(2) = S(1,3);
w1(3) = -S(1,2);

Jw1 = sym('Jw1',[3,6]);  
for n=1:3
    for i=1:6
        Jw1(n,i) = diff(w1(n),qdot(i));
    end
end
simplify(Jw1);

J1 = simplify(cat(1,Jv1,Jw1),"Steps",150);

%% Computing Jacobian with Alternative method

J2 = sym('J2',[6,6]);

for n=1:6
    J2(1:3,n) = cross(T0(1:3,3,n),(T0(1:3,4,6)-T0(1:3,4,n)));
end
for n=1:6
    J2(4:6,n) = T0(1:3,3,n);
end

J2 = simplify(J2,"Steps",150);

%% Computing Jacobian with Angular & Liner Velocity

% w(i+1) = R'(i+1)*w(i) + diff(tta(i+1))*z(i+1)

w = sym('w',[3,1,6]);
w0 = [0; 0; 0];
w(:,:,1) = transpose(R(:,:,1))*w0 + [0; 0; qdot(1)];
for n=2:6
    w(:,:,n) = transpose(R(:,:,n))*w(:,:,n-1) + [0; 0; qdot(n)];
end

% v(i+1) = R'(i+1)*(v(i) + cross(w(i),p(i+1)))

v = sym('v',[3,1,6]);
v0 = [0; 0; 0];
v(:,:,1) = transpose(R(:,:,1)) * (v0 + cross(w0,p(:,:,1)));
for n=2:6
    v(:,:,n) = transpose(R(:,:,n)) * (v(:,:,n-1) + cross(w(:,:,n-1),p(:,:,n)));
end

% Jacobian

Jw3 = sym('Jw3',[3,6]);
for n=1:3
    for i=1:6
        Jw3(n,i) = diff(w(n,1,6),qdot(i));
    end
end
simplify(Jw3);

Jv3 = sym('Jv3',[3,6]);
for n=1:3
    for i=1:6
        Jv3(n,i) = diff(v(n,1,6),qdot(i));
    end
end
simplify(Jv3);

% J3 = simplify(cat(1,Jv3,Jw3),"Steps",50);

Jv0 = R_ee * Jv3;
Jw0 = R_ee * Jw3;
J0 = simplify(cat(1,Jv0,Jw0),"Steps",200);
%% Jacobian validation
val12 = J1 - J2;
% val12 = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]
val23 = J2 - J3;
% val23 = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]