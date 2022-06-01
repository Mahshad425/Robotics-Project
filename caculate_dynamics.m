clear; clc; close all;
%% computing transformation matrix

syms tta1(t) tta2(t) tta3(t) tta4(t) tta5(t) tta6(t);

tta1dot = diff(tta1(t),t,1);
tta2dot = diff(tta2(t),t,1);
tta3dot = diff(tta3(t),t,1);
tta4dot = diff(tta4(t),t,1);
tta5dot = diff(tta5(t),t,1);
tta6dot = diff(tta6(t),t,1);

tta1dotdot = diff(tta1dot,t,1);
tta2dotdot = diff(tta2dot,t,1);
tta3dotdot = diff(tta3dot,t,1);
tta4dotdot = diff(tta4dot,t,1);
tta5dotdot = diff(tta5dot,t,1);
tta6dotdot = diff(tta6dot,t,1);

q = [tta1(t); tta2(t); tta3(t); tta4(t); tta5(t); tta6(t)];
qdot = [tta1dot; tta2dot; tta3dot; tta4dot; tta5dot; tta6dot];
qdotdot = [tta1dotdot; tta2dotdot; tta3dotdot; tta4dotdot; tta5dotdot; tta6dotdot];

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

%w(i+1) = R'(i+1)*w(i) + diff(tta(i+1))*z(i+1)

w = sym('w',[3,1,6]);
w0 = [0; 0; 0];
w(:,:,1) = transpose(R(:,:,1))*w0 + [0; 0; qdot(1)];
for n=2:6
    w(:,:,n) = transpose(R(:,:,n))*w(:,:,n-1) + [0; 0; qdot(n)];
end

%v(i+1) = R'(i+1)*(v(i) + cross(w(i),p(i+1)))

v = sym('v',[3,1,6]);
v0 = [0; 0; 0];
v(:,:,1) = transpose(R(:,:,1)) * (v0 + cross(w0,p(:,:,1)));
for n=2:6
    v(:,:,n) = transpose(R(:,:,n)) * (v(:,:,n-1) + cross(w(:,:,n-1),p(:,:,n)));
end

%% Calculate Dynamic with Lagrange method

% Inertia Matrix
I = zeros(3,3,6);
I(:,:,1) = [335574671.26068  -25418821.76495  1193465.35991    ;-25418821.76495  429333178.02072   697646.34938   ;1193465.35991     697646.34938    324561265.61203];
I(:,:,2) = [561741324.02744 30852061.47573 462263673.51030;30852061.47573 1459457732.75253 27242765.48285;462263673.51030 27242765.48285 1006916175.70601];
I(:,:,3) = [137683420.48675  13243829.21878   -5923243.78632   ;13243829.21878   223903903.74202   -635534.81033  ;-5923243.78632    -635534.81033   256878102.25448];
I(:,:,4) = [219729519.56683  86670.13507      -1286682.42000   ;86670.13507      212690436.41858   -5883314.56115 ;1286682.42000     -5883314.56115  28969112.55938];
I(:,:,5) = [2383969.05153 -25.92426 -18.74015;-25.92426 1773748.83004 19.57785;-18.74015 19.57785 2769649.43788];
I(:,:,6) = [317521.64754 0.00000 0.00001;0.00000 318795.19255 -59.09778;0.00001 -59.09778 269236.26617];

% center of mass from solid
cm = zeros(3,1,6);
cm(:,:,1) = [37.33974;-7.58732;-43.55393];
cm(:,:,2) = [114.63352;6.19339;115.87354];
cm(:,:,3) = [45.20000;5.45123;-7.14639];
cm(:,:,4) = [0.64604;5.18455;-100.54730];
cm(:,:,5) = [-0.00067;6.24562;0.49862];
cm(:,:,6) = [0.00000;0.03014;-16.69195];  

% Mass parametrs
m = zeros(6,1);
m(1) = 44.07843482;
m(2) = 35954.94593;
m(3) = 33.50138982;
m(4) = 14.85437505;
m(5) = 2.65150879;
m(6) = 0.56015778;

%Jacobian of Masscenter
% T10 = vpa(subs(T10,[tta1 tta2 tta3 tta4 tta5 tta6],[0 ,-PI/2,PI,0,0,0]),3);
% T20 = vpa(subs(T20,[tta1 tta2 tta3 tta4 tta5 tta6],[0 ,-PI/2,PI,0,0,0]),3);
% T30 = vpa(subs(T30,[tta1 tta2 tta3 tta4 tta5 tta6],[0 ,-PI/2,PI,0,0,0]),3);
% T40 = vpa(subs(T40,[tta1 tta2 tta3 tta4 tta5 tta6],[0 ,-PI/2,PI,0,0,0]),3);
% T50 = vpa(subs(T50,[tta1 tta2 tta3 tta4 tta5 tta6],[0 ,-PI/2,PI,0,0,0]),3);
% T60 = vpa(subs(T60,[tta1 tta2 tta3 tta4 tta5 tta6],[0 ,-PI/2,PI,0,0,0]),3);

R0 = sym('R0',[3,3,6]);
p0 = sym('p0',[3,1,6]);
zc = sym('zc',[3,1,6]);
for n=1:6
    R0(:,:,n) = T0(1:3,1:3,n);
    p0(:,:,n) = T0(1:3,4,n);
    zc(:,:,n) = T0(1:3,3,n);    % T_cito0 = T_citoi * T_ito0
end

c0 = sym('c0',[3,1,6]);
for n=1:6
    c0(:,:,n) = R0(:,:,n)*cm(:,:,n) + p0(:,:,n);
end

Jvc = sym('Jvc',[3,6,6]);
Jwc = sym('Jwc',[3,6,6]);
for i=1:3
    for j=1:6
        for k=1:6
            Jvc(i,j,k) = 0;
            Jwc(i,j,k) = 0;
        end
    end
end

n = 1;
while n<7
    for i=1:n
        Jvc(:,i,n) = cross(zc(:,:,i),(c0(:,:,n)-p0(:,:,i)));
    end
    n = n+1;
end

n = 1;
while n<7
    for i=1:n
        Jwc(:,i,n) = zc(:,:,i);
    end
    n = n+1;
end

Jvc = simplify(Jvc,"Steps",150);
Jwc = simplify(Jwc,"Steps",150);

D = sym('D',[6,6]);
for i=1:6
    for j=1:6
        D(i,j) = 0;
    end
end
for i=1:6
    D = D + (m(i) * (transpose(Jvc(:,:,i))*Jvc(:,:,i))) + (transpose(Jwc(:,:,i))*I(:,:,i)*Jwc(:,:,i));
end

Cf = sym('Cf',[6,6,6]);
for k=1:6
    for i=1:6
        for j=1:6
            Cf(k,i,j) = 0.5 * ( diff(D(i,k),q(j)) + diff(D(j,k),q(i)) - diff(D(i,j),q(k)) );
        end
    end
end

C = sym('C',[6,6]);
for i=1:6
    for j=1:6
        C(i,j) = 0;
    end
end
for k=1:6
    for i=1:6
        for j=1:6
            C(k,i) = C(k,i) + Cf(k,i,j)*qdot(j);
        end
    end
end

gravity = [0; 0; 9.780327];
syms P
P = 0;
for i=1:6
    P = P + transpose(gravity)*c0(:,:,i)*m(i);
end
g = sym('g',[6,1]);
for i=1:6
    g(i) = diff(P,q(i));
end

% tau_L = sym('tau_L',[6,1]);
% tau_L = simplify(D*qdotdot + C*qdot + g, "steps", 50);
tau_L = D*qdotdot + C*qdot + g;
%% Calculate Dynamic with Newton-Euler method

% Inertia Matrix for Center of Mass
Ic = zeros(3,3,6);
Ic(:,:,1) = [249422846.19748  -12931021.40358  72877881.91444  ;-12931021.40358  284262207.25676  -13868409.42389  ;72877881.91444  -13868409.42389  260567153.68314];
Ic(:,:,2) = [77606666.13731 5325137.87555 -15325585.25146;5325137.87555 504223928.41652 1439708.87877;-15325585.25146 1439708.87877 533058706.52778];
Ic(:,:,3) = [134976948.88672  4989227.62644    4898269.14798   ;4989227.62644    153748275.34069  669566.98141     ;4898269.14798   669566.98141     187437893.85186];
Ic(:,:,4) = [69156078.82273   36916.60552      -321781.05587   ;36916.60552      62510075.61450   1860165.59197    ;-321781.05587   1860165.59197    28563633.22075];
Ic(:,:,5) = [2279880.33742 -14.87518 -17.85805;-14.87518 1773089.61115 -8237.69507;-17.85805 -8237.69507 2666219.94030];
Ic(:,:,6) = [161449.36513 0.00000 0.00000;0.00000 162723.41911 222.74706;0.00000 222.74706 269235.75720];

wdot = sym('wdot',[3,1,6]);
w0dot = [0; 0; 0];
wdot(:,:,1) = transpose(R(:,:,1))*w0dot + cross(transpose(R(:,:,1))*w0 , [0; 0; qdot(1)]) + [0; 0; qdotdot(1)];
for n=2:6
    wdot(:,:,n) = transpose(R(:,:,n))*wdot(:,:,n-1) + cross(transpose(R(:,:,n))*w(:,:,n-1) , [0; 0; qdot(n)]) + [0; 0; qdotdot(n)];
end

vdot = sym('vdot',[3,1,6]);
v0dot = [0; 0; 9.780327];
vdot(:,:,1) = transpose(R(:,:,1)) * ( v0dot + cross(w0,cross(w0,p(:,:,1))) + cross(w0dot,p(:,:,1)) );
for n=2:6
    vdot(:,:,n) = transpose(R(:,:,n)) * ( vdot(:,:,n-1) + cross(w(:,:,n-1),cross(w(:,:,n-1),p(:,:,n))) + cross(wdot(:,:,n-1),p(:,:,n)) );
end

vcdot = sym('vcdot',[3,1,6]);
for n=1:6
    vcdot(:,:,n) = cross(wdot(:,:,n),cm(:,:,n)) + cross(w(:,:,n),cross(w(:,:,n),cm(:,:,n))) + vdot(:,:,n);
end

F_capt = sym('F_capt',[3,1,6]);
for n=1:6
    F_capt(:,:,n) = m(n) * vcdot(:,:,n);
end

N_capt = sym('N_capt',[3,1,6]);
for n=1:6
    N_capt(:,:,n) = transpose(Ic(:,:,n))*wdot(:,:,n) + cross(w(:,:,n),transpose(Ic(:,:,n))*w(:,:,n));
end

f = sym('f',[3,1,6]);
f(:,:,6) = F_capt(:,:,6);
for n=5:-1:1
    f(:,:,n) = R(:,:,n+1)*f(:,:,n+1) + F_capt(:,:,n);
end

ni = sym('ni',[3,1,6]);
ni(:,:,6) = N_capt(:,:,6) + cross(cm(:,:,6),F_capt(:,:,6));
for n=5:-1:1
    ni(:,:,n) = N_capt(:,:,n) + R(:,:,n+1)*ni(:,:,n+1) + cross(cm(:,:,n),F_capt(:,:,n)) + cross(p(:,:,n+1),R(:,:,n+1)*f(:,:,n+1));
end

tau_NE = sym('tau_NE',[6,1]);
for n=1:6
    tau_NE(n) = ni(3,1,n);
end

%% Dynamics validation
e = tau_NE - tau_L;
% e1 = vpa(subs(e,[tta1(t) tta2(t) tta3(t) tta4(t) tta5(t) tta6(t)],[1,1,1,1,1,1]),3);
% e1 = [0; -1.91e-6; 0; 7.45e-9; 5.59e-9; 1.14e-13]     --->     Both of Euler-Lagrange & Newton-Euler methods gives same torques
% e2 = vpa(subs(e,[tta1(t) tta2(t) tta3(t) tta4(t) tta5(t) tta6(t)],[0,0,0,0,0,0]),3)
% e3 = vpa(subs(e,[tta1(t) tta2(t) tta3(t) tta4(t) tta5(t) tta6(t)],[0,-PI/2,PI,0,0,0]),3)
% e2,e3 = [0; 0; 0; 0; 0; 0]
% e4 = vpa(subs(e,[tta1(t) tta2(t) tta3(t) tta4(t) tta5(t) tta6(t)],[3*PI/4,-PI/2,PI,PI/3,2*PI/3,-PI]),3)
% e4 = [-1.22e-7; 2.38e-7; 2.38e-7; 0; 0; 0]






