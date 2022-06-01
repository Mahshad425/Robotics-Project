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

%% Computing Jacobian with Alternative method

J = sym('J',[6,6]);

for n=1:6
    J(1:3,n) = cross(T0(1:3,3,n),(T0(1:3,4,6)-T0(1:3,4,n)));
end
for n=1:6
    J(4:6,n) = T0(1:3,3,n);
end

J = simplify(J,"Steps",150);
%% Singularity of Robot

Jdet = simplify(det(J),"Steps",50);

% Acording to datasheet:
% theta2 [-190 , 45]Degree
% theta3 [-29 , 256]Degree
% theta5 [-120 , 120]Degree

% det(J) = 3375*sin(tta5)*(3805*sin(tta2) - 3805*cos(atan(3157/2124) + tta2 + 2*tta3) + 54*3805^(1/2)*cos(atan(18/59) + tta2 - tta3) + 54*3805^(1/2)*cos(tta2 - atan(18/59) + tta3) + 30*3805^(1/2)*cos(atan(18/59) - tta3))
% as we can see. In the determinant, only tta2, tta3 & tta5 appear.
% Suppose that : x=atan(18/59)    &    y=atan(3157/2124)
% det(J) = 3375*3805^(1/2)  sin(tta5)  ( 27  cos(x-tta3)  (cos(tta2) + (10/9)) + 3805^(1/2) * (sin(tta2) + cos(tta2)*cos(y + 2*tta3) - sin(tta2)*sin(y + 2*tta3))
% As a result:
% sin(tta5) = 0 ---> tta5 = 0

syms theta2 theta3

singularity = 27 * cos(atan(18/59)-theta3) * (cos(theta2) + (10/9)) + 3805^(1/2) * (sin(theta2) + cos(theta2 + 2*theta3 + atan(3157/2124)));

solve(singularity==0,theta2,"real",1)
solve1=[-2*atan(((108*3805^(1/2)*cos(theta3 + 22962614765628757/18014398509481984) - 342*cos(2*theta3 - 5334336242485005/9007199254740992) - 30440*sin(2*theta3 + 2203534815392969/2251799813685248) + 108*3805^(1/2)*cos(3*theta3 + 12293942280658747/18014398509481984) + 30098)^(1/2) + 2*3805^(1/2) - 2*3805^(1/2)*sin(2*theta3 + 2203534815392969/2251799813685248))/(6*cos(theta3 - 5334336242485005/18014398509481984) - 2*3805^(1/2)*cos(2*theta3 + 2203534815392969/2251799813685248))) 2*atan(((108*3805^(1/2)*cos(theta3 + 22962614765628757/18014398509481984) - 342*cos(2*theta3 - 5334336242485005/9007199254740992) - 30440*sin(2*theta3 + 2203534815392969/2251799813685248) + 108*3805^(1/2)*cos(3*theta3 + 12293942280658747/18014398509481984) + 30098)^(1/2) - 2*3805^(1/2) + 2*3805^(1/2)*sin(2*theta3 + 2203534815392969/2251799813685248))/(6*cos(theta3 - 5334336242485005/18014398509481984) - 2*3805^(1/2)*cos(2*theta3 + 2203534815392969/2251799813685248)))]

eval = subs(singularity,[theta2],[solve1(1)])
eval2 = subs(eval,[theta3],[0]) %0 or other number
double(eval2) %eval2 ~= 0


%solve(singularity==0,theta2,'ReturnConditions',1)
 %struct with fields:
    %theta2: [2×1 sym]
    %parameters: k
    %conditions: [2×1 sym]
%ans.conditions
%in(k, 'integer') & 27*exp(theta3*1i)*exp(5334336242485005i/18014398509481984) + 27*exp(theta3*3i)*exp(-5334336242485005i/18014398509481984) + 2*3805^(1/2)*exp(theta3*4i)*exp(2203534815392969i/2251799813685248) ~= 3805^(1/2)*exp(theta3*2i)*2i & (- 30098*exp(theta3*3i) + exp(theta3*1i)*exp(-2203534815392969i/2251799813685248)*15220i - exp(theta3*5i)*exp(2203534815392969i/2251799813685248)*15220i - 54*3805^(1/2)*exp(-12293942280658747i/18014398509481984) + 171*exp(theta3*1i)*exp(5334336242485005i/9007199254740992) + 171*exp(theta3*5i)*exp(-5334336242485005i/9007199254740992) - 54*3805^(1/2)*exp(theta3*2i)*exp(-22962614765628757i/18014398509481984) - 54*3805^(1/2)*exp(theta3*4i)*exp(22962614765628757i/18014398509481984) - 54*3805^(1/2)*exp(theta3*6i)*exp(12293942280658747i/18014398509481984))^(1/2) ~= 30*exp((theta3*1i)/2)*exp(5334336242485005i/18014398509481984) + 30*exp((theta3*5i)/2)*exp(-5334336242485005i/18014398509481984) & (- 30098*exp(theta3*3i) + exp(theta3*1i)*exp(-2203534815392969i/2251799813685248)*15220i - exp(theta3*5i)*exp(2203534815392969i/2251799813685248)*15220i - 54*3805^(1/2)*exp(-12293942280658747i/18014398509481984) + 171*exp(theta3*1i)*exp(5334336242485005i/9007199254740992) + 171*exp(theta3*5i)*exp(-5334336242485005i/9007199254740992) - 54*3805^(1/2)*exp(theta3*2i)*exp(-22962614765628757i/18014398509481984) - 54*3805^(1/2)*exp(theta3*4i)*exp(22962614765628757i/18014398509481984) - 54*3805^(1/2)*exp(theta3*6i)*exp(12293942280658747i/18014398509481984))^(1/2) + 30*exp((theta3*1i)/2)*exp(5334336242485005i/18014398509481984) + 30*exp((theta3*5i)/2)*exp(-5334336242485005i/18014398509481984) ~= 0
%in(k, 'integer') & 27*exp(theta3*1i)*exp(5334336242485005i/18014398509481984) + 27*exp(theta3*3i)*exp(-5334336242485005i/18014398509481984) + 2*3805^(1/2)*exp(theta3*4i)*exp(2203534815392969i/2251799813685248) ~= 3805^(1/2)*exp(theta3*2i)*2i & (- 30098*exp(theta3*3i) + exp(theta3*1i)*exp(-2203534815392969i/2251799813685248)*15220i - exp(theta3*5i)*exp(2203534815392969i/2251799813685248)*15220i - 54*3805^(1/2)*exp(-12293942280658747i/18014398509481984) + 171*exp(theta3*1i)*exp(5334336242485005i/9007199254740992) + 171*exp(theta3*5i)*exp(-5334336242485005i/9007199254740992) - 54*3805^(1/2)*exp(theta3*2i)*exp(-22962614765628757i/18014398509481984) - 54*3805^(1/2)*exp(theta3*4i)*exp(22962614765628757i/18014398509481984) - 54*3805^(1/2)*exp(theta3*6i)*exp(12293942280658747i/18014398509481984))^(1/2) ~= 30*exp((theta3*1i)/2)*exp(5334336242485005i/18014398509481984) + 30*exp((theta3*5i)/2)*exp(-5334336242485005i/18014398509481984) & (- 30098*exp(theta3*3i) + exp(theta3*1i)*exp(-2203534815392969i/2251799813685248)*15220i - exp(theta3*5i)*exp(2203534815392969i/2251799813685248)*15220i - 54*3805^(1/2)*exp(-12293942280658747i/18014398509481984) + 171*exp(theta3*1i)*exp(5334336242485005i/9007199254740992) + 171*exp(theta3*5i)*exp(-5334336242485005i/9007199254740992) - 54*3805^(1/2)*exp(theta3*2i)*exp(-22962614765628757i/18014398509481984) - 54*3805^(1/2)*exp(theta3*4i)*exp(22962614765628757i/18014398509481984) - 54*3805^(1/2)*exp(theta3*6i)*exp(12293942280658747i/18014398509481984))^(1/2) + 30*exp((theta3*1i)/2)*exp(5334336242485005i/18014398509481984) + 30*exp((theta3*5i)/2)*exp(-5334336242485005i/18014398509481984) ~= 0
%ans.theta2
%2*pi*k - log(-((- 30098*exp(theta3*3i) + exp(theta3*1i)*exp(-2203534815392969i/2251799813685248)*15220i - exp(theta3*5i)*exp(2203534815392969i/2251799813685248)*15220i - 54*3805^(1/2)*exp(-12293942280658747i/18014398509481984) + 171*exp(theta3*1i)*exp(5334336242485005i/9007199254740992) + 171*exp(theta3*5i)*exp(-5334336242485005i/9007199254740992) - 54*3805^(1/2)*exp(theta3*2i)*exp(-22962614765628757i/18014398509481984) - 54*3805^(1/2)*exp(theta3*4i)*exp(22962614765628757i/18014398509481984) - 54*3805^(1/2)*exp(theta3*6i)*exp(12293942280658747i/18014398509481984))^(1/2) + 30*exp((theta3*1i)/2)*exp(5334336242485005i/18014398509481984) + 30*exp((theta3*5i)/2)*exp(-5334336242485005i/18014398509481984))/(27*exp((theta3*1i)/2)*exp(5334336242485005i/18014398509481984) + 27*exp((theta3*5i)/2)*exp(-5334336242485005i/18014398509481984) - 3805^(1/2)*exp((theta3*3i)/2)*2i + 2*3805^(1/2)*exp((theta3*7i)/2)*exp(2203534815392969i/2251799813685248)))*1i
%2*pi*k - log(-(- (- 30098*exp(theta3*3i) + exp(theta3*1i)*exp(-2203534815392969i/2251799813685248)*15220i - exp(theta3*5i)*exp(2203534815392969i/2251799813685248)*15220i - 54*3805^(1/2)*exp(-12293942280658747i/18014398509481984) + 171*exp(theta3*1i)*exp(5334336242485005i/9007199254740992) + 171*exp(theta3*5i)*exp(-5334336242485005i/9007199254740992) - 54*3805^(1/2)*exp(theta3*2i)*exp(-22962614765628757i/18014398509481984) - 54*3805^(1/2)*exp(theta3*4i)*exp(22962614765628757i/18014398509481984) - 54*3805^(1/2)*exp(theta3*6i)*exp(12293942280658747i/18014398509481984))^(1/2) + 30*exp((theta3*1i)/2)*exp(5334336242485005i/18014398509481984) + 30*exp((theta3*5i)/2)*exp(-5334336242485005i/18014398509481984))/(27*exp((theta3*1i)/2)*exp(5334336242485005i/18014398509481984) + 27*exp((theta3*5i)/2)*exp(-5334336242485005i/18014398509481984) - 3805^(1/2)*exp((theta3*3i)/2)*2i + 2*3805^(1/2)*exp((theta3*7i)/2)*exp(2203534815392969i/2251799813685248)))*1i


% theta2 = 0.047659495549047566029893493506547 & theta3 = 0.79051484366232575942234064481128
% subs(solve1,[theta3],[0.79051484366232575942234064481128])
% double(ans)
%-1.618980055346616
%0.047659495549047
