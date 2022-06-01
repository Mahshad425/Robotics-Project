%% computing transformation matrix
clear;clc;close all;
syms tta1 tta2 tta3 tta4 tta5 tta6
PI = sym(pi);
T1 = compute_DH_mod_tf(0,0,0.335,tta1);
T2 = compute_DH_mod_tf(0.075,-PI/2,0,tta2);
T3 = compute_DH_mod_tf(0.27,0,0,tta3-PI);
T4 = compute_DH_mod_tf(0.090,-PI/2,0.295,tta4);
T5 = compute_DH_mod_tf(0,PI/2,0,tta5);
T6 = compute_DH_mod_tf(0,-PI/2,0.080,tta6);

T10 = T1;
T20 = T1*T2;
T30 = T20 * T3;
T40 = T30 * T4;
T50 = T40 * T5;
T60 = T50 * T6;
T_ee = (simplify( T1*T2*T3*T4*T5*T6))
%% animation using just code
param_Js_position = [T10(1:3,4),T20(1:3,4),T30(1:3,4),T40(1:3,4),T50(1:3,4),T60(1:3,4)];
NoP = 1000;
XJ = zeros(1000,6);
YJ = zeros(1000,6);
ZJ = zeros(1000,6);
i = 1;
for i1 = linspace(-170/180*pi,170/180*pi,NoP)
    second_joint_iter = pi/3*sin(4*(2*pi/(-190/180*pi-45/180*pi))*i1);
    current_config = [i1 ,-pi/2+second_joint_iter, pi, 0, 0 ,0 ];
    current_joints_position = subs(param_Js_position,[tta1 tta2 tta3 tta4 tta5 tta6],current_config);
    XJ(i,:) = current_joints_position(1,:);
    YJ(i,:) = current_joints_position(2,:);
    ZJ(i,:) = current_joints_position(3,:);
    i=i+1
end 

figure(20)
framesPerSecond = 50;
r = rateControl(framesPerSecond);
for index = 1:NoP
    plot3(XJ(1:index,6),YJ(1:index,6),ZJ(1:index,6),'.r')
    grid on
    hold on
    plot3([0,XJ(index,1:6)],[0,YJ(index,1:6)],[0,ZJ(index,1:6)],'b',"LineWidth",1)
    view([-125.02 17.07])
    axis([-0.6 0.6 -0.6 0.6 -0.25 1])
    waitfor(r)
    hold off
end

%% building robot using toolbox
%tta3_prime = tta3-PI;
% home_pose = [ 90, -155 , 245-180 , -90, -90, 0] / 180 * pi
home_pose = [ 0, -90, 180-180 , 0,0, 0] / 180 * pi

mdh_params = [0,    0,      335/1000,    0;
              75/1000,   -pi/2,  0,      0;
              270/1000,  0,      0,      0;
              90/1000,   -pi/2,  295/1000,    0;
              0,    pi/2,   0,      0;
              0,    -pi/2,  80/1000,     0];
lnk1 = rigidBody('Link1');
jnt1 = rigidBodyJoint('joint1','revolute');
jnt1.HomePosition = home_pose(1);
jnt1.PositionLimits = [-170 170] / 180 * pi;
setFixedTransform(jnt1,mdh_params(1,:),'mdh');
lnk1.Joint = jnt1;

lnk2 = rigidBody('Link2');
jnt2 = rigidBodyJoint('joint2','revolute');
jnt2.HomePosition = home_pose(2);
jnt2.PositionLimits = [-190 45] / 180 * pi;
setFixedTransform(jnt2,mdh_params(2,:),'mdh');
lnk2.Joint = jnt2;

lnk3 = rigidBody('Link3');
jnt3 = rigidBodyJoint('joint3','revolute');
jnt3.PositionLimits = [-29-180 256-180] / 180 * pi;
jnt3.HomePosition = home_pose(3);
setFixedTransform(jnt3,mdh_params(3,:),'mdh');
lnk3.Joint = jnt3;

lnk4 = rigidBody('Link4');
jnt4 = rigidBodyJoint('joint4','revolute');
jnt4.HomePosition = home_pose(4);
jnt4.PositionLimits = [-190 190] / 180 * pi;
setFixedTransform(jnt4,mdh_params(4,:),'mdh');
lnk4.Joint = jnt4;

lnk5 = rigidBody('Link5');
jnt5 = rigidBodyJoint('joint5','revolute');
jnt5.HomePosition = home_pose(5);
jnt5.PositionLimits = [-120 120]/ 180 * pi;
setFixedTransform(jnt5,mdh_params(5,:),'mdh');
lnk5.Joint = jnt5;

lnk6 = rigidBody('Link6');
jnt6 = rigidBodyJoint('joint6','revolute');
jnt6.HomePosition = home_pose(6);
jnt6.PositionLimits = [-360 360] / 180 * pi;
setFixedTransform(jnt6,mdh_params(6,:),'mdh');
lnk6.Joint = jnt6;

robot = rigidBodyTree('DataFormat','row');
addBody(robot,lnk1,'base')
addBody(robot,lnk2,'Link1')
addBody(robot,lnk3,'Link2')
addBody(robot,lnk4,'Link3')
addBody(robot,lnk5,'Link4')
addBody(robot,lnk6,'Link5')
%% showing robot at transportation configuration
figure(1)
show(robot,'Frames',"on");
axis([-0.5 0.5 -0.5 0.5 -0.25 1])
view([-125.02 17.07])
%% an animation of robot rotating and painting
figure(2)
axis manual
% axis([-0.6 0.6 -0.6 0.6 -0.25 1])
% view([-125.02 17.07])
framesPerSecond = 20;
r = rateControl(framesPerSecond);
% show(robot,[jnt1.PositionLimits(1) -pi/2 pi-pi 0 0 0 ]);
for i1 = linspace(jnt1.PositionLimits(1),jnt1.PositionLimits(2),1000)
    second_joint_iter = pi/3*sin(4*(2*pi/(jnt1.PositionLimits(2)-jnt1.PositionLimits(1)))*i1);
    current_config = [i1 ,-pi/2+second_joint_iter, pi-pi, 0, 0 ,0 ];
    
    ax= show(robot,current_config,'FastUpdate',1,'PreservePlot',0);
    current_Tee = getTransform(robot,current_config,'Link6','base');
    current_x = current_Tee(1,4);
    current_y = current_Tee(2,4);
    current_z = current_Tee(3,4);
    hold on
    plot3(ax,current_x,current_y,current_z,'.r')
     drawnow
    view([-125.02 17.07])
    axis([-0.6 0.6 -0.6 0.6 -0.25 1])
%     waitfor(r);
end 
%% plotting workspace with toolbox
figure(3)
axis manual
axis([-0.6 0.6 -0.6 0.6 -0.25 1])
view([-125.02 17.07])
framesPerSecond = 20;
i=1;
for i2 = linspace(jnt2.PositionLimits(1),jnt2.PositionLimits(2),30)
for i3 = linspace(jnt3.PositionLimits(1),jnt3.PositionLimits(2),30)
for i4 = linspace(jnt4.PositionLimits(1),jnt4.PositionLimits(2),10)
for i5 = linspace(jnt5.PositionLimits(1),jnt5.PositionLimits(2),30)
    current_config = [0 ,i2, i3-pi, i4, i5 ,0 ];
%     ax= show(robot,current_config,'FastUpdate',1,'PreservePlot',0);
    current_Tee = getTransform(robot,current_config,'Link6','base');
    current_x(i) = current_Tee(1,4);
    current_y(i) = current_Tee(2,4);
    current_z(i) = current_Tee(3,4);
    i = i+1
end 
end
end
end
k = boundary([current_x',current_y',current_z']);
trisurf(k,current_x,current_y,current_z,'LineStyle',':')
axis equal
%% checking Tee
config1 = [25 10 20 20 10 20]/180 * pi ;

ans1 = vpa(subs(T_ee,[tta1 tta2 tta3 tta4 tta5 tta6],config1),5)
config1 = [25 10 20-180 20 10 20]/180 * pi ;
getTransform(robot,config1,'Link6')
%% computing jacobian matrix
% using straight diff
Jv = simplify(jacobian(T_ee(1:3,4),[tta1 tta2 tta3 tta4 tta5 tta6]));
R_ee = T_ee(1:3,1:3);
syms 
%% plotting workspace using code
figure(4)
axis manual
axis([-0.6 0.6 -0.6 0.6 -0.25 1])
view([-125.02 17.07])
i=1;
T_ee_func = matlabFunction(T_ee, 'vars', [tta1 tta2 tta3 tta4 tta5 tta6]);
for i1 = linspace(-170/180*pi,170/180*pi,30)
for i2 = linspace(-190/180*pi,45/180*pi,30)
for i3 = linspace(-29/180*pi,256/180*pi,30)
for i4 = linspace(-190/180*pi,190/180*pi,10)
for i5 = linspace(-120/180*pi,120/180*pi,30)
    current_Tee  = T_ee_func(i1,i2,i3,i4,i5,0);
    current_x(i) = current_Tee(1,4);
    current_y(i) = current_Tee(2,4);
    current_z(i) = current_Tee(3,4);
    i = i+1
end 
end
end
end
end
% k = boundary([current_x',current_y',current_z']);
% trisurf(k,current_x,current_y,current_z,'LineStyle',':')
plot3(current_x,current_y,current_z,'.')
axis equal

