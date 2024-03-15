%% If you want to use this code, the citation of our paper is needed

clc
close all
clear all

%% Geometrical parameters
R0 = 200;         % [mm] The radius of the centroid axis of the ring model
E2 = 10000;         % [MPa] The elastic modulus of the ring in the longitudinal direaction
G = 5;              % [MPa] The shear modulus of the ring
b = 60;           % [mm] The width of the ring
h = 20;           % [mm] The thickness of the ring
I2 = b*h^3/12;    % The moment of intertia of the cross-section of the ring
A2 = b*h;
% delta0 = 20;      % The vertical deflection of the plate
delta1 = 20;
R_c = 300;        % The radius of the curved road

% Error = abs(ur_bar_EA_PairForce(0,0,R0,E2,I2,A2))*1e-03;
Error = 8e-07;
k_gen = 0.1;   % Assuming an relatively small stiffness between the road and ring
N_ang = 1800; % The number of dividing regions (Due to the symmetry, N_ang should be even)
F_i_0 = zeros(N_ang+1,length(delta1));
ang_i_0 = zeros(N_ang+1,length(delta1));
theta_di_0 = zeros(N_ang+1,length(delta1));
F_ver = zeros(length(delta1),1);
X_Rc0_0 = zeros(N_ang+1,length(delta1));
Y_Rc0_0 = zeros(N_ang+1,length(delta1));

% Iteration starting
for DD = 1:length(delta1)
delta0 = delta1(DD);

%% Define the origin of the concave surface
d_c = (R0+h/2)-R_c-delta0; % The distance between the road center and ring center
x_c = d_c; 
y_c = 0;          % The position of the circular center of the road

%% Calculating the force applied to the contact region

Cnt_ang0 = pi/3;

% N_ang = 1000; % The number of dividing regions (Due to the symmetry, N_ang should be even)
D_ang = 2*Cnt_ang0/N_ang;

ang_i = (-Cnt_ang0:D_ang:Cnt_ang0)'; %The rotation angle of the force pair

ang_i_0(:,DD) = ang_i;

% Road function for calculating geometry errors (using polynomial fitting)
d1 = zeros(length(ang_i),1);
ang_R0 = zeros(length(ang_i),1);

x_Rc0 = zeros(length(ang_i),1);
y_Rc0 = zeros(length(ang_i),1);

% Undeformed ring centroid
r0 = (R0+h/2)*ones(length(ang_i),1);
Cnt_ur0 = zeros(length(ang_i),1);
Cnt_ur = zeros(length(ang_i),1);
Cnt_UR = zeros(length(ang_i),1);
Cnt_ut0 = zeros(length(ang_i),1);
Cnt_ut = zeros(length(ang_i),1);
Cnt_UT = zeros(length(ang_i),1);
Cnt_fi0 = zeros(length(ang_i),1);
Cnt_fi = zeros(length(ang_i),1);
Cnt_FI = zeros(length(ang_i),1);
D_Fr = zeros(length(ang_i),1);
D_Fr0 = zeros(length(ang_i),1);
F_i = zeros(length(ang_i),1);
r_deform = r0;
theta_di = zeros(length(ang_i),1) + ang_i;
Road = zeros(length(theta_di),1);
theta_Road = zeros(length(Road),1);

for j0 = 1:20000 % The maximum iteration step

for i0 = 1:length(theta_di)
    d1(i0) = abs(d_c)*cos(pi-abs(theta_di(i0)))+sqrt(d_c^2*cos(pi-abs(theta_di(i0)))^2+R_c^2-d_c^2);
    ang_R0(i0) = sign(theta_di(i0))*acos((R_c^2+d_c^2-d1(i0)^2)/(2*R_c*abs(d_c)));
    x_Rc0(i0) = x_c + R_c*cos(ang_R0(i0));
    y_Rc0(i0) = y_c + R_c*sin(ang_R0(i0));
end

[theta_Road, Road] = cart2pol(x_Rc0,y_Rc0);

% The distributing force along the contact region
for i2 = 1:N_ang+1
    D_Fr0(i2) = k_gen*(r_deform(i2) - Road(i2));
    Delta_F = F_i(i2) + D_Fr0(i2);
    if Delta_F > 0
        D_Fr(i2) = D_Fr0(i2); % This step ensures the direction of the pressure
    else
        D_Fr(i2) = 0;
    end
end

F_i = F_i + D_Fr;

% The deformation of the penetration part due to the distributed load
for i3 = 1:length(theta_di)
    for i4 = 1:length(D_Fr)
        Cnt_ur0(i4) = D_Fr(i4)*(R0+h/2)*D_ang*ur_bar_EA_GA_PairForce(theta_di(i3),theta_di(i4),R0,E2,G,I2,A2);
        Cnt_ut0(i4) = D_Fr(i4)*(R0+h/2)*D_ang*ut_bar_EA_GA_PairForce(theta_di(i3),theta_di(i4),R0,E2,G,I2,A2);
        Cnt_fi0(i4) = D_Fr(i4)*(R0+h/2)*D_ang*ufi_bar_EA_GA_PairForce(theta_di(i3),theta_di(i4),R0,E2,G,I2,A2);
    end
    Cnt_ur(i3) = sum(Cnt_ur0);
    Cnt_ut(i3) = sum(Cnt_ut0);
    Cnt_fi(i3) = sum(Cnt_fi0);
end

Cnt_UR = Cnt_UR + Cnt_ur;
Cnt_UT = Cnt_UT + Cnt_ut;
Cnt_FI = Cnt_FI + Cnt_fi;

% The deformed cylindrical coordinates of the penetration part % The total deformation of the contact region
for i0 = 1:N_ang+1
    r_deform(i0) = sqrt((R0+h/2+Cnt_UR(i0))^2+(Cnt_UT(i0)+h/2*Cnt_FI(i0))^2); % Define the initial deformed matrix
    theta_di(i0) = ang_i(i0) + atan((Cnt_UT(i0)+h/2*Cnt_FI(i0))/(R0+h/2+Cnt_UR(i0)));
end

% The error between the road and the ring
D_err = zeros(length(ang_i),1);
for i5 = 1:N_ang+1
    D_err(i5) = r_deform(i5) - Road(i5);  
end

% Stop the iteration if the change of the deformation is quite small
if max(abs(Cnt_ur)) < Error      %8e-7 2.5e-6 need to be modified
    break
else
    Cnt_ur = zeros(length(ang_i),1);
    Cnt_ut = zeros(length(ang_i),1);
    Cnt_fi = zeros(length(ang_i),1);
    continue
end

end

F_i_0(:,DD) = F_i;
X_Rc0_0(:,DD) = x_Rc0;
Y_Rc0_0(:,DD) = y_Rc0;

% Calculate the total vertical load
F_v = zeros(length(F_i),1);

for i6 = 1:length(F_i)
    F_v(i6) = F_i(i6)*cos(ang_i(i6))*(R0+h/2)*D_ang;
end

F_ver(DD) = sum(F_v);  % The vertical reaction force on the plate
theta_di_0(:,DD) = theta_di;

end