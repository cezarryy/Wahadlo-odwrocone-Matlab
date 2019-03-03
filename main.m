%Projekt wahad³a odwróconego

%Parametry wahad³a:
M = 0.5; %Masa wózka
m = 0.2; %Masa wahad³a
l = 0.3; %Dlugoœæ od mocowania do œrodka ciê¿koœci wahad³a
I = 0.006; %Moment bezw³adnoœci wahad³a
b = 0.1; %Wspó³czynnik tarcia wózka
g = 9.80665; %Przyspieszenie ziemskie

% parametry symulacji
tSim = 5;
global h;
h = 0.01;
t = (0:h:tSim)';
tt = numel(t);
setTheta = [t,ones(tt,1) * 0];
set_x= [t,ones(tt,1) * 0];


%Warunki pocz¹tkowe:
x0 = 0;
dx0 = 0;
theta0 = pi/10;
dtheta0 = 0;

%Macierze stanu:
A = [0 1 0 0;
    0 -b/M -(m*g)/M 0;
    0 0 0 1;
    0 b/(M*l) ((M+m)*g)/(M*l) 0];
B = [0; 1/M; 0; -1/(M*l)];
C = [1 0 0 0;
    0 1 0 0
    0 0 1 0
    0 0 0 1];
D = [0; 0; 0; 0];

%Utworzenie uk³adu zmiennych stanu;
s_system = ss(A,B,C,D);


%Utworzenie macierzy warunków pocz¹tkowych: x, dx, theta, dtheta
state0 = [x0, dx0, theta0, dtheta0];
angle0 = [dtheta0, theta0];
%Transmitancja:
[L,Mian] = ss2tf(A,B,C,D);

%Dyskretyzacja
d_system = c2d(s_system,0.01);
[Ad,Bd,Cd,Dd] = ssdata(d_system);

[Ld,Md] = ss2tf(Ad,Bd,Cd,Dd);

Lpd = Ld(1,:);
Lwd = Ld(3,:);
Gp = tf(Lpd,Md,0.01);
[L_pol,M_pol] = tfdata(Gp,'v');
Gw = tf(Lwd,Md,0.01,'variable','z^-1')
[L_wych,M_wych] = tfdata(Gw,'v');

%Regulator LQR:
Q = diag([100 0 3600 0]);
R = 1;
K = dlqr(Ad,Bd,Q,R);
%sim('pendulum_PID.slx');

%Transmitancja z kartki:
a1 = (m.*l - (M+m)*l);
a2 = -1.*b.*l;
a3 = (M+m).*g;
a4 = b.*g;

Num = [-1 0];
Den = [a1 a2 a3 a4];
disc_g = tf(Num,Den,0.01,'variable','z^-1');
[Num_d, Den_d] = tfdata(disc_g,'v');

%Synteza deadbeat:
q0 = 1./(Lwd(2) + Lwd(3) + Lwd(4) + Lwd(5));
L_db(1)=q0;
for i=2:(size(Md,2)-1)
    L_db(i) = q0.*Md(i+1);
end
M_db(1)=1;
for i=2:(size(Lwd,2)-1)
    M_db(i) = -1.*q0.*Lwd(i+1);
end
G_db = tf(L_db,M_db,0.01,'variable','z^-1')