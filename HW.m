clear;
close all;
%% Data
rho = 7600;                     % kg/m^3
E = 210e9;                      % Pa (GPa to Pa)
nu = 0.3;                       % Poisson's ratio
alpha = 1.2e-5;                 % 1/K
k = 45;                         % W/mK
T0=20;                          % °C

% A-1
a = 100e-3;                     % m (converted from mm)
b = 250e-3;                     % m
Ti = 40;                        % deg C
qi = 620;                       % W/m^2

% B-5
p0 = 10e5;                      % Pa (bar to Pa)
ho = 8;                         % W/m^2K
To = 25;                        % deg C
qo = 660;                       % W/m^2

% C-4
hi = 10;                        % W/m^2K

%% #2
ab=(b-a)*1000+1;
r = linspace(a,b,ab);  

% Temperature distribution function
term1 = (b * ho * (Ti - To) * log(a)) / (k + b * ho * log(b/a));
term2 = (b * ho * (Ti - To)) / (-k + b * ho * log(a/b));
theta = Ti + term1 + term2*log(r);

% Plot
figure();
hold on;
plot(r, theta, 'LineWidth', 2)
xlabel('Radius r [m]')
ylabel('T(r) [°C]')
title('Radial Temperature Distribution')
grid on
r_loop=r;
filename = 'temperature.dat';
fileID = fopen(filename, 'r');
data = fscanf(fileID, '%f', [2, Inf]);
fclose(fileID);

data = data';

x = data(:, 1);
y = data(:, 2);
plot(x,y,'o');
legend("Analytical","Numerical")
%% #3
syms A_prime B_prime r
E_prime=E/(1-nu^2);
fun=@(r) r.*(Ti+(b .* ho .* (Ti - To) .* log(a)) ./ (k + b .* ho .* log(b/a))+(b .* ho .* (Ti - To)) ./ (-k + b .* ho .* log(a/b)) .* log(r));
eq1= A_prime-B_prime/a^2==-p0;
eq2= A_prime-B_prime/b^2==((E_prime*alpha)/b^2)*integral(fun,a,b);
[M1,M2] = equationsToMatrix([eq1, eq2], [A_prime,B_prime]);
sol=vpa(linsolve(M1,M2));
A_prime=double(sol(1));
B_prime=double(sol(2));

for i=1:1:ab
sigma_r(i)=-((E_prime*alpha)/r_loop(i)^2)*integral(fun,a,r_loop(i))+A_prime-B_prime/(r_loop(i)^2);
%sigma_teta(i)=((E_prime*alpha)/r_loop(i)^2)*integral(fun,a,r_loop(i))+A_prime+B_prime/(r_loop(i)^2)-E_prime*alpha*theta(i);
sigma_teta(i)=-sigma_r(i)-E_prime*alpha*theta(i)+2*A_prime;
sigma_z(i)=nu*(sigma_r(i)+sigma_teta(i));
end

%filename = 'stress.dat'; %1
%filename= 'stress_closed.dat'; %2
filename = 'stress_what.dat'; %1

fileID = fopen(filename, 'r');
data = fscanf(fileID, '%f', [4, Inf]);
fclose(fileID);

% Transpose the matrix to have columns properly
data = data';

% Separate into two arrays
r = data(:, 1);
sigma_r_ansys = data(:, 2);
sigma_theta_ansys = data(:, 3);
sigma_z_ansys = data(:, 4);

%% Plot
figure()
grid on;
title("Stress components")
plot(r_loop*1000,sigma_r/10^6,'r-','LineWidth',1.5);
xlabel("r [mm]"); ylabel("σ [MPa]")
hold on;
plot(r_loop*1000,sigma_teta/10^6,'g-','LineWidth',1.5);
plot(r_loop*1000,sigma_z/10^6,'b-','LineWidth',1.5);

plot(r*1000,sigma_r_ansys,'o','Color','r');
plot(r*1000,sigma_theta_ansys,'o','Color','g');
plot(r*1000,sigma_z_ansys,'o','Color','b');

legend("σ_r analytical","σ_θ analytical","σ_z analytical","σ_r numerical","σ_θ numerical","σ_z numerical")
grid on;
