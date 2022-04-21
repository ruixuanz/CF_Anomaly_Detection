%%

clear all
clc

%% Define the parameter

alpha = [0.7, 0.2, 0.1];
beta = [0.7, 0.2, 0.1];
tau_1 = 0.2;
tau_2 = 0.5;
v_e = 22.47;
s_e = 30;
s_0 = 2;
v_0 = 33.33;
T = 1.1;
a = 1;
b = 2;

%% Define the intermediate terms

f_vn = -4*a*v_e^3/v_0^4 - 2*a*T*(s_0+T*v_e)/s_e^2;
f_dvn = -beta(1)*v_e*(s_0+T*v_e)*sqrt(a/b)/s_e^2;
f_sn = 2*alpha(1)*a*(s_0+T*v_e)^2/s_e^3;

f_normal = @(lambda, omega) lambda^2 - lambda*(f_vn*exp(-lambda*tau_1) + f_dvn*beta(1)*(exp(-lambda*tau_1)-exp(-1i*omega-lambda*tau_2)) ...
    + (exp(1i*omega)-1)*f_dvn*(beta(2)*exp(-2i*omega-lambda*tau_2)+beta(3)*exp(-3i*omega-lambda*tau_2))) ...
    - (exp(-1i*omega)-1)*f_sn*alpha(1)*exp(-lambda*tau_1) - (exp(-1i*omega)-1)*f_sn*(alpha(2)*exp(-1i*omega-lambda*tau_2)+alpha(3)*exp(-2i*omega-lambda*tau_2));

% noise terms
A = 0;
B = 25;
C = -10;

f_vn_tilde = -4*a*(v_e+A)^3/v_0^4 - 2*a*(T+C/sqrt(a*b))*(s_0+T*(v_e+A)+C*(v_e+A)/(2*sqrt(a*b)))/(s_e+B)^2;
f_dvn_tilde = -beta(1)*sqrt(a/b)*(v_e+A)*(s_0+T*(v_e+A)+C*(v_e+A)/2*sqrt(a*b))/(s_e+B)^2;
f_sn_tilde = 2*alpha(1)*a*(s_0+T*(v_e+A)+C*(v_e+A)/2*sqrt(a*b))^2/(s_e+B)^3;

f_noise = @(lambda, omega) lambda^2 - lambda*(f_vn_tilde*exp(-lambda*tau_1) + f_dvn_tilde*beta(1)*(exp(-lambda*tau_1)-exp(-1i*omega-lambda*tau_2)) ...
    + (exp(1i*omega)-1)*f_dvn_tilde*(beta(2)*exp(-2i*omega-lambda*tau_2)+beta(3)*exp(-3i*omega-lambda*tau_2))) ...
    - (exp(-1i*omega)-1)*f_sn_tilde*alpha(1)*exp(-lambda*tau_1) - (exp(-1i*omega)-1)*f_sn_tilde*(alpha(2)*exp(-1i*omega-lambda*tau_2)+alpha(3)*exp(-2i*omega-lambda*tau_2));

%% Given an omega, solve for lambda

vpasolve(f_normal(sym('lambda'),0))

%% 

N = 500; % number of steps
Omega = linspace(0, pi, N);
result_normal = zeros(1, N);
result_noise = zeros(1, N);
syms lambda

i = 1;
for omega = Omega
    result_normal(i) = real(vpasolve(f_normal(lambda, omega)));
    result_noise(i) = real(vpasolve(f_noise(lambda, omega)));
    i = i+1;
end

%% Plot

plot(Omega, result_normal)
hold on
plot(Omega, result_noise)
grid on
grid minor

%% Plot

p = 0.9;
plot(Omega, p*result_normal+(1-p)*result_noise)
grid on
grid minor
disp(max(p*result_normal+(1-p)*result_noise))