%% Define the parameter

alpha = [0.7, 0.2, 0.1];
beta = [0.7, 0.2, 0.1];
tau_1 = 0.5;
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

f = @(lambda, omega, tau_1, tau_2, v_e, s_e, s_0, v_0, T, a, b) lambda^2 - lambda*(f_vn*exp(-lambda*tau_1) + f_dvn*beta(1)*(exp(-lambda*tau_1)-exp(-1i*omega-lambda*tau_2)) ...
    + (exp(1i*omega)-1)*f_dvn*(beta(2)*exp(-2i*omega-lambda*tau_2)+beta(3)*exp(-3i*omega-lambda*tau_2))) ...
    - (exp(-1i*omega)-1)*f_sn*alpha(1)*exp(-lambda*tau_1) - (exp(-1i*omega)-1)*f_sn*(alpha(2)*exp(-1i*omega-lambda*tau_2)+alpha(3)*exp(-2i*omega-lambda*tau_2));

%% Given an omega, solve for lambda

vpasolve(f(sym('lambda'),0, tau_1, tau_2, v_e, s_e, s_0, v_0, T, a, b))

%% 

N = 1000; % number of steps
Omega = linspace(0,2*pi, N);
result=zeros(1, N);
syms lambda

i = 1;
for omega = Omega
    result(i) = real(vpasolve(f(lambda, omega, tau_1, tau_2, v_e, s_e, s_0, v_0, T, a, b)));
    i = i+1;
end

%% Plot

plot(Omega, result)