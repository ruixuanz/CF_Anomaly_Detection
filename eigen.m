% %
% clear all
% clc
% 
% % Define the parameter
% 
% alpha = [0.7, 0.2, 0.1];
% beta = [0.7, 0.2, 0.1];
% tau_1 = 0;
% tau_2 = 0.5;
% v_e = 22.47;
% s_e = 30;
% s_0 = 2;
% v_0 = 33.33;
% T = 1.1;
% a = 1;
% b = 2;
% 
% % Solve for eigenvalues and their norms
% 
% N = 10000; % number of steps
% Omega = linspace(0, pi, N);
% result=zeros(4, N);
% result_tilde=zeros(4, N);
% syms lambda
% 
% i = 1;
% for omega = Omega
%     [f_vn, f_dvn, f_sn] = derivative(0,0,0,alpha,beta,v_e,s_e,v_0,s_0,T,a,b);
%     [f_vn_tilde, f_dvn_tilde, f_sn_tilde] = derivative(0,25,-10,alpha,beta,v_e,s_e,v_0,s_0,T,a,b);
%     result(:,i) = eig(TF(1i*omega, f_vn, f_dvn, f_sn, alpha, beta, tau_1, tau_2));
%     result_tilde(:,i) = eig(TF(1i*omega, f_vn_tilde, f_dvn_tilde, f_sn_tilde, alpha, beta, tau_1, tau_2));
%     i = i+1;
% end
% 
% norm_val = sort(abs(result),'descend');
% norm_val_tilde = sort(abs(result_tilde),'descend');
% [max_norm_val,index]=max(norm_val);
% [max_norm_val_tilde,index_tilde]=max(norm_val_tilde);
% 
% % find the critical detection rate
% 
% P = linspace(0, 1, 10000);
% 
% for p = P
% 
%     if max(p*max_norm_val + (1-p)*max_norm_val_tilde) < 1
%         continue;
%     else
%         critical_p = p;
%     end
% end
% plot(Omega, max_norm_val, '--', 'color','b', 'LineWidth',2)
% hold on
% plot(Omega, max_norm_val_tilde, '--', 'color','r', 'LineWidth',2)
% hold on
% plot(Omega, critical_p*max_norm_val+(1-critical_p)*max_norm_val_tilde, '--', 'color','g', 'LineWidth',2)
% 
% legend('Successfully detect','Fail to detect','Stochastic max eigenvalue');
% 
% grid on
% grid minor
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% A_vector = linspace(-15, 15, 16);
% B_vector = linspace(15, -15, 16);
% C_vector = linspace(-15, 15, 16);
% Result = zeros(length(A_vector),length(B_vector),length(C_vector));
% 
% for i = 1:length(B_vector)
%     B = B_vector(i);
%     for j = 1:length(C_vector)
%         C = C_vector(j);
%         for k = 1:length(A_vector)
%             A = A_vector(k);
%             result_tilde = zeros(4, N);
%             count = 1;
%             for omega = Omega
%                 [f_vn_tilde, f_dvn_tilde, f_sn_tilde] = derivative(A,B,C,alpha,beta,v_e,s_e,v_0,s_0,T,a,b);
%                 result_tilde(:,count) = eig(TF(1i*omega, f_vn_tilde, f_dvn_tilde, f_sn_tilde, alpha, beta, tau_1, tau_2));
%                 count = count+1;
%             end
%             norm_val_tilde = sort(abs(result_tilde),'descend');
%             [max_norm_val_tilde,index_tilde]=max(norm_val_tilde);
%             Result(i,j,k) = max(max_norm_val_tilde);
%         end
%     end
% end

% Sub-plot

Result_plot = importdata('Result.mat');
figure(1);
for fig  = 1:16
    subplot(4,4,fig);
    matrix_plot = Result_plot(:,:,fig);
    matrix_plot(matrix_plot>1.0001)=0 ;
    h = heatmap(matrix_plot,'XData', linspace(-15, 15, 16),'fontsize',10,'YData', linspace(15, -15, 16));
    h.XLabel = '$\Delta\tilde{v}$ m/s';
    h.YLabel = '$\tilde{s}$ m';
    h.Title = ['$\tilde{v}_n$ = ',num2str(A_vector(fig)), 'm/s'];
    h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
    h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
    h.NodeChildren(3).Title.Interpreter = 'latex';
end

% result_tilde = zeros(4, N);
% count = 1;
% for omega = Omega
%     [f_vn_tilde, f_dvn_tilde, f_sn_tilde] = derivative(-3,7,1,alpha,beta,v_e,s_e,v_0,s_0,T,a,b);
%     result_tilde(:,count) = eig(TF(1i*omega, f_vn_tilde, f_dvn_tilde, f_sn_tilde, alpha, beta, tau_1, tau_2));
%     count = count+1;
% end
% norm_val_tilde = sort(abs(result_tilde),'descend');
% [max_norm_val_tilde,index_tilde]=max(norm_val_tilde);
% max(max_norm_val_tilde)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the functions to compute the intermediate terms

function [f_vn, f_dvn, f_sn] = derivative(A,B,C,alpha,beta,v_e,s_e,v_0,s_0,T,a,b)
    f_vn = -4*a*(v_e+A)^3/v_0^4 - 2*a*(T+C/sqrt(a*b))*(s_0+T*(v_e+A)+C*(v_e+A)/(2*sqrt(a*b)))/(s_e+B)^2;
    f_dvn = -beta(1)*sqrt(a/b)*(v_e+A)*(s_0+T*(v_e+A)+C*(v_e+A)/2*sqrt(a*b))/(s_e+B)^2;
    f_sn = 2*alpha(1)*a*(s_0+T*(v_e+A)+C*(v_e+A)/2*sqrt(a*b))^2/(s_e+B)^3;
end

function P_hat = TF(s, f_vn, f_dvn, f_sn, alpha, beta, tau_1, tau_2)

    T_1 = (s*(beta(2)-beta(1))*f_dvn+(alpha(1)-alpha(2))*f_sn)*exp(-s*tau_2)/(s^2-s*(f_vn+beta(1)*f_dvn)*exp(-s*tau_1)+alpha(1)*f_sn*exp(-s*tau_1));
    T_2 = (s*(beta(3)-beta(2))*f_dvn+(alpha(2)-alpha(3))*f_sn)*exp(-s*tau_2)/(s^2-s*(f_vn+beta(1)*f_dvn)*exp(-s*tau_1)+alpha(1)*f_sn*exp(-s*tau_1));
    T_3 = (-s*beta(3)*f_dvn+alpha(3)*f_sn)*exp(-s*tau_2)/(s^2-s*(f_vn+beta(1)*f_dvn)*exp(-s*tau_1)+alpha(1)*f_sn*exp(-s*tau_1));
    
    P_hat = [T_1 T_2 T_3 0;
             1   0   0   0
             0   1   0   0
             0   0   1   0];
end



