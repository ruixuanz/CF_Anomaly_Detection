
% clear all
% clc
% 
% N = 31;
% 
% alpha = [0.7, 0.2, 0.1];
% beta = [0.7, 0.2, 0.1];
% v_e = 22.47;
% s_e = 30;
% s_0 = 2;
% v_0 = 33.33;
% T = 1.1;
% a = 1;
% b = 2;
% 
% % define attack parameters A, B, C
% 
% % attack parameter
% 
% Attack_para = [0, 0, 0
%                -3, -5, -11];
% 
% 
% 
% % define the frequency step
% N_w = 10000; % number of steps
% Omega = linspace(0, pi, N_w);
% 
% 
% Tau_1 = linspace(0, 3, N);
% Tau_2 = linspace(3, 0, N);
% 
% Tau_analysis_result = zeros(length(Tau_1),length(Tau_1), 2);
% 
% 
% for i = 1:length(Tau_2)
%     tau_2 = Tau_2(i);
%     
%     for j = 1:length(Tau_1)
%         tau_1 = Tau_1(j);
%         
%         
%         for k = 1:2
%             
%             A =  Attack_para(k,1);
%             B =  Attack_para(k,2);
%             C =  Attack_para(k,3);
%             
%             result_tilde = zeros(4, N);
%             count = 1;
%             for omega = Omega
%                 [f_vn_tilde, f_dvn_tilde, f_sn_tilde] = derivative(A,B,C,alpha,beta,v_e,s_e,v_0,s_0,T,a,b);
%                 result_tilde(:,count) = eig(TF(1i*omega, f_vn_tilde, f_dvn_tilde, f_sn_tilde, alpha, beta, tau_1, tau_2));
%                 count = count+1;
%             end
%             norm_val_tilde = sort(abs(result_tilde),'descend');
%             [max_norm_val_tilde,index_tilde]=max(norm_val_tilde);
%             Tau_analysis_result(i,j,k) = max(max_norm_val_tilde);
%         end
%     end
%     
%     display(i)
% end


Result_plot = Tau_analysis_result;
% Result_plot = importdata('Result.mat'); 
[m ,n, p] = size(Result_plot);
minColorLimit = 1;                   % determine colorbar limits from data
sorted_result = sort(Result_plot(:));
maxColorLimit = sorted_result(round(length(sorted_result)*0.95));

title_list = ["Originally stable" "Originally unstable"];

fig = figure(1);
for fig_num  = 1:p
    sph{fig_num} = subplot(1,2,fig_num,'Parent',fig);
    R = Result_plot(:,:,fig_num);
    R(R<=1.0001) = NaN;
    h = imagesc(sph{fig_num}, R);
    
    set(sph{fig_num}, 'xtick', 1:n, 'xticklabel', ["0" " " " " " " " " "0.5" " " " " " " " " "1" " " " " " " " " "1.5" " " " " " " " " "2" " " " " " " " " "2.5" " " " " " " " " "3"], ...
    'ytick', 1:m, 'yticklabel', ["3" " " " " " " " " "2.5" " " " " " " " " "2" " " " " " " " " "1.5" " " " " " " " " "1" " " " " " " " " "0.5" " " " " " " " " "0"])
    set(h,'alphadata',~isnan(R));
    
    title(title_list(p), 'interpreter','latex');

    caxis(sph{fig_num},[minColorLimit,maxColorLimit]);
end

h = axes(fig,'visible','off'); 
h.Title.Visible = 'on';
h.XLabel.Visible = 'on';
h.YLabel.Visible = 'on';
xlabel(h,['$\tau_1$', ' /s'],'FontWeight','bold', 'interpreter','latex');
ylabel(h,['$\tau_2$', ' /s'],'FontWeight','bold', 'interpreter','latex', 'Position', [-0.8, 0.5, 0]);
title(h,'The largest eigenvalue $\lambda$ of transfer matrix $\tilde{P}_n$','interpreter','latex', 'fontsize', 12, 'Position', [0.5, 1.2, 0]);

c = colorbar(h,'Position',[0.93 0.168 0.022 0.7]);  % attach colorbar to h
colormap(c,'jet')
caxis(h,[minColorLimit,maxColorLimit]);
axis equal


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

