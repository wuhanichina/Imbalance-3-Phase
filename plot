%% Codes for "Analytical Energy Loss Estimation for Arbitrarily Line Currents"
% ����ÿ������ÿ��ʱ�̵ĵ�������ͬһ���ֲ����÷ֲ�������GMM��ϣ�����һ������������ʱ�̵ĵ����Ϳ�����һ���������x��ʾ��
% ע��������ﲢ����˵�������x���Բ��Ϊk�����и�˹�ֲ�������������������,E(x)=pi1*g1+pi2*g2,...,+pik*gk.
% ʵ������������õ���GMM�ڷ��������ϵĽ��ͣ�Ҳ����˵�Ѳɼ����ĵ����ֳ�k�ݣ�ÿ�ݶ�����ĳһ����˹�ֲ���
% �����ڼ���I^2Rʱ��������ɢ�ģ����������ɸ������������I^2R��Ӽ��ɵõ�������ʵ��ֵ��
%
% A part of my  "Project C", which contains many works that relates to the  convolution and copula.
% Han Wu, Hohai University, China
% E-mail: wuhan@hhu.edu.cn 
% Working on this item during <27th Mar 2020, 28th Mar 2020 >% 
% Copyrights reserved
clc;clear;
%% ʵ�����ݲ��� Now we focus on a real case
load("loaddata.mat");% ���ݴ��ھ���three_phace�� load data
options = statset('Display','final');
m = 2;% Gaussian component�������� 
ra = 0.01; rb = 0.02; rc = 0.03; % ����·���� load line resistance
del_t = 10*60; % sec ʵ������Ϊ10����һ���㣬����Ϊ��λ����������������λJ
% ��GMM����������� fit data use GMM
gmm_fit_3P=fitgmdist(three_phase,m,'Options',options);

% �������ʱ��ʵ�ʵ�����ƽ�� square sum of feeder current Phase A
Ia_square= three_phase(:,1).*three_phase(:,1);
Ia_square_sum = sum(Ia_square);

% square sum of feeder current Phase B
Ib_square= three_phase(:,2).*three_phase(:,2);
Ib_square_sum = sum(Ib_square);

%square sum of feeder current Phase C
Ic_square= three_phase(:,3).*three_phase(:,3);
Ic_square_sum = sum(Ic_square);

% �����ʱ����ʵ��������Actual energy loss
Energy_loss = del_t*(Ia_square_sum*ra+Ib_square_sum*rb+Ic_square_sum*rc);
Energy_loss_A = del_t*Ia_square_sum*ra;
Energy_loss_B = del_t*Ib_square_sum*rb;
Energy_loss_C = del_t*Ic_square_sum*rc;
%% ������忨���ֲ��µľ�ֵCompute expected energy loss for three-phase unbalanced network
% ����һ���� First for the whole energy loss
mean_all_com1 = gmm_fit_3P.mu(1,:);
mean_all_com2 = gmm_fit_3P.mu(2,:);
Three_phase_cov1 = gmm_fit_3P.Sigma(:,:,1);
Three_phase_cov2 = gmm_fit_3P.Sigma(:,:,2);
com_pi = gmm_fit_3P.ComponentProportion;
A = [ra,rb,rc].*del_t;
chi_loss_mean_1 = trace(A.*Three_phase_cov1)+A*(mean_all_com1.*mean_all_com1)'; 
chi_loss_mean_2 = trace(A.*Three_phase_cov2)+A*(mean_all_com2.*mean_all_com2)';
Chi_Energy_loss = (chi_loss_mean_1*com_pi(1,1)+chi_loss_mean_2*com_pi(1,2))*length(three_phase);
% ��� Estermation error
Chi_Error = abs(Chi_Energy_loss-Energy_loss)/Energy_loss;
disp Erros_real_All=  
disp(Chi_Error*100);

% ���൥���� Then for energy loss on each phase independently
gmm_fit_Pa=fitgmdist(three_phase(:,1),m,'Options',options);
gmm_fit_Pb=fitgmdist(three_phase(:,2),m,'Options',options);
gmm_fit_Pc=fitgmdist(three_phase(:,3),m,'Options',options);

% A�� Phase A
mean_A_com1 = gmm_fit_Pa.mu(1,:);
mean_A_com2 = gmm_fit_Pa.mu(2,:);
phase_A_cov1 = gmm_fit_Pa.Sigma(:,:,1);
phase_A_cov2 = gmm_fit_Pa.Sigma(:,:,2);
com_pi_A = gmm_fit_Pa.ComponentProportion;

chi_loss_mean_A_1 = (trace(ra.*del_t.*phase_A_cov1)+ra.*del_t*(mean_A_com1.*mean_A_com1)');
chi_loss_mean_A_2 = (trace(ra.*del_t.*phase_A_cov2)+ra.*del_t*(mean_A_com2.*mean_A_com2)');
Chi_Energy_loss_A = (chi_loss_mean_A_1*com_pi_A(1,1)+chi_loss_mean_A_2*com_pi_A(1,2))*length(three_phase);
% ���Estimated energy loss on Phase A
Chi_Error_A = abs(Chi_Energy_loss_A-Energy_loss_A)/Energy_loss_A;
disp Erros_real_A=  
disp(Chi_Error_A*100);

% B�� Phase B
mean_B_com1 = gmm_fit_Pb.mu(1,:);
mean_B_com2 = gmm_fit_Pb.mu(2,:);
phase_B_cov1 = gmm_fit_Pb.Sigma(:,:,1);
phase_B_cov2 = gmm_fit_Pb.Sigma(:,:,2);
com_pi_B = gmm_fit_Pb.ComponentProportion;

chi_loss_mean_B_1 = (trace(rb.*del_t.*phase_B_cov1)+rb.*del_t*(mean_B_com1.*mean_B_com1)');
chi_loss_mean_B_2 = (trace(rb.*del_t.*phase_B_cov2)+rb.*del_t*(mean_B_com2.*mean_B_com2)');
Chi_Energy_loss_B = (chi_loss_mean_B_1*com_pi_B(1,1)+chi_loss_mean_B_2*com_pi_B(1,2))*length(three_phase);
% ��� Estimated energy loss on Phase B
Chi_Error_B = abs(Chi_Energy_loss_B-Energy_loss_B)/Energy_loss_B;
disp Erros_real_B=  
disp(Chi_Error_B*100);

% C�� Phase C
mean_C_com1 = gmm_fit_Pc.mu(1,:);
mean_C_com2 = gmm_fit_Pc.mu(2,:);
phase_C_cov1 = gmm_fit_Pc.Sigma(:,:,1);
phase_C_cov2 = gmm_fit_Pc.Sigma(:,:,2);
com_pi_C = gmm_fit_Pc.ComponentProportion;

chi_loss_mean_C_1 = (trace(rc.*del_t.*phase_C_cov1)+rc.*del_t*(mean_C_com1.*mean_C_com1)');
chi_loss_mean_C_2 = (trace(rc.*del_t.*phase_C_cov2)+rc.*del_t*(mean_C_com2.*mean_C_com2)');
Chi_Energy_loss_C = (chi_loss_mean_C_1*com_pi_C(1,1)+chi_loss_mean_C_2*com_pi_C(1,2))*length(three_phase);
% ��� Estimated energy loss on Phase C
Chi_Error_C = abs(Chi_Energy_loss_C-Energy_loss_C)/Energy_loss_C;
disp Erros_real_C=  
disp(Chi_Error_C*100);

% �����͵���� Estimated energy loss on for simple sum
Chi_Error_ABC = abs(Chi_Energy_loss_A+Chi_Energy_loss_B+Chi_Energy_loss_C-Energy_loss)/Energy_loss;
disp Erros_real_ABC=  
disp(Chi_Error_ABC*100);

%% ��ͼ PLOTS
% A�� Phase A
figure;
[hist_a_y,hist_a_x]=hist(three_phase(:,1),50); %��Ϊ50������ͳ��
hist_a_y=hist_a_y/length(three_phase(:,1))/mean(diff(hist_a_x));   %��������ܶ� ��Ƶ�����������������������
bar(hist_a_x,hist_a_y,1);    
hold on 
gmPDF_A = pdf(gmm_fit_Pa,hist_a_x'); 
plot(hist_a_x,gmPDF_A,'linewidth',3);
legend('Measured data','GMM')
hold off
set(gca,'fontsize',14,'fontname','Times','box','off');
xlabel('Line currents (A)','Fontname', 'Times New Roman','FontSize',14);
ylabel('PDF','Fontname', 'Times New Roman','FontSize',14);
title('{\bf Phase A}')
% B�� Phase B
figure;
[hist_b_y,hist_b_x]=hist(three_phase(:,2),50); %��Ϊ50������ͳ��
hist_b_y=hist_b_y/length(three_phase(:,2))/mean(diff(hist_b_x));   %��������ܶ� ��Ƶ�����������������������
bar(hist_b_x,hist_b_y,1);    
hold on 
gmPDF_B = pdf(gmm_fit_Pb,hist_b_x'); 
plot(hist_b_x,gmPDF_B,'linewidth',3);
legend('Measured data','GMM')
hold off
set(gca,'fontsize',14,'fontname','Times','box','off');
xlabel('Line currents (A)','Fontname', 'Times New Roman','FontSize',14);
ylabel('PDF','Fontname', 'Times New Roman','FontSize',14);
title('{\bf Phase B}')

% C�� Phase C
figure;
[hist_c_y,hist_c_x]=hist(three_phase(:,3),50); %��Ϊ50������ͳ��
hist_c_y=hist_c_y/length(three_phase(:,3))/mean(diff(hist_c_x));   %��������ܶ� ��Ƶ�����������������������
bar(hist_c_x,hist_c_y,1);    
hold on 
gmPDF_C = pdf(gmm_fit_Pc,hist_c_x'); 
plot(hist_c_x,gmPDF_C,'linewidth',3);
legend('Measured data','GMM')
hold off
set(gca,'fontsize',14,'fontname','Times','box','off');
xlabel('Line currents (A)','Fontname', 'Times New Roman','FontSize',14);
ylabel('PDF','Fontname', 'Times New Roman','FontSize',14);
title('{\bf Phase C}')
%% Moment-matching use Gamma distribution
figure;

xrange = 1:12000;
GammaPhaseA = gampdf(xrange,1.9744, 492.136);
GammaPhaseB = gampdf(xrange,2.63188, 1400.81);
GammaPhaseC = gampdf(xrange,4.79132, 836.764);

plot(xrange,GammaPhaseA,xrange,GammaPhaseB,'--',xrange,GammaPhaseC,':','linewidth',3);
legend({'Phase A','Phase B','Phase C'})
legend('boxoff')
set(gca,'fontsize',14,'fontname','Times','box','off');
xlabel('Energy Loss (kWh)','Fontname', 'Times New Roman','FontSize',14);
ylabel('Probability density','Fontname', 'Times New Roman','FontSize',14);

%% Fitting true energy loss distribution
I2RPhaseA = Ia_square*ra*del_t*14400/(3600*1000);
I2RPhaseB = Ib_square*rb*del_t*14400/(3600*1000);
I2RPhaseC = Ic_square*rc*del_t*14400/(3600*1000);

nbins = 100;

[fPhaseA,xlableA] = ecdf(I2RPhaseA);
[fPhaseB,xlableB] = ecdf(I2RPhaseB);
[fPhasec,xlableC] = ecdf(I2RPhaseC);

plot(xlableA,fPhaseA,xlableB,fPhaseB,'--',xlableC,fPhasec,':','linewidth',3);
legend({'Phase A','Phase B','Phase C'})
legend('boxoff')
set(gca,'fontsize',14,'fontname','Times','box','off');
xlabel('Energy Loss (kWh)','Fontname', 'Times New Roman','FontSize',14);
ylabel('Probability density','Fontname', 'Times New Roman','FontSize',14);
