


%% Experiment 1 General Experiments 
clc
clearvars 
clear all
% path='OrderSelection_4th/';

%% Just change this part


name='forth_unkown';

AICc_Arr=[0;0;0;0;0]; %This is the value holder that help us decide the order 
AIC_Arr=[0;0;0;0;0]; %This is the value holder that help us decide the order 
SBC_Arr=[0;0;0;0;0]; %This is the value holder that help us decide the order 

Order=[2;3;4;5;6];




%% 4th order
x0 = [0;0;0; 1];

A=[1;5;5;0];
A_buff=[A;1];

%low noise 
% Noise_var=0.5;
% fs = 500 ; % sampling rate  120000

% 
% large noise 
Noise_var=1.5;
fs = 5000 ; % sampling rate  120000

% % large noise 
Noise_var=2;
fs =10000 ; % sampling rate  120000


% % Medium noise 
Noise_var=1;
fs = 5000 ; % sampling rate  120000


% path='OrderSelection_4th_low_many_run_write/';
% path='OrderSelection_4th_High_many_run_write/';
path='OrderSelection_4th_VERY_High_many_run_write/';
% path='OrderSelection_4th_medium_many_run_write/';




if ~exist(path, 'dir')
     mkdir(path)
end
write = 1; % to save the results and graphs







a = 0;
b=6;
t_array = a:(1/fs):b;
save('A_matrix.mat','A'); %wirte the estimate to the matrix to be used by ode
% We need to write A to the Matriz A %general case 
[t_array, x_array] = ode45(@utilities.sys_new_4th_order_LTI_editable,t_array,x0);
x_array = x_array';

y = [1 0 0 0]*x_array; %



y = y';
Y_deri_ture = zeros(length(y),4);
Y_deri_ture(:,1)=y;
Y_deri_ture(:,2)=([0 1 0 0]*x_array)'; %first derivative
Y_deri_ture(:,3)=([0 0 1 0]*x_array)'; %2nd
Y_deri_ture(:,4)=([0 0 0 1]*x_array)'; %3rd


y_est=zeros(length(y),1);
y_est2=zeros(length(y),1);
Num = 0; 

%% Generate and add noise 
%addition of noise%
y_noisy = utilities.White_noise_adder(y,0,Noise_var);
SNR = snr (y,y_noisy-y) ;

y_n = y_noisy;
y_old=y_n;

% %% smoothed 
window_size=fs/2;
y_smooth= smoothdata(y_n,'loess',window_size);
y_n=y_smooth;


figure
hold on 
plot(t_array,y)
plot(t_array,y_smooth)

grid on
legend('True Y','Smoothed Y','Location','northwest')
xlabel('Time t') 
ylabel('Value') 
F = getframe(gcf);
if write==1
    imwrite(F.cdata, fullfile(path,strcat(name,'true_Y_smoothed.png')),'png')
end 


immse( y_smooth , y)



figure
hold on 
plot(t_array(1:20:end),y_old(1:20:end),'o','MarkerSize',3)
plot(t_array,y)

grid on
legend('Noisy Y','True Y','Location','northwest')
xlabel('Time t') 
ylabel('Value') 
F = getframe(gcf);
if write==1
    imwrite(F.cdata, fullfile(path,strcat(name,'Noisy_True_Y.png')),'png')
end 

%% %%%%%%%%%%%%


%%%%
tic
%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% % second order 
n=2;
G_mat_regression_parallel=zeros(length(y_n),3); %assume 3th order 

% #for parallrl 
G_mat_regression_0=zeros(length(y_n),1);
G_mat_regression_1=zeros(length(y_n),1);
G_mat_regression_2=zeros(length(y_n),1);


parfor i = 1:1:length(t_array)
    if(mod(i,100)==0)
     fprintf('%d\n',i);
    end
    t_i = t_array(i);
    alfa_function=((t_i-a).^n +(b-t_i).^n);
    G_mat_regression_0(i)=(integral(@(t,v)utilities.K_fy_i(t,t_array,y_n,n,0,t_i,a),a,t_i,'RelTol',1e-6,'AbsTol',1e-10) + integral(@(t,v)utilities.K_By_i(t,t_array,y_n,n,0,t_i,b),t_i,b,'RelTol',1e-6,'AbsTol',1e-10))/alfa_function;
    G_mat_regression_1(i)=(integral(@(t,v)utilities.K_fy_i(t,t_array,y_n,n,1,t_i,a),a,t_i,'RelTol',1e-6,'AbsTol',1e-10) + integral(@(t,v)utilities.K_By_i(t,t_array,y_n,n,1,t_i,b),t_i,b,'RelTol',1e-6,'AbsTol',1e-10))/alfa_function;
    G_mat_regression_2(i)=(integral(@(t,v)utilities.K_fy_i(t,t_array,y_n,n,2,t_i,a),a,t_i,'RelTol',1e-6,'AbsTol',1e-10) + integral(@(t,v)utilities.K_By_i(t,t_array,y_n,n,2,t_i,b),t_i,b,'RelTol',1e-6,'AbsTol',1e-10))/alfa_function;

end




G_mat_regression_parallel = [G_mat_regression_0 G_mat_regression_1];
G= [G_mat_regression_parallel G_mat_regression_2];

Y_dummy=y_old; 
[Beta,FitInfo] = lasso(G,Y_dummy,'Alpha',0.9,'CV',12,'Options',statset('UseParallel',true));  % 1 is ridgit vs 0.01 is lasso
lassoPlot(Beta,FitInfo,'PlotType','CV');
legend('show') % Show legend
if write==1
    F = getframe(gcf);
    imwrite(F.cdata, fullfile(path,strcat(name,'lasso_lambda_2nd.png')),'png')
end 
index1 = FitInfo.Index1SE;
coef = Beta(:,index1)
index2 = FitInfo.IndexMinMSE;
coef = Beta(:,index2)
A_hat=coef;
A_estimate=A_hat;
disp(A_estimate);
%fundimental solutions:
x0_fund = [[1; 0] [0; 1]];
t_array = a:(1/fs):b;
y_fund=zeros(length(y_n),2);

A=A_estimate;
save('A_matrix.mat','A'); 
for i= 1:2
    [t_array, x_array] = ode45(@utilities.sys_new_2th_order_LTI_editable,t_array,x0_fund(i,:));
    x_array = x_array';
    y_fund(:,i) = ([1 0]*x_array)'; %ISSUE Is it X0 or X1 
   
end
Q=utilities.gram_schmidt(y_fund,n );
y_m =y_n; % the reason for this is to be similar to the notes
% Now we need to find the values of Beta 
buffer=y_m .* Q;
%integral
Beta= trapz(buffer);
y_rep=zeros(length(y_m),1);
%y_reproduced
%tic
parfor i = 1:1:length(t_array) %PARALLEL 

    y_rep(i) = Q(i,:) * Beta';
end
% %toc
SBC_Arr(1)=utilities.Aic(y_old,y_rep,n,1);
AIC_Arr(1)=utilities.Aic(y_old,y_rep,n,2);
AICc_Arr(1)=utilities.Aic(y_old,y_rep,n,0);

figure 
hold on 
plot(t_array,y)
plot(t_array,y_rep)
grid on
legend('True Y','Y estimate','Location','northwest')
xlabel('Time t') 
ylabel('Value') 
title('2 Order Estimation')
F = getframe(gcf);
if write==1
    imwrite(F.cdata, fullfile(path,strcat(name,'Estimate_True_Y_2nd.png')),'png')
end 

err = immse( y_rep , y);

if write==1
    dlmwrite( fullfile(path,strcat(name,'MSE_2_error.txt')),err,'delimiter','\t','precision',8)
    dlmwrite( fullfile(path,strcat(name,'A_2_estimate.txt')),A_estimate,'delimiter','\t','precision',8)
%     dlmwrite( fullfile(path,strcat(name,'SNR.txt')),SNR,'delimiter','\t','precision',8)

end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% % third order 
n=3;
% G_mat_regression_parallel=zeros(length(y_n),length(A_buff));
G_mat_regression_parallel=zeros(length(y_n),4); %assume 3th order 

% #for parallrl 
G_mat_regression_0=zeros(length(y_n),1);
G_mat_regression_1=zeros(length(y_n),1);
G_mat_regression_2=zeros(length(y_n),1);
G_mat_regression_3=zeros(length(y_n),1);

parfor i = 1:1:length(t_array)
    if(mod(i,100)==0)
     fprintf('%d\n',i);
    end
    t_i = t_array(i);
    alfa_function=((t_i-a).^n +(b-t_i).^n);
    G_mat_regression_0(i)=(integral(@(t,v)utilities.K_fy_i(t,t_array,y_n,n,0,t_i,a),a,t_i,'RelTol',1e-6,'AbsTol',1e-10) + integral(@(t,v)utilities.K_By_i(t,t_array,y_n,n,0,t_i,b),t_i,b,'RelTol',1e-6,'AbsTol',1e-10))/alfa_function;
    G_mat_regression_1(i)=(integral(@(t,v)utilities.K_fy_i(t,t_array,y_n,n,1,t_i,a),a,t_i,'RelTol',1e-6,'AbsTol',1e-10) + integral(@(t,v)utilities.K_By_i(t,t_array,y_n,n,1,t_i,b),t_i,b,'RelTol',1e-6,'AbsTol',1e-10))/alfa_function;
    G_mat_regression_2(i)=(integral(@(t,v)utilities.K_fy_i(t,t_array,y_n,n,2,t_i,a),a,t_i,'RelTol',1e-6,'AbsTol',1e-10) + integral(@(t,v)utilities.K_By_i(t,t_array,y_n,n,2,t_i,b),t_i,b,'RelTol',1e-6,'AbsTol',1e-10))/alfa_function;
    G_mat_regression_3(i)=(integral(@(t,v)utilities.K_fy_i(t,t_array,y_n,n,3,t_i,a),a,t_i,'RelTol',1e-6,'AbsTol',1e-10) + integral(@(t,v)utilities.K_By_i(t,t_array,y_n,n,3,t_i,b),t_i,b,'RelTol',1e-6,'AbsTol',1e-10))/alfa_function;

end

G_mat_regression_parallel = [G_mat_regression_0 G_mat_regression_1 G_mat_regression_2 ];
G= [G_mat_regression_parallel G_mat_regression_3];

Y_dummy=y_old; 
[Beta,FitInfo] = lasso(G,Y_dummy,'Alpha',0.9,'CV',12,'Options',statset('UseParallel',true));  % 1 is ridgit vs 0.01 is lasso
lassoPlot(Beta,FitInfo,'PlotType','CV');
legend('show') % Show legend
if write==1
    F = getframe(gcf);
    imwrite(F.cdata, fullfile(path,strcat(name,'lasso_lambda_3rd.png')),'png')
end 

index1 = FitInfo.Index1SE;
coef = Beta(:,index1)
index2 = FitInfo.IndexMinMSE;
coef = Beta(:,index2)
A_hat=coef;

A_estimate=A_hat;
% writematrix(err,'MSE_error.txt','Delimiter','tab') % worksis in R2019
disp(A_estimate);
%fundimental solutions:
x0_fund = [[1; 0; 0] [0; 1; 0] [0; 0; 1]];
t_array = a:(1/fs):b;
y_fund=zeros(length(y_n),3);

A=A_estimate;
save('A_matrix.mat','A'); 
for i= 1:3
    [t_array, x_array] = ode45(@utilities.sys_new_3th_order_LTI_editable,t_array,x0_fund(i,:));
    x_array = x_array';
    y_fund(:,i) = ([1 0 0]*x_array)'; %ISSUE Is it X0 or X1 
   
end
Q=utilities.gram_schmidt(y_fund,n );
y_m =y_n; % the reason for this is to be similar to the notes
% Now we need to find the values of Beta 
buffer=y_m .* Q;
%integral
Beta= trapz(buffer);
y_rep=zeros(length(y_m),1);
%y_reproduced
%tic
parfor i = 1:1:length(t_array) %PARALLEL 

    y_rep(i) = Q(i,:) * Beta';
end
% %toc
SBC_Arr(2)=utilities.Aic(y_old,y_rep,n,1);
AIC_Arr(2)=utilities.Aic(y_old,y_rep,n,2);
AICc_Arr(2)=utilities.Aic(y_old,y_rep,n,0);







figure 
hold on 
plot(t_array,y)
plot(t_array,y_rep)

grid on
legend('True Y','Y estimate','Location','northwest')
xlabel('Time t') 
ylabel('Value') 
title('Third Order Estimation')
F = getframe(gcf);
if write==1
    imwrite(F.cdata, fullfile(path,strcat(name,'Estimate_True_Y_3rd.png')),'png')
end 

err = immse( y_rep , y);

if write==1
    dlmwrite( fullfile(path,strcat(name,'MSE_3_error.txt')),err,'delimiter','\t','precision',8)
    dlmwrite( fullfile(path,strcat(name,'A_3_estimate.txt')),A_estimate,'delimiter','\t','precision',8)

end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 5th order 
n=5; 
G_mat_regression_parallel=zeros(length(y_n),6); %assume 2th order 
G_mat_regression_0=zeros(length(y_n),1);
G_mat_regression_1=zeros(length(y_n),1);
G_mat_regression_2=zeros(length(y_n),1);
G_mat_regression_3=zeros(length(y_n),1);
G_mat_regression_4=zeros(length(y_n),1);
G_mat_regression_5=zeros(length(y_n),1);

tic %timer 
parfor i = 1:1:length(t_array)
    if(mod(i,100)==0)
     fprintf('%d\n',i);
    end
    t_i = t_array(i);
    alfa_function=((t_i-a).^n +(b-t_i).^n);
    
    G_mat_regression_0(i)=(integral(@(t,v)utilities.K_fy_i(t,t_array,y_n,n,0,t_i,a),a,t_i,'RelTol',1e-6,'AbsTol',1e-10) + integral(@(t,v)utilities.K_By_i(t,t_array,y_n,n,0,t_i,b),t_i,b,'RelTol',1e-6,'AbsTol',1e-10))/alfa_function;
    G_mat_regression_1(i)=(integral(@(t,v)utilities.K_fy_i(t,t_array,y_n,n,1,t_i,a),a,t_i,'RelTol',1e-6,'AbsTol',1e-10) + integral(@(t,v)utilities.K_By_i(t,t_array,y_n,n,1,t_i,b),t_i,b,'RelTol',1e-6,'AbsTol',1e-10))/alfa_function;
    G_mat_regression_2(i)=(integral(@(t,v)utilities.K_fy_i(t,t_array,y_n,n,2,t_i,a),a,t_i,'RelTol',1e-6,'AbsTol',1e-10) + integral(@(t,v)utilities.K_By_i(t,t_array,y_n,n,2,t_i,b),t_i,b,'RelTol',1e-6,'AbsTol',1e-10))/alfa_function;
    G_mat_regression_3(i)=(integral(@(t,v)utilities.K_fy_i(t,t_array,y_n,n,3,t_i,a),a,t_i,'RelTol',1e-6,'AbsTol',1e-10) + integral(@(t,v)utilities.K_By_i(t,t_array,y_n,n,3,t_i,b),t_i,b,'RelTol',1e-6,'AbsTol',1e-10))/alfa_function;
    G_mat_regression_4(i)=(integral(@(t,v)utilities.K_fy_i(t,t_array,y_n,n,4,t_i,a),a,t_i,'RelTol',1e-6,'AbsTol',1e-10) + integral(@(t,v)utilities.K_By_i(t,t_array,y_n,n,4,t_i,b),t_i,b,'RelTol',1e-6,'AbsTol',1e-10))/alfa_function;
    G_mat_regression_5(i)=(integral(@(t,v)utilities.K_fy_i(t,t_array,y_n,n,5,t_i,a),a,t_i,'RelTol',1e-6,'AbsTol',1e-10) + integral(@(t,v)utilities.K_By_i(t,t_array,y_n,n,5,t_i,b),t_i,b,'RelTol',1e-6,'AbsTol',1e-10))/alfa_function;

 end
toc
G_mat_regression_parallel = [G_mat_regression_0 G_mat_regression_1 G_mat_regression_2 G_mat_regression_3 G_mat_regression_4];
G= [G_mat_regression_parallel G_mat_regression_5];


Y_dummy=y_old;
[Beta,FitInfo] = lasso(G,Y_dummy,'Alpha',0.9,'CV',12,'Options',statset('UseParallel',true));  % 1 is ridgit vs 0.01 is lasso
lassoPlot(Beta,FitInfo,'PlotType','CV');
legend('show') 
if write==1
    F = getframe(gcf);
    imwrite(F.cdata, fullfile(path,strcat(name,'lasso_lambda_5th.png')),'png')
end 

index1 = FitInfo.Index1SE;
coef = Beta(:,index1)
index2 = FitInfo.IndexMinMSE;
coef = Beta(:,index2)
A_hat=coef;
A_estimate=A_hat;
disp(A_estimate);
%fundimental solutions:
x0_fund = [[1; 0; 0; 0;0]          [0; 1; 0; 0;0]          [0; 0; 1; 0;0]     [0; 0; 0; 1;0]      [0; 0; 0; 0;1]];
t_array = a:(1/fs):b;
y_fund=zeros(length(y_n),5);

A=A_estimate;
save('A_matrix.mat','A'); 
for i= 1:5
    [t_array, x_array] = ode45(@utilities.sys_new_5th_order_LTI_editable,t_array,x0_fund(i,:));
    x_array = x_array';
    y_fund(:,i) = ([1 0 0 0 0]*x_array)'; %ISSUE Is it X0 or X1 
   
end
Q=utilities.gram_schmidt(y_fund,n );
y_m =y_n; % the reason for this is to be similar to the notes
% Now we need to find the values of Beta 
buffer=y_m .* Q;
%integral
Beta= trapz(buffer);
y_rep=zeros(length(y_m),1);
%y_reproduced
%tic
parfor i = 1:1:length(t_array) %PARALLEL 
    y_rep(i) = Q(i,:) * Beta';
end
% %toc
SBC_Arr(4)=utilities.Aic(y_old,y_rep,n,1);
AIC_Arr(4)=utilities.Aic(y_old,y_rep,n,2);
AICc_Arr(4)=utilities.Aic(y_old,y_rep,n,0);


figure 
hold on 
plot(t_array,y)
plot(t_array,y_rep)

grid on
legend('True Y','Y estimate','Location','northwest')
xlabel('Time t') 
ylabel('Value') 
title('Fifth Order Estimation')
F = getframe(gcf);
if write==1
    imwrite(F.cdata, fullfile(path,strcat(name,'Estimate_True_Y_5th.png')),'png')
end 

err = immse( y_rep , y);


if write==1
    dlmwrite( fullfile(path,strcat(name,'MSE_5_error.txt')),err,'delimiter','\t','precision',8)
    dlmwrite( fullfile(path,strcat(name,'A_5_estimate.txt')),A_estimate,'delimiter','\t','precision',8)
end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 6th order 
n=6; 
G_mat_regression_parallel=zeros(length(y_n),7); %assume 2th order 
G_mat_regression_0=zeros(length(y_n),1);
G_mat_regression_1=zeros(length(y_n),1);
G_mat_regression_2=zeros(length(y_n),1);
G_mat_regression_3=zeros(length(y_n),1);
G_mat_regression_4=zeros(length(y_n),1);
G_mat_regression_5=zeros(length(y_n),1);
G_mat_regression_6=zeros(length(y_n),1);


tic %timer 
parfor i = 1:1:length(t_array)
    if(mod(i,100)==0)
     fprintf('%d\n',i);
    end
    t_i = t_array(i);
    alfa_function=((t_i-a).^n +(b-t_i).^n);
    
    G_mat_regression_0(i)=(integral(@(t,v)utilities.K_fy_i(t,t_array,y_n,n,0,t_i,a),a,t_i,'RelTol',1e-6,'AbsTol',1e-10) + integral(@(t,v)utilities.K_By_i(t,t_array,y_n,n,0,t_i,b),t_i,b,'RelTol',1e-6,'AbsTol',1e-10))/alfa_function;
    G_mat_regression_1(i)=(integral(@(t,v)utilities.K_fy_i(t,t_array,y_n,n,1,t_i,a),a,t_i,'RelTol',1e-6,'AbsTol',1e-10) + integral(@(t,v)utilities.K_By_i(t,t_array,y_n,n,1,t_i,b),t_i,b,'RelTol',1e-6,'AbsTol',1e-10))/alfa_function;
    G_mat_regression_2(i)=(integral(@(t,v)utilities.K_fy_i(t,t_array,y_n,n,2,t_i,a),a,t_i,'RelTol',1e-6,'AbsTol',1e-10) + integral(@(t,v)utilities.K_By_i(t,t_array,y_n,n,2,t_i,b),t_i,b,'RelTol',1e-6,'AbsTol',1e-10))/alfa_function;
    G_mat_regression_3(i)=(integral(@(t,v)utilities.K_fy_i(t,t_array,y_n,n,3,t_i,a),a,t_i,'RelTol',1e-6,'AbsTol',1e-10) + integral(@(t,v)utilities.K_By_i(t,t_array,y_n,n,3,t_i,b),t_i,b,'RelTol',1e-6,'AbsTol',1e-10))/alfa_function;
    G_mat_regression_4(i)=(integral(@(t,v)utilities.K_fy_i(t,t_array,y_n,n,4,t_i,a),a,t_i,'RelTol',1e-6,'AbsTol',1e-10) + integral(@(t,v)utilities.K_By_i(t,t_array,y_n,n,4,t_i,b),t_i,b,'RelTol',1e-6,'AbsTol',1e-10))/alfa_function;
    G_mat_regression_5(i)=(integral(@(t,v)utilities.K_fy_i(t,t_array,y_n,n,5,t_i,a),a,t_i,'RelTol',1e-6,'AbsTol',1e-10) + integral(@(t,v)utilities.K_By_i(t,t_array,y_n,n,5,t_i,b),t_i,b,'RelTol',1e-6,'AbsTol',1e-10))/alfa_function;
    G_mat_regression_6(i)=(integral(@(t,v)utilities.K_fy_i(t,t_array,y_n,n,6,t_i,a),a,t_i,'RelTol',1e-6,'AbsTol',1e-10) + integral(@(t,v)utilities.K_By_i(t,t_array,y_n,n,6,t_i,b),t_i,b,'RelTol',1e-6,'AbsTol',1e-10))/alfa_function;

 end
toc
G_mat_regression_parallel = [G_mat_regression_0 G_mat_regression_1 G_mat_regression_2 G_mat_regression_3 G_mat_regression_4 G_mat_regression_5];
G= [G_mat_regression_parallel G_mat_regression_6];


Y_dummy=y_old;
[Beta,FitInfo] = lasso(G,Y_dummy,'Alpha',0.9,'CV',12,'Options',statset('UseParallel',true));  % 1 is ridgit vs 0.01 is lasso
lassoPlot(Beta,FitInfo,'PlotType','CV');
legend('show') 
if write==1
    F = getframe(gcf);
    imwrite(F.cdata, fullfile(path,strcat(name,'lasso_lambda_6th.png')),'png')
end 

index1 = FitInfo.Index1SE;
coef = Beta(:,index1)
index2 = FitInfo.IndexMinMSE;
coef = Beta(:,index2)
A_hat=coef;
A_estimate=A_hat;
disp(A_estimate);
%fundimental solutions:
x0_fund = [[1; 0; 0; 0;0;0]          [0; 1; 0; 0;0;0]          [0; 0; 1; 0;0;0]     [0; 0; 0; 1;0;0]      [0; 0; 0; 0;1;0] [0; 0; 0; 0;0;1]];
t_array = a:(1/fs):b;
y_fund=zeros(length(y_n),6);

A=A_estimate;
save('A_matrix.mat','A'); 
for i= 1:6
    [t_array, x_array] = ode45(@utilities.sys_new_6th_order_LTI_editable,t_array,x0_fund(i,:));
    x_array = x_array';
    y_fund(:,i) = ([1 0 0 0 0 0]*x_array)'; %ISSUE Is it X0 or X1 
   
end
Q=utilities.gram_schmidt(y_fund,n );
y_m =y_n; % the reason for this is to be similar to the notes
% Now we need to find the values of Beta 
buffer=y_m .* Q;
%integral
Beta= trapz(buffer);
y_rep=zeros(length(y_m),1);
%y_reproduced
%tic
parfor i = 1:1:length(t_array) %PARALLEL 
    y_rep(i) = Q(i,:) * Beta';
end
% %toc
SBC_Arr(5)=utilities.Aic(y_old,y_rep,n,1);
AIC_Arr(5)=utilities.Aic(y_old,y_rep,n,2);
AICc_Arr(5)=utilities.Aic(y_old,y_rep,n,0);


figure 
hold on 
plot(t_array,y)
plot(t_array,y_rep)

grid on
legend('True Y','Y estimate','Location','northwest')
xlabel('Time t') 
ylabel('Value') 
title('Six Order Estimation')
F = getframe(gcf);
if write==1
    imwrite(F.cdata, fullfile(path,strcat(name,'Estimate_True_Y_6th.png')),'png')
end 

err = immse( y_rep , y);


if write==1
    dlmwrite( fullfile(path,strcat(name,'MSE_6_error.txt')),err,'delimiter','\t','precision',8)
    dlmwrite( fullfile(path,strcat(name,'A_6_estimate.txt')),A_estimate,'delimiter','\t','precision',8)
end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% % 4th order 
n=4;% the order of our system  

% G_mat_regression_parallel=zeros(length(y_n),length(A_buff));
G_mat_regression_parallel=zeros(length(y_n),5); %assume 4th order 

% #for parallrl 
G_mat_regression_0=zeros(length(y_n),1);
G_mat_regression_1=zeros(length(y_n),1);
G_mat_regression_2=zeros(length(y_n),1);
G_mat_regression_3=zeros(length(y_n),1);
G_mat_regression_4=zeros(length(y_n),1);
%first we need to find G
% we also update y_n
% The first run is slower than the second run, because the parallel pool has to be started

t_array_copy_1=t_array; % to avoid overhead comunicatiion in matlab
max_count=10000;
parfor i = 1:1:length(t_array)
    if(mod(i,100)==0)
     fprintf('%d\n',i);
    end
    t_i = t_array(i);
    alfa_function=((t_i-a).^n +(b-t_i).^n);
    G_mat_regression_0(i)=(integral(@(t,v)utilities.K_fy_i(t,t_array,y_n,n,0,t_i,a),a,t_i,'RelTol',1e-6,'AbsTol',1e-10) + integral(@(t,v)utilities.K_By_i(t,t_array,y_n,n,0,t_i,b),t_i,b,'RelTol',1e-6,'AbsTol',1e-10))/alfa_function;
    G_mat_regression_1(i)=(integral(@(t,v)utilities.K_fy_i(t,t_array,y_n,n,1,t_i,a),a,t_i,'RelTol',1e-6,'AbsTol',1e-10) + integral(@(t,v)utilities.K_By_i(t,t_array,y_n,n,1,t_i,b),t_i,b,'RelTol',1e-6,'AbsTol',1e-10))/alfa_function;
    G_mat_regression_2(i)=(integral(@(t,v)utilities.K_fy_i(t,t_array,y_n,n,2,t_i,a),a,t_i,'RelTol',1e-6,'AbsTol',1e-10) + integral(@(t,v)utilities.K_By_i(t,t_array,y_n,n,2,t_i,b),t_i,b,'RelTol',1e-6,'AbsTol',1e-10))/alfa_function;
    G_mat_regression_3(i)=(integral(@(t,v)utilities.K_fy_i(t,t_array,y_n,n,3,t_i,a),a,t_i,'RelTol',1e-6,'AbsTol',1e-10) + integral(@(t,v)utilities.K_By_i(t,t_array,y_n,n,3,t_i,b),t_i,b,'RelTol',1e-6,'AbsTol',1e-10))/alfa_function;
    G_mat_regression_4(i)=(integral(@(t,v)utilities.K_fy_i(t,t_array,y_n,n,4,t_i,a),a,t_i,'RelTol',1e-6,'AbsTol',1e-10) + integral(@(t,v)utilities.K_By_i(t,t_array,y_n,n,4,t_i,b),t_i,b,'RelTol',1e-6,'AbsTol',1e-10))/alfa_function;

end



figure 
hold on 
plot(t_array,G_mat_regression_0)
plot(t_array,G_mat_regression_1)
plot(t_array,G_mat_regression_2)
plot(t_array,G_mat_regression_3)
plot(t_array,G_mat_regression_4)


forth_G_mat_regression_0=G_mat_regression_0;
forth_G_mat_regression_1=G_mat_regression_1;
forth_G_mat_regression_2=G_mat_regression_2;
forth_G_mat_regression_3=G_mat_regression_3;
forth_G_mat_regression_4=G_mat_regression_4;


G_mat_regression_parallel = [G_mat_regression_0 G_mat_regression_1 G_mat_regression_2 G_mat_regression_3];
G= [G_mat_regression_parallel G_mat_regression_4];


Y_dummy=y_old;
[Beta,FitInfo] = lasso(G,Y_dummy,'Alpha',0.9,'CV',12,'Options',statset('UseParallel',true));  % 1 is ridgit vs 0.01 is lasso
lassoPlot(Beta,FitInfo,'PlotType','CV');
legend('show') % Show legend
if write==1
    F = getframe(gcf);
    imwrite(F.cdata, fullfile(path,strcat(name,'lasso_lambda_4th.png')),'png')
end 


index1 = FitInfo.Index1SE;
coef = Beta(:,index1)
index2 = FitInfo.IndexMinMSE;
coef = Beta(:,index2)

A_hat=coef;

% A_hat=regress(Y_dummy,G_mat_regression_parallel);


A_estimate=A_hat;
% writematrix(err,'MSE_error.txt','Delimiter','tab') % worksis in R2019
disp(A_estimate);
%fundimental solutions:
x0_fund = [[1; 0; 0; 0]           [0; 1; 0; 0]          [0; 0; 1; 0]           [0; 0; 0; 1]];
t_array = a:(1/fs):b;
y_fund=zeros(length(y_n),4);

A=A_estimate;
save('A_matrix.mat','A'); 
for i= 1:4
    [t_array, x_array] = ode45(@utilities.sys_new_4th_order_LTI_editable,t_array,x0_fund(i,:));
    x_array = x_array';
    y_fund(:,i) = ([1 0 0 0]*x_array)'; %ISSUE Is it X0 or X1 
   
end
Q=utilities.gram_schmidt(y_fund,n );
y_m =y_n; % the reason for this is to be similar to the notes
% Now we need to find the values of Beta 
buffer=y_m .* Q;
%integral
Beta= trapz(buffer);
y_rep=zeros(length(y_m),1);
%y_reproduced
%tic
parfor i = 1:1:length(t_array) %PARALLEL 
    y_rep(i) = Q(i,:) * Beta';
end
% %toc
SBC_Arr(3)=utilities.Aic(y_old,y_rep,n,1);
AIC_Arr(3)=utilities.Aic(y_old,y_rep,n,2);
AICc_Arr(3)=utilities.Aic(y_old,y_rep,n,0);



figure 
hold on 
plot(t_array,y)
plot(t_array,y_rep)

grid on
legend('True Y','Y estimate','Location','northwest')
xlabel('Time t') 
ylabel('Value') 
title('Fourth  Order Estimation')
F = getframe(gcf);
if write==1
    imwrite(F.cdata, fullfile(path,strcat(name,'Estimate_True_Y_4th.png')),'png')
end 



err = immse( y_rep , y);

if write==1
    dlmwrite( fullfile(path,strcat(name,'MSE_4_error.txt')),err,'delimiter','\t','precision',8)
    dlmwrite( fullfile(path,strcat(name,'A_4_estimate.txt')),A_estimate,'delimiter','\t','precision',8)
end 










%% % decide the order 

%% % Compare AICc then decide the one we want to go 
figure
hold on
plot(Order,AICc_Arr)
legend('AICc')
grid on
xlabel('The Estimate Model Order ') 
ylabel('Value AICc') 
F = getframe(gcf);
if write==1
    imwrite(F.cdata, fullfile(path,strcat(name,'AICc.png')),'png')
end 



figure
hold on
plot(Order,SBC_Arr)
legend('SBC')
grid on
xlabel('The Estimate Model Order ') 
ylabel('Value SBC') 
F = getframe(gcf);
if write==1
    imwrite(F.cdata, fullfile(path,strcat(name,'SBC.png')),'png')
end

figure
hold on
plot(Order,AIC_Arr)
legend('AIC')
grid on
xlabel('The Estimate Model Order ') 
ylabel('Value AIC') 
F = getframe(gcf);
if write==1
    imwrite(F.cdata, fullfile(path,strcat(name,'AIC.png')),'png')
end 




if write==1
    dlmwrite( fullfile(path,strcat(name,'AICc_.txt')),AICc_Arr,'delimiter','\t','precision',8)
    dlmwrite( fullfile(path,strcat(name,'AIC_.txt')),AIC_Arr,'delimiter','\t','precision',8)
    dlmwrite( fullfile(path,strcat(name,'SBC_.txt')),SBC_Arr,'delimiter','\t','precision',8)

end 



%% it is the 4th order 


%% Retrain and evaluate the selected order and derivatives


if write==1

    dlmwrite( fullfile(path,strcat(name,'MSE_error.txt')),err,'delimiter','\t','precision',8)
    dlmwrite( fullfile(path,strcat(name,'A_estimate.txt')),A_estimate,'delimiter','\t','precision',8)
    dlmwrite( fullfile(path,strcat(name,'SNR.txt')),SNR,'delimiter','\t','precision',8)
    dlmwrite( fullfile(path,strcat(name,'coef_fromPen.txt')),coef,'delimiter','\t','precision',8)
    dlmwrite( fullfile(path,strcat(name,'_fs.txt')),fs,'delimiter','\t','precision',8)

end 

%% Derivatives 
Y_deri = zeros(length(y),4);
Y_deri(:,1)=y_rep;
for i=1:1:3
    %%change the y_rep
    Y_deri(:,i+1)=utilities.deriv_helper(t_array,y_rep,A_estimate,4,a,b, i ,Y_deri);
%     Y_deri(:,i+1)=utilities.deriv_helper(t_array,y_n,A_estimate,4,a,b, i ,Y_deri);

    %plot the  true Y and estimated
    figure 
    hold on 
    plot(t_array,Y_deri(:,i+1))
    plot(t_array,Y_deri_ture(:,i+1))   
    grid on
    buff1 = sprintf('True Y ^{(%d)} ', i);
    buff2 = sprintf('Estimated Y ^{(%d)}', i);
    legend(buff1,buff2,'Location','northwest')
    xlabel('Time t') 
    ylabel('Value') 
    
    
    if write==1
        F = getframe(gcf);
        baseFileName = sprintf('deriv_Y_%d.png', i);
        imwrite(F.cdata, fullfile(path,strcat(name,baseFileName)),'png')
    end 
end

