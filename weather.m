opts = detectImportOptions('data/temperature.csv');
T = readtable('data/temperature.csv',opts);
y_noisy=T.SanFrancisco;
id2 = find(isnan(y_noisy));
y_noisy([id2]) = [];
y;
x=T.datetime;
x([id2]) = [];
t_array = 1:1:size(y_noisy,1);
% figure; hold on; grid on;
% plot(t_array,y_noisy)



y_n = y_noisy;
y_old=y_n;
% window_size=fs+100 ;
window_size=fs/2;
% window_size=fs/8 -100;
% y_smooth= smoothdata(y_n,'lowess',window_size);
y_smooth= smoothdata(y_n,'loess',window_size);



figure
hold on 
plot(t_array,y_noisy)
% plot(t_array,y)
plot(t_array,y_smooth)


%% 



% % % % % % % % % % % % % % % % % % % % % % 5555555555
grid on
legend('True Y','Smoothed Y','Location','northeast')
xlabel('Time t') 
ylabel('Value') 
% F = getframe(gcf);
% % if write==1
% %     imwrite(F.cdata, fullfile(path,strcat(name,'true_Y_smoothed.png')),'png')%eps
% % end 



y_n=y_smooth;
% y_n=y;
figure
hold on 
plot(t_array,y_n)


figure
hold on 
plot(t_array(1:1:end),y_old(1:1:end),'o','MarkerSize',3)

grid on
legend('Noisy Y','True Y','Location','northeast')
xlabel('Time t') 
ylabel('Value') 
% 
% if write==1
%     F = getframe(gcf);
%     imwrite(F.cdata, fullfile(path,strcat(name,'Noisy_True_Y.png')),'png')
% end 


%%%%%
tic
%%%%%



G_mat_regression_parallel=zeros(length(y_n),4);

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

% y_n=y_old;

max_count=10000;
tic %timer 
parfor i = 1:1:length(t_array)
    if(mod(i,100)==0)
     fprintf('%d\n',i);
    end
    t_i = t_array(i);
    alfa_function=((t_i-a).^n +(b-t_i).^n);
    G_mat_regression_0(i)=(integral(@(t,v)utilities.K_fy_i(t,t_array,y_n,n,0,t_i,a),a,t_i,'RelTol',1e-2,'AbsTol',1e-6) + integral(@(t,v)utilities.K_By_i(t,t_array,y_n,n,0,t_i,b),t_i,b,'RelTol',1e-2,'AbsTol',1e-6))/alfa_function;
    G_mat_regression_1(i)=(integral(@(t,v)utilities.K_fy_i(t,t_array,y_n,n,1,t_i,a),a,t_i,'RelTol',1e-2,'AbsTol',1e-6) + integral(@(t,v)utilities.K_By_i(t,t_array,y_n,n,1,t_i,b),t_i,b,'RelTol',1e-2,'AbsTol',1e-6))/alfa_function;
    G_mat_regression_2(i)=(integral(@(t,v)utilities.K_fy_i(t,t_array,y_n,n,2,t_i,a),a,t_i,'RelTol',1e-2,'AbsTol',1e-6) + integral(@(t,v)utilities.K_By_i(t,t_array,y_n,n,2,t_i,b),t_i,b,'RelTol',1e-2,'AbsTol',1e-6))/alfa_function;
    G_mat_regression_3(i)=(integral(@(t,v)utilities.K_fy_i(t,t_array,y_n,n,3,t_i,a),a,t_i,'RelTol',1e-2,'AbsTol',1e-6) + integral(@(t,v)utilities.K_By_i(t,t_array,y_n,n,3,t_i,b),t_i,b,'RelTol',1e-2,'AbsTol',1e-6))/alfa_function;
    G_mat_regression_4(i)=(integral(@(t,v)utilities.K_fy_i(t,t_array,y_n,n,4,t_i,a),a,t_i,'RelTol',1e-2,'AbsTol',1e-6) + integral(@(t,v)utilities.K_By_i(t,t_array,y_n,n,4,t_i,b),t_i,b,'RelTol',1e-2,'AbsTol',1e-6))/alfa_function;

end
toc
 




G_mat_regression_parallel = [G_mat_regression_0 G_mat_regression_1 G_mat_regression_2 G_mat_regression_3];

G= [G_mat_regression_parallel G_mat_regression_4];

%%


figure
hold on 
plot(t_array,y_n)
plot(t_array,G_mat_regression_4)



%% 
Y_dummy=y_old;
[Beta,FitInfo] = lasso(G,Y_dummy,'Alpha',0.9,'CV',10);  % 1 is ridgit vs 0.01 is lasso
lassoPlot(Beta,FitInfo,'PlotType','CV');
legend('show') % Show legend

index1 = FitInfo.Index1SE;
coef = Beta(:,index1)
index2 = FitInfo.IndexMinMSE;
coef = Beta(:,index2)
% 
idxLambda1SE=index2+floor(abs(index2-index1)/2);
coef = Beta(:,idxLambda1SE)
coef0 = FitInfo.Intercept(idxLambda1SE);
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
Y_dummy=y_old-G_mat_regression_4;


[Beta,FitInfo] = lasso(G_mat_regression_parallel,Y_dummy,'Alpha',0.9,'CV',12,'Options',statset('UseParallel',true));  % 1 is ridgit vs 0.01 is lasso
lassoPlot(Beta,FitInfo,'PlotType','CV');
legend('show') % Show legend
if write==1
    F = getframe(gcf);
    imwrite(F.cdata, fullfile(path,strcat(name,'lasso_lambda.png')),'png')
end 


index1 = FitInfo.Index1SE;
coef = Beta(:,index1)
index2 = FitInfo.IndexMinMSE;
coef = Beta(:,index2)

% idxLambda1SE=index2+floor(abs(index2-index1)/2);
% coef = Beta(:,idxLambda1SE)
% coef0 = FitInfo.Intercept(idxLambda1SE);
% 
% 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% now we are doing the fundimental solutions:
%%%%%%%%%%%%%%%%
% Here we will create fundamental solution usinging  different starting
% points  # true value of our system

x0_fund = [[1; 0; 0; 0]           [0; 1; 0; 0]          [0; 0; 1; 0]           [0; 0; 0; 1]];
t_array = a:(1/fs):b;
y_fund=zeros(length(y_n),4);
%general case 



%% What if I use Coef from lasso 
A=coef;
%% What if I use Coef from lasso 
A_estimate=coef;



save('A_matrix.mat','A'); %wirte the estimate to the matrix to be used by ode
% We need to write A to the Matriz A 
% options = odeset('A', A); 
for i= 1:4
    [t_array, x_array] = ode45(@utilities.sys_new_4th_order_LTI_editable,t_array,x0_fund(i,:));
    x_array = x_array';
    y_fund(:,i) = ([1 0 0 0]*x_array)'; %ISSUE Is it X0 or X1 
   
end



% # gram-schmidt
Q=utilities.gram_schmidt(y_fund,n );


y_m =y_n; % the reason for this is to be similar to the notes
y_m=y_old;  %Not sure about this 

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


figure 
hold on 
plot(t_array,y)
plot(t_array,y_rep)

grid on
legend('True Y','Y estimate','Location','northeast')
xlabel('Time t') 
ylabel('Value') 
F = getframe(gcf);
if write==1
    imwrite(F.cdata, fullfile(path,strcat(name,'Estimate_True_Y.png')),'png')
end 

err = immse( y_rep , y);


if write==1

    dlmwrite( fullfile(path,strcat(name,'MSE_error.txt')),err,'delimiter','\t','precision',3)
    dlmwrite( fullfile(path,strcat(name,'A_estimate.txt')),A_estimate,'delimiter','\t','precision',3)
    % dlmwrite( fullfile(path,strcat(name,'J_A.txt')),minVal,'delimiter','\t','precision',3)
    dlmwrite( fullfile(path,strcat(name,'SNR.txt')),SNR,'delimiter','\t','precision',3)
    dlmwrite( fullfile(path,strcat(name,'coef_fromPen.txt')),coef,'delimiter','\t','precision',3)
    dlmwrite( fullfile(path,strcat(name,'_fs.txt')),fs,'delimiter','\t','precision',3)

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
    legend(buff1,buff2,'Location','northeast')
    xlabel('Time t') 
    ylabel('Value') 
    
    
    if write==1
        F = getframe(gcf);
        baseFileName = sprintf('deriv_Y_%d.png', i);
        imwrite(F.cdata, fullfile(path,strcat(name,baseFileName)),'png')
    end 
end


%% LETS TRY USING THE DERIVITATIVES
% 
% y_4=   Y_deri_ture*A_buff(1:4,:);
% plot(t_array,y_4)
% 
% 
% y_4_est= Y_deri*A_estimate;
% plot(t_array,y_4_est)
% 
% 
% A_hat=regress(y_n,Y_deri)
% 
% Y_deri_ture(:,)
% 
% A_hat=regress(y,Y_deri_ture)


