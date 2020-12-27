classdef utilities
    % Matlab doesn't let us to create a library so I made a class with all
    % the functions I will be using this way I can call them in different
    % locattions and as well maintane them in one locations 
    methods (Static)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%start up 

%to do 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function v = K_fy(t,t_array,y_array,A,n,t_i,a)
           % n is the degree of our model set to 4 usually
           % t is the time 
           % t is Theta 
           % a is the  the integral regions boundries 
           % A is the coefficent of the system 
           % t_array is the array of t values 
           % t_i is t  value 
           % y_array is the values of y 
           
           %This funtion calculate the value of the Kernal   KF;y
           y = interp1(t_array,y_array,t);
           buff1=0;
           buff2 = 0;
           for j=1:n
               buff1 = buff1+ ((-1).^(j+1)).* nchoosek(n,j) .* ( factorial(n).*((t_i-t).^(j-1)) .* ((t-a).^(n-j))) /(factorial(n-j).*factorial(j-1));
           end
           %matlab arries are created from 1 to n that is why I added the i+1
           for i = 0:n-1 
               buff3=0;
               for j=0:i
                    buff3 = buff3+((-1).^(j+1)).* nchoosek(i,j) .* ( factorial(n).*((t_i-t).^(n-i+j-1)) .* ((t-a).^(n-j))) /(factorial(n-j).*factorial(n-i+j-1));                    
               end                
               buff2=buff2+ A(i+1) .* buff3;
           end
           v=buff1 +buff2;
           v=v.*y;
%            v=vpa(v); %so it does the divide
       end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       function v = K_By(t,t_array,y_array,A,n,t_i,b)
           % n is the degree of our model set to 4 usually
           % t is the time 
           % T is Theta 
           % b is  the integral regions boundries 
           % A is the coefficent of the system 
           % B is the coefficent of the input 
           % t_i is t  value 
           % y_array is the values of y 
           
           y = interp1(t_array,y_array,t);
           buff1=0;
           buff2 =0;
           for j=1:n
               buff1 = buff1+  nchoosek(n,j) .* ( factorial(n).*((t_i-t).^(j-1)) .* ((b-t).^(n-j))) /(factorial(n-j).*factorial(j-1));
           end
          
           for i = 0:n-1 
               buff3=0;
               for j=0:i
                    buff3 = buff3+ nchoosek(i,j) .* ( factorial(n).*((t_i-t).^(n-i+j-1)) .* ((b-t).^(n-j))) /(factorial(n-j).*factorial(n-i+j-1));                    
               end                
               buff2=buff2+ A(i+1) .* buff3; %matlab arries are created from 1 to n that is why I added the i+1
           end
           v=buff1 +buff2;
           v=v.*y;

       end
       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
       function v = K_fu(t,t_array,u_array,B,n,t_i,a)
           % n is the degree of our model set to 4 usually
           % t is the time 
           % t is Theta 
           % a is the  the integral regions boundries 
           % B is the coefficent of the Input 
           % t_array is the array of t values 
           % t_i is t  value 
           % u_array is the values of u 
           
           %This funtion calculate the value of the Kernal   KF;u
           u = interp1(t_array,u_array,t);
           buff2=0;
           %matlab arries are created from 1 to n that is why I added the i+1
           for i = 0:n-1 
               buff1=0;           
               for j=0:i
                    buff1 = buff1+((-1).^(j+1)).* nchoosek(i,j) .* ( factorial(n).*((t_i-t).^(n-i+j-1)) .* ((t-a).^(n-j))) /(factorial(n-j).*factorial(n-i+j-1));    
               end
               buff2=buff2+ buff1.*B(i);
           end
           v=buff2;
           v=v.*u;
       end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       function v = K_Bu(t,t_array,u_array,B,n,t_i,b)
           % n is the degree of our model set to 4 usually
           % t is the time 
           % t is Theta 
           % a is the  the integral regions boundries 
           % B is the coefficent of the Input 
           % t_array is the array of t values 
           % t_i is t  value 
           % u_array is the values of u 
           
           %This funtion calculate the value of the Kernal   KF;u
           u = interp1(t_array,u_array,t);
           buff2=0;
          
           %matlab arries are created from 1 to n that is why I added the i+1
           for i = 0:n-1 
               buff1=0;                          
               for j=0:i %matlab arries are created from 1 to n that is why I added the i+1
                    buff1 = buff1+ nchoosek(i,j) .* ( factorial(n).*((t_i-t).^(n-i+j-1)) .* ((b-t).^(n-j))) /(factorial(n-j).*factorial(n-i+j-1));                    
               end
                buff2=buff2+ buff1.*B(i);
           end
           v=buff2;
           v=v.*u;

       end
       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function v = K_v_i(t,t_array,y_array,n,i,t_i,a,b)
     % n is the degree of our model set to 4 usually
           % t is the time 
           % t is Theta 
           % a is the  the integral regions boundries 
           % A is the coefficent of the system 
           % t_array is the array of t values 
           % t_i is t  value 
           % y_array is the values of y 
           
           %This funtion calculate the value of the Kernal   Kv
           
           
           
          
    
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function dy = sys_new_4th_order_LTI_y(t, y)
        %this function has been coppied from older code and it is base on
        %this example :
        % % a0=1154.8332;a1=0;a2=156.96;a3=0;

        dy = zeros(4,1);
        dy(1) = y(2);
        dy(2) = y(3);
        dy(3) = y(4);
        dy(4) = -156.96.*y(3) - 1154.8332.*y(1);
    end
    
    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function dy = sys_new_4th_order_LTI_y_general (t,y)
        %I tried to make it a bit more robust 

        A=[-150;-125;-31;-5];
%         A=[-121.6388,-47.8382,-30.0131,-2.4396];
        %         A=[-126.1684;-112.1701;-29.9769;-4.5383];
        dy = zeros(4,1);
        dy(1) = y(2);
        dy(2) = y(3);
        dy(3) = y(4);
        dy(4) = A(4).*y(4)+ A(3).*y(3)+ A(2).*y(2) + A(1).*y(1);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
    function dy = sys_new_4th_order_LTI_y_general_with_U (t,y,U)
        %this Function uses x and U 
        A=[-1154.8332;0.0;-156.96;0.0];
        dy = zeros(4,1);
        dy(1) = y(2);
        dy(2) = y(3);
        dy(3) = y(4);
        dy(4) = A(4).*y(4)+ A(3).*y(3)+ A(2).*y(2) + A(1).*y(1);
    end

% 
%     
%  function dy = sys_new_3th_order_LTI_editable (t,y) %not working 
%         %I tried to make it a bit more robust 
%         % we read the coefficents from a matrix
%         values = load('A_matrix.mat');
%         A=values.A;
%         dy = zeros(3,1);
%         dy(1) = y(2);
%         dy(2) = y(3);
%         dy(3) =  (-1)*A(3).*y(3)+ (-1)*A(2).*y(2) + (-1)*A(1).*y(1);
%         
%         
%  end
%  

%% six order
 function dy = sys_new_6th_order_LTI_editable (t,y) %not working 
        %I tried to make it a bit more robust 
        % we read the coefficents from a matrix
        values = load('A_matrix.mat');
        A=values.A;
        dy = zeros(6,1);
        dy(1) = y(2);
        dy(2) = y(3);
        dy(3) = y(4);
        dy(4) = y(5);
        dy(5)= y(6);
        dy(6)= (-1)* A(6).*y(6)+ (-1)* A(5).*y(5)+(-1)* A(4).*y(4)+ (-1)*A(3).*y(3)+ (-1)*A(2).*y(2) + (-1)*A(1).*y(1);                 
 end



%% fifth order
 function dy = sys_new_5th_order_LTI_editable (t,y) %not working 
        %I tried to make it a bit more robust 
        % we read the coefficents from a matrix
        values = load('A_matrix.mat');
        A=values.A;
        dy = zeros(5,1);
        dy(1) = y(2);
        dy(2) = y(3);
        dy(3) = y(4);
        dy(4) = y(5);
        dy(5)=(-1)* A(5).*y(5)+(-1)* A(4).*y(4)+ (-1)*A(3).*y(3)+ (-1)*A(2).*y(2) + (-1)*A(1).*y(1);                 
 end


%%

     
 function dy = sys_new_3th_order_LTI_editable (t,y) %not working 
        %I tried to make it a bit more robust 
        % we read the coefficents from a matrix
        values = load('A_matrix.mat');
        A=values.A;
        dy = zeros(3,1);
        dy(1) = y(2);
        dy(2) = y(3);
        dy(3) =  (-1)*A(3).*y(3)+ (-1)*A(2).*y(2) + (-1)*A(1).*y(1);
        
        
 end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function dy = sys_new_2th_order_LTI_editable (t,y) %not working 
        %I tried to make it a bit more robust 
        % we read the coefficents from a matrix
        values = load('A_matrix.mat');
        A=values.A;
        dy = zeros(2,1);
        dy(1) = y(2);
        dy(2) =(-1)*A(2).*y(2) + (-1)*A(1).*y(1);
    end
    
 %% 
  function dy = sys_new_4th_order_LTI_editable (t,y) %not working 
        %I tried to make it a bit more robust 
        % we read the coefficents from a matrix
        values = load('A_matrix.mat');
        A=values.A;
        dy = zeros(4,1);
        dy(1) = y(2);
        dy(2) = y(3);
        dy(3) = y(4);
        dy(4) = (-1)* A(4).*y(4)+ (-1)*A(3).*y(3)+ (-1)*A(2).*y(2) + (-1)*A(1).*y(1);
    end
    
 function dy = sys_new_4th_order_LTI_editable_General (t,y) %not working 
        %I tried to make it a bit more robust 
        % we read the coefficents from a matrix
        values = load('A_matrix.mat');
        A=values.A;
        dy = zeros(4,1);
        dy(1) = y(2);
        dy(2) = y(3);
        dy(3) = y(4);
        dy(4) = (-1)* A(4).*y(4)+ (-1)*A(3).*y(3)+ (-1)*A(2).*y(2) + (-1)*A(1).*y(1);
        dy(4)= A(5).* dy(4);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = K_DS_i_Integral(t,t_array,y_array,n,i,t_i,a,b)
           % n is the degree of our model set to 4 usually
           % t is the time 
           % t is Theta 
           % a is the  the integral regions boundries 
           % A is the coefficent of the system 
           % t_array is the array of t values 
           % t_i is t  value 
           % y_array is the values of y 
           
           %This funtion calculate the value of the Kernal   KF;y used for
           %normal  integration
           y = interp1(t_array,y_array,t);
           buff1=utilities.K_DS_i_helper(t,n,i,t_i,a,b);

           v=buff1.*y;
       end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 %% %AIC Calculator    
 function rtn=Aic(y_old,Yhat,n,flag)
     %if flag is 0 return AICc
     SSE=immse(y_old, Yhat); %Sum of Squared Errors for the training set
     k = length(Yhat); % Number of training cases
     p = n ; % Number of parameters (weights and biases)
    
     if (flag==1)         % Schwarz's Bayesian criterion (or BIC) (Schwarz, 1978)
        SBC = k * log(SSE/k) + p * log(k);
        rtn=SBC;
     elseif (flag==2) % Akaike's information criterion (Akaike, 1969)
        AIC = k * log(SSE/k) + 2 * p;
        rtn=AIC;
     else % Corrected AIC (Hurvich and Tsai, 1989)
        AICc = k * log(SSE/k) + (k + p) / (1 - (p + 2) / k);
        rtn=AICc;
     end 
 
 end
 %% 
%%%%%%%%%%%



function v = K_DS_i_trapz(n,i,t_array,y_array,t_i,a,b)
           % n is the degree of our model set to 4 usually
           % t is the time 
           % t is Theta 
           % a is the  the integral regions boundries 
           % A is the coefficent of the system 
           % t_array is the array of t values 
           % t_i is t  value 
           % y_array is the values of y 
           
           %This funtion calculate the value of the Kernal   KF;y used for
           %Trapz integration
           
           buffer=zeros(length(y_array),1);
           for j = 1:1:size(y_array)  
               t= t_array(j);  %#T
               buffer(j)=utilities.K_DS_i_helper(t,n,i,t_i,a,b)* y_array(j);
           end
           
           v=buffer;
       end

%            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
function v = K_DS_i_helper(t,n,i,t_i,a,b)
           % n is the degree of our model set to 4 usually
           % t is the time 
           % t_i is Theta 
           % a is the  the integral regions boundries 
           % A is the coefficent of the system 
           % t_array is the array of t values 
           % t_i is t  value 
           % y_array is the values of y 
           
           %This funtion calculate the value of the Kernal   KF;y used for
           %Trapz integration
                      
           buff1=0;
           if(t<=t_i) %KFy_i
                if (i==n)
                   for j=1:n
                       buff1 = buff1+ ((-1).^(j+1)).* nchoosek(n,j) .* ( factorial(n).*((t_i-t).^(j-1)) .* ((t-a).^(n-j))) /(factorial(n-j).*factorial(j-1));
                   end                         
               elseif(i<0)
                   disp('error: i is negative');
               elseif(i < n)     
                   for j=0:i
                        buff1 = buff1+((-1).^(j+1)).* nchoosek(i,j) .* ( factorial(n).*((t_i-t).^(n-i+j-1)) .* ((t-a).^(n-j))) /(factorial(n-j).*factorial(n-i+j-1));    
                   end                   
               else
                   disp('error');              
               end
           else %KBY_i
               if (i==n)
                   for j=1:n
                        buff1 = buff1+  nchoosek(n,j) .* ( factorial(n).*((t_i-t).^(j-1)) .* ((b-t).^(n-j))) /(factorial(n-j).*factorial(j-1));
                   end                         
               %should also check if I is negative
               elseif(i<0)
                   disp('error: i is negative');
               elseif(i < n)     
                   for j=0:i %matlab arries are created from 1 to n that is why I added the i+1
                        buff1 = buff1+ nchoosek(i,j) .* ( factorial(n).*((t_i-t).^(n-i+j-1)) .* ((b-t).^(n-j))) /(factorial(n-j).*factorial(n-i+j-1));                    
                   end     
               else
                   disp('error');              
               end               
           end
           v=buff1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = K_DS_i_trapz_repo(n,i,t_array,y_array,t_i,a,b,A)
           % n is the degree of our model set to 4 usually
           %WE have A matrix known !!
           % t is the time 
           % t is Theta 
           % a is the  the integral regions boundries 
           % A is the coefficent of the system 
           % t_array is the array of t values 
           % t_i is t  value 
           % y_array is the values of y 
           
           %This funtion calculate the value of the Kernal   KF;y used for
           %Trapz integration
           
           buffer=zeros(length(y_array),1);
           for j = 1:1:size(y_array)  
               t= t_array(j);  %#T
               buffer(j)=utilities.K_DS_i_helper_repo(t,n,i,t_i,a,b,A)* y_array(j);
           end
           
           v=buffer;
       end

%            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
function v = K_DS_i_helper_repo(t,n,i,t_i,a,b,A)
           % n is the degree of our model set to 4 usually
           % t is the time 
           % t_i is Theta 
           % a is the  the integral regions boundries 
           % A is the coefficent of the system 
           % t_array is the array of t values 
           % t_i is t  value 
           % y_array is the values of y 
           
           %This funtion calculate the value of the Kernal   KF;y used for
           %Trapz integration
                      
           buff1=0;
           if(t<=t_i) %KFy_i
                if (i==n)
                   for j=1:n
                       buff1 = buff1+ ((-1).^(j+1)).* nchoosek(n,j) .* ( factorial(n).*((t_i-t).^(j-1)) .* ((t-a).^(n-j))) /(factorial(n-j).*factorial(j-1));
                   end                         
               elseif(i<0)
                   disp('error: i is negative');
               elseif(i < n)     
                   for j=0:i
                        buff1 = buff1+((-1).^(j+1)).* nchoosek(i,j) .* ( factorial(n).*((t_i-t).^(n-i+j-1)) .* ((t-a).^(n-j))) /(factorial(n-j).*factorial(n-i+j-1));    
                        buff1= buff1.*A(j+1); %%%%%%%%%%%%%%%%%%%%%%%%ADDED THIS
                   end    
                   
               else
                   disp('error');              
               end
           else %KBY_i
               if (i==n)
                   for j=1:n
                        buff1 = buff1+  nchoosek(n,j) .* ( factorial(n).*((t_i-t).^(j-1)) .* ((b-t).^(n-j))) /(factorial(n-j).*factorial(j-1));
                   end                         
               %should also check if I is negative
               elseif(i<0)
                   disp('error: i is negative');
               elseif(i < n)     
                   for j=0:i %matlab arries are created from 1 to n that is why I added the i+1
                        buff1 = buff1+ nchoosek(i,j) .* ( factorial(n).*((t_i-t).^(n-i+j-1)) .* ((b-t).^(n-j))) /(factorial(n-j).*factorial(n-i+j-1));                    
                        buff1= buff1.*A(j+1);%%%%%%%%%%%%%%%%%%%%%%%%ADDED THIS
                   end     
               else
                   disp('error');              
               end               
           end
           v=buff1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%%%%%%%%%%%
function v = K_DS_i_trapz_small_horizen(n,i,t_array,y_array,t_i,a,b)
           % n is the degree of our model set to 4 usually
           % t is the time 
           % t is Theta 
           % a is the  the integral regions boundries 
           % A is the coefficent of the system 
           % t_array is the array of t values 
           % t_i is t  value 
           % y_array is the values of y 
           
           %This funtion calculate the value of the Kernal   KF;y used for
           %Trapz integration
           buffer=zeros(length(y_array),1);
           for j = 1:1:size(y_array)  
               t = t_array(j);  %#T
               if (t<6)
                   a=1;
                   b=10;
               elseif(t>size(y_array)-5)
                   a= size(y_array)-5; 
                   b= size(y_array);
               else
                   a=t-5;
                   b=t+5;                                
               end
               
               buffer(j)=utilities.K_DS_i_helper(t,n,i,t_i,a,b)* y_array(j);
           end
  
          
           v=buffer;
       end

%            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    















function v = K_DS_i_trapz_deri(n,k,t_array,y_array,t_i,a,b, A)  %TO_DO no A
           % n is the degree of our model set to 4 usually
           % t is the time 
           % t is Theta 
           % a is the  the integral regions boundries 
           % A is the coefficent of the system 
           % t_array is the array of t values 
           % t_i is t  value 
           % y_array is the values of y 
           
           %This funtion calculate the value of the Kernal   KF;y used for
           %Trapz integration
           
           buffer=zeros(length(y_array),1);
           for j = 1:1:size(y_array)  
               t= t_array(j);  %#T
               buffer(j)=utilities.K_DS_i_derivative_helper(t,n,k,t_i,a,b,A)* y_array(j);
           end
           
           v=buffer;
       end

%            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %problem with i 
    
function v = K_DS_i_derivative_helper(t,n,k,t_i,a,b,A)
           % n is the degree of our model set to 4 usually
           % t is the time 
           % t is Theta 
           % a is the  the integral regions boundries 
           % A is the coefficent of the system 
           % t_array is the array of t values 
           % t_i is t  value 
           % y_array is the values of y 
           % K is the k-th derivastive 
           
           %This funtion calculate the value of the Kernal   KF;y used for
           
           
           p= n-k;
           if p<=0
               p=n;
               disp('K has to be smaller than n');
           end 
  
           
           buff1=0;
           if(t<=t_i) %KFy_i
               
               for j=1:p
                    buff1 = buff1+ ((-1).^(j+k+1)).* nchoosek(n,k+j) .* ( factorial(n).*((t_i-t).^(j-1)) .* ((t-a).^(p-j))) /(factorial(p-j).*factorial(j-1));
               end  
               
               buff2=0;
               for i=0:(p-1)
                   a=A(i+1);
                   dummy=0;
                   for j=0:i
                        dummy = dummy+((-1).^(j+1)).* nchoosek(i,j) .* ( factorial(n).*((t_i-t).^(p-i+j-1)) .* ((t-a).^(n-j))) /(factorial(n-j).*factorial(p-i+j-1));    
                   end 
                   buff2=buff2+a.*dummy;
               end
               buff1=buff1+buff2;
                        
               buff3=0;
               for i=p:(n-1)
                   a=A(i+1);
                   dummy=0;
                   for j=1:p
                        dummy = dummy+((-1).^(j+i-p+1)).* nchoosek(i,i-p+j) .* ( factorial(n).*((t_i-t).^(j-1)) .* ((t-a).^(n-i+p-j))) /(factorial(n-i+p-j).*factorial(j-1));    
                   end 
                   buff3=buff3+a.*dummy;
               end
               
               buff1=buff1+buff3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
           else %KBY_i
               for j=1:p
                    buff1 = buff1+  nchoosek(n,k+j) .* ( factorial(n).*((t_i-t).^(j-1)) .* ((b-t).^(p-j))) /(factorial(p-j).*factorial(j-1));
               end  
               
               buff2=0;
               for i=0:(p-1)
                   a=A(i+1);
                   dummy=0;
                   for j=0:i
                        dummy = dummy+ nchoosek(i,j) .* ( factorial(n).*((t_i-t).^(p-i+j-1)) .* ((b-t).^(n-j))) /(factorial(n-j).*factorial(p-i+j-1));    
                   end 
                   buff2=buff2+a.*dummy;
               end
               buff1=buff1+buff2;
                        
               buff3=0;
               for i=p:(n-1)
                   a=A(i+1);
                   dummy=0;
                   for j=1:p
                        dummy = dummy+ nchoosek(i,i-p+j) .* ( factorial(n).*((t_i-t).^(j-1)) .* ((b-t).^(n-i+p-j))) /(factorial(n-i+p-j).*factorial(j-1));    
                   end 
                   buff3=buff3+a.*dummy;
               end
               
               buff1=buff1+buff3;
                   
                   
                   
                   
                   
                   
                   

           end
           
           
           
           v=buff1;
       end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
function v = K_fy_i(t,t_array,y_array,n,i,t_i,a)
           % n is the degree of our model set to 4 usually
           % t is the time 
           % t is Theta 
           % a is the  the integral regions boundries 
           % A is the coefficent of the system 
           % t_array is the array of t values 
           % t_i is t  value 
           % y_array is the values of y 
           
           %This funtion calculate the value of the Kernal   KF;y
           y = interp1(t_array,y_array,t);
           buff1=0;
           if (i==n)
               for j=1:n
                   buff1 = buff1+ ((-1).^(j+1)).* nchoosek(n,j) .* ( factorial(n).*((t_i-t).^(j-1)) .* ((t-a).^(n-j))) /(factorial(n-j).*factorial(j-1));
               end                         
           elseif(i<0)
               disp('error: i is negative');
           elseif(i < n)     
               for j=0:i
                    buff1 = buff1+((-1).^(j+1)).* nchoosek(i,j) .* ( factorial(n).*((t_i-t).^(n-i+j-1)) .* ((t-a).^(n-j))) /(factorial(n-j).*factorial(n-i+j-1));    
               end                   
           else
               disp('error');              
           end

           v=buff1.*y;
       end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function v = K_By_i(t,t_array,y_array,n,i,t_i,b)
           % n is the degree of our model set to 4 usually
           % t is the time 
           % T is Theta 
           % b is  the integral regions boundries 
           % A is the coefficent of the system 
           % B is the coefficent of the input 
           % t_i is t  value 
           % y_array is the values of y 
           %i is the i-th order used for gi 
           buff1=0;
           y = interp1(t_array,y_array,t);           
           if (i==n)
               for j=1:n
                    buff1 = buff1+  nchoosek(n,j) .* ( factorial(n).*((t_i-t).^(j-1)) .* ((b-t).^(n-j))) /(factorial(n-j).*factorial(j-1));
               end                         
           %should also check if I is negative
           elseif(i<0)
               disp('error: i is negative');
           elseif(i < n)     
               for j=0:i %matlab arries are created from 1 to n that is why I added the i+1
                    buff1 = buff1+ nchoosek(i,j) .* ( factorial(n).*((t_i-t).^(n-i+j-1)) .* ((b-t).^(n-j))) /(factorial(n-j).*factorial(n-i+j-1));                    
               end     
           else
               disp('error');              
           end
           v=buff1.*y;
       end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y_noisy = White_noise_adder(y,M,Sdev)
            % %addition of noise%
            % This function add white noise (Guassian)
            % y is the signal to be passed in 
            % M is the mean of the noise (should be zero)
            % Sdev is the standard devation of the noise
            WN_noise = random('norm',M,Sdev,  length(y),1); % mean zero , STD=2.5
            y_noisy = y+WN_noise;
        end 
        
        
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     

    function y_est=EvalFunction(t_array,a,b,n,y_n,max_count,A_buff)
            % this function uses equation 12 to estimate Y
            %A_buff  A is known 
            % Eq 9 assumption is that U is zero here
            
            % (a,b) integral boundaries is known
            % t_array is the time array 
            % n is the order of the system 
            %  y_n
            % max_count
            G_mat=zeros(length(A_buff),1);
            y_est=zeros(length(y_n),1);


                %  do a dot product 
            tic %timer 
            for i = 1:1:length(t_array)
                t_i = t_array(i);    
                alfa_function=((t_i-a).^n +(b-t_i).^n);
                for j=0:n  
                    G_mat(j+1)=( quadgk(@(t,v)utilities.K_fy_i(t,t_array,y_n,n,j,t_i,a),a,t_i,'MaxIntervalCount',max_count) + quadgk(@(t,v)utilities.K_By_i(t,t_array,y_n,n,j,t_i,b),t_i,b,'MaxIntervalCount',max_count));
                end
                y_est(i) = dot(A_buff,G_mat) /alfa_function; 
            end
            toc
    end
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%inspried by https://github.com/makintunde/gram-schmidt
%https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process
%http://web.mit.edu/18.06/www/Essays/gramschmidtmat.pdf

    function[result] = gram_schmidt(A, n)
    % MATLAB implementation of Classical Gram-Schmidt Algorithm.
    % A is the matrix 
    % n is the order of the system 
    R=zeros( size(A) );
    Q=zeros( size(A) );

    
    for j=1:n
       v= A(:, j); 
       for i=1:j-1
           R(i,j)=Q(:,i)' * A(:,j);  %transpose and multiply
           v=v-R(i,j)*Q(:,i);
       end 
        R(j,j)=norm(v);
        Q(:,j)=v/R(j,j);
    end 
    
    
    
    result = Q;
end
 function[result] = gram_schmidt_L2(A, n)
    %L2 normed MATLAB implementation of Classical Gram-Schmidt Algorithm. 
    % A is the matrix 
    % n is the order of the system 
    R=zeros( size(A) );
    Q=zeros( size(A) );

    
    for j=1:n
       v= A(:, j); 
       for i=1:j-1
           R(i,j)=Q(:,i)' * A(:,j);  %transpose and multiply
           v=v-R(i,j)*Q(:,i);
       end 
        R(j,j)=norm(v);
        Q(:,j)=v/R(j,j);
    end 
    
    
    
    result = Q;
 end
%% 

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
 function v = K_fy_K_deri(t,t_array,y_array,A,n,t_i,a,k)
           % n is the degree of our model set to 4 usually
           % t is the time 
           % t is Theta 
           % a is the  the integral regions boundries 
           % A is the coefficent of the system 
           % t_array is the array of t values 
           % t_i is t  value 
           % y_array is the values of y 
           
           %This funtion calculate the value of the Kernal   KF;y
           y = interp1(t_array,y_array,t);
           p=n-k;           
           buff1=0;
           buff2 = 0;
           buff4=0;
           
           for j=1:p
               buff1 = buff1+ ((-1).^(j+n-p+1)).* nchoosek(n,n-p+j) .* ( factorial(n).*((t_i-t).^(j-1)) .* ((t-a).^(p-j))) /(factorial(p-j).*factorial(j-1));
           end
           %matlab arries are created from 1 to n that is why I added the i+1
           for i = 0:p-1 
               buff3=0;
               for j=0:i
                    buff3 = buff3+((-1).^(j+1)).* nchoosek(i,j) .* ( factorial(n).*((t_i-t).^(p-i+j-1)) .* ((t-a).^(n-j))) /(factorial(n-j).*factorial(p-i+j-1));                    
               end                
               buff2=buff2+ A(i+1) .* buff3;
           end
           
          for i = p:n-1 
               buff5=0;
               for j=1:p
                    buff5 = buff5+((-1).^(j+i-p+1)).* nchoosek(i,i-p+j) .* ( factorial(n).*((t_i-t).^(j-1)) .* ((t-a).^(n-i+p-j))) /(factorial(n-i+p-j).*factorial(j-1));                    
               end                
               buff4=buff4+ A(i+1) .* buff5;
           end
           
           v=buff1 +buff2+buff4;
           v=v.*y;           
           
       end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function v = K_By_i_deri(t,t_array,y_array,A,n,t_i,b,k)
           % n is the degree of our model set to 4 usually
           % t is the time 
           % T is Theta 
           % b is  the integral regions boundries 
           % A is the coefficent of the system 
           % B is the coefficent of the input 
           % t_i is t  value 
           % y_array is the values of y 
           %i is the i-th order used for gi 
           p=n-k;           
           buff1=0;
           buff2 = 0;
           buff4=0;
           y = interp1(t_array,y_array,t);           
           for j=1:p
               buff1 = buff1+  nchoosek(n,n-p+j) .* ( factorial(n).*((t_i-t).^(j-1)) .* ((b-t).^(p-j))) /(factorial(p-j).*factorial(j-1));
           end
          
           for i = 0:p-1 
               buff3=0;
               for j=0:i
                    buff3 = buff3+ nchoosek(i,j) .* ( factorial(n).*((t_i-t).^(p-i+j-1)) .* ((b-t).^(n-j))) /(factorial(n-j).*factorial(p-i+j-1));                    
               end                
               buff2=buff2+ A(i+1) .* buff3; %matlab arries are created from 1 to n that is why I added the i+1
           end
           
           for i = p:n-1 
               buff5=0;
               for j=1:p
                    buff5 = buff5+ nchoosek(i,i-p+j) .* ( factorial(n).*((t_i-t).^(j-1)) .* ((b-t).^(n-i+p-j))) /(factorial(n-i+p-j).*factorial(j-1));                    
               end                
               buff4=buff4+ A(i+1) .* buff5; %matlab arries are created from 1 to n that is why I added the i+1
           end
           
           v=buff1 +buff2+buff4;
           v=v.*y;
       end
       
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        function v=AICc
           
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function v = deriv_helper(t_array,y_array,A,n,a,b,k,Y_deri)
           % n is the degree of our model set to 4 usually
           % t is the time 
           % T is Theta 
           % b is  the integral regions boundries 
           % A is the coefficent of the system 
           % B is the coefficent of the input 
           % t_i is t  value 
           % y_array is the values of y 
           %i is the i-th order used for gi 
           % Y_deri is the matrix hoding the derivatives %
           % Y_deri(t_i+1,k-j+1) (row ,column) the columns are shifted by 1 in matlab ,  starts at 1 
           % % Y_deri(t_i+1,1)= y0(t_i)
           p=n-k;           
           y_dummy =zeros(length(y_array),1);
           y_n=y_array;
           
           parfor c = 1:1:length(t_array) %PARALLEL 
                if(mod(c,100)==0)
                        fprintf('%d\n',c);
                end
                t_i = t_array(c);    %t_i is the t in the notes
                alfa_function=((t_i-a).^n +(b-t_i).^n);
                
                % a
                buff1=0;
                for j=1:k  %in the notes it says i but we use j 
                        buff1 = buff1+ (((-1).^(j+1)).* nchoosek(p+j-1,j) .* (factorial(n).*((t_i-a).^(n-j))) /(factorial(n-j))).* Y_deri(c,k-j+1); 
                end
                
                buff2 = 0;
                
                
                for i = p:n-1 
                       buff3=0;
                       for j=0:i-p
                            buff3 = buff3+(((-1).^(j+1)).* nchoosek(p+j-1,j) .* ( factorial(n).*((t_i-a).^(n-j))) /(factorial(n-j))).*Y_deri(c,i-j-p+1);                    
                       end                
                       buff2=buff2+ A(i+1) .* buff3;
                end
                
                
                % b 
                buff4=0;
                for j=1:k  %in the notes it says i but we use j 
                        buff4 = buff4+ ( nchoosek(p+j-1,j) .* (factorial(n).*((b-t_i).^(n-j))) /(factorial(n-j))).* Y_deri(c,k-j+1); 
                end
                
                buff5=0;
                for i = p:n-1 
                       buff6=0;
                       for j=0:i-p
                            buff6 = buff6+( nchoosek(p+j-1,j) .* ( factorial(n).*((b-t_i).^(n-j))) /(factorial(n-j))).*Y_deri(c,i-j-p+1);                    
                       end                
                       buff5=buff5+ A(i+1) .* buff6;
                end
                %
                dummy= buff1+buff2-buff4-buff5;
                y_dummy(c)=(dummy+ integral(@(t,v)utilities.K_fy_K_deri(t,t_array,y_n,A,n,t_i,a,k),a,t_i,'RelTol',1e-2,'AbsTol',1e-6) + integral(@(t,v)utilities.K_By_i_deri(t,t_array,y_n,A,n,t_i,b,k),t_i,b,'RelTol',1e-2,'AbsTol',1e-6))/alfa_function;
           end
                      
           v=y_dummy;                                         
       end
       

    
        end % static methods
    end % classdef


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

