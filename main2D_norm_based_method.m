
 clc;clear;close all;

%% Initialization



% initial state
n = 2;

% X0_1 = [0.9,1.1];
% X0_2 = [-0.1,0.1];

x0 = [10.5,-5.5]';
xc = [10,-5]';
xb = 1;

% input signal
Tend = 20;
step = 0.01;
t = 0:step:Tend;

U = 0.5;
figure()
hold on
plot(t,U*sin(t))
xlabel('time (s)')
ylabel('input')
% legend('u')
title('Input signal')
hold off
A = [-1 -4; 4, -1];
B = [1;1]; 
C = [1 0];
sys = ss(A,B,C,0);


% Ud = -0.3+0.6*rand(length(t),1);
ud = 0.05;
Ud = -ud+2*ud*rand(length(t),1);

Uu = U*sin(t)';
% plant trajectory

tspan = [0 Tend];

[y_ode,TTTT,x_ode] = lsim(sys,Uu+Ud,t,x0);
% [T_ode,x_ode] = ode45(@(t,x) nlinearF(t,x,A,H), tspan, x0);

x_ode = x_ode';
figure()
hold on
plot(t,x_ode(1,:),'color',[.5 0 .5])
plot(t,x_ode(2,:),'color','green')
xlabel('time')
ylabel('state')
legend('$x_1$','$x_2$','$x_3$','Interpreter','latex','FontSize',20)
title('Plant trajectories')
hold off

%% solve LMI for P and nu_1 nu_2
iii=1;
% for iter = 0:0.001:0.80
% E= [1;1];
% gamma_old = 0;
tune = 7.90;
tune2 = 20;
% for iii = 1:100

setlmis([]) 
p = lmivar(1,[2 1]);
nu1 = lmivar(1,[1 0]);
nu2 = lmivar(1,[1 0]);
pk = lmivar(2,[2 1]);

lmiterm([1 1 1 p],1,A,'s');   % A^TP+ PA


lmiterm([1 1 1 pk],1,C,'s');   % A^TP+ PA


lmiterm([1 1 1 nu1],1,eye(2));   
lmiterm([1 2 1 p],1,1);   
lmiterm([1 1 2 p],1,1);  


lmiterm([1 2 2 nu2],-1,eye(2));
lmiterm([-2 1 1 p],1,1);         % P>0
lmiterm([2 1 1 0],0)            % P>0


lmiterm([-3 1 1 nu1],1,1);      % nu1>0
lmiterm([3 1 1 0],tune);      % nu1>tune


lmiterm([ -4 1 1 nu2],1,1);      % nu2>0

lmiterm([5 1 1 nu2],1,1);      % nu2 < tune2
lmiterm([-5 1 1 0],tune2);      %

% lmiterm([5 1 1 0],0);      % nu2>0



lmis = getlmis;
% options = [1e-5,0,0,0,0] 
% [copt,xopt] = mincx(lmis,[0 0 0 1 1],options)

options = [0,0,14,0,0];
[tmin,xfeas] = feasp(lmis,options,0);


P = dec2mat(lmis,xfeas,p);
Nu1 = dec2mat(lmis,xfeas,nu1);
Nu2 = dec2mat(lmis,xfeas,nu2);
PK = dec2mat(lmis,xfeas,pk);

K = inv(P)*PK;
L = -K;

A_bar = A + K*C;

% Nu1=100
% Nu2 = 10000
%check LMI
LMI = [A_bar'*P + P*A_bar + Nu1*eye(2) P;P -Nu2*eye(2)];
isnegdef = all(eig(LMI)<0);
check  = 0;
if isnegdef ~= 1
    warning('LMI failed')
    check = 1;
%     break
end 



E_p = eig(P);
Lam_max = E_p(2);
Lam_min = E_p(1);
 
Lam_e = Nu1./Lam_max/n;


gamma_hat = sqrt(n*Nu2./Lam_min./Lam_e)*ud


% tune = tune + 70*(gamma_hat-gamma_old);
% tune2 = tune2 - 50*(gamma_hat-gamma_old);
% gamma_old = gamma_hat;
% 
% gamma_save(iii)=gamma_hat + check;
% iii = iii+1;
% end

% end

% plot(1:length(gamma_save),gamma_save)
%%

% assume model used for observer is identical to the plant
Aobs = A;
Bobs = B;
Cobs = C;
% obtain both plant and observer in a single system called sysNew
% first 2 state are plant states, last 2 states are observer states
Anew = [A      zeros(2);
        L*C    Aobs-L*Cobs]; % output of the plant affects the observer 
                             %    states due to the way observers work
%first input to the system with observer is the unknown disturbance, 
%   which only affects the plant not the observer, 
%   second input is the control input
%   both affect the plant the same way (i.e.,force acting on the mass)
Bnew = [B          B
        zeros(2,1) Bobs]; 
% first output is plant output,  
%    last 2  are observer state estimates    
Cnew = [C          zeros(1,2); 
        zeros(2,2) eye(2)];
sysNew = ss(Anew,Bnew,Cnew,0);


A_bar =  Aobs-L*Cobs;

figure()
hold on
[y,~,xobs] = lsim(sysNew,[Ud Uu],t,[x0;xc]);
plot(t,xobs(:,1),t,xobs(:,3)) % first state of the observer is the position of the mass
legend('plant output','observer output estimate')
hold off

figure()
hold on
[y,~,xobs] = lsim(sysNew,[Ud Uu],t,[x0;xc]);
plot(t,xobs(:,2),t,xobs(:,4)) % first state of the observer is the position of the mass
legend('plant output','observer output estimate')
hold off





for i=1:length(xobs(:,1))
    err_o_1(i) = abs(xobs(i,1)-xobs(i,3));
    err_o_2(i) = abs(xobs(i,2)-xobs(i,4));
%     beta_hat(i) = sqrt(2*Lam_max./Lam_min)*exp(-Lam_e*t(i)./2)*xb;    % 0.1 from x_b
end

beta_hat = @(t) sqrt(2*Lam_max./Lam_min)*exp(-Lam_e*t./2)*xb;


figure()
hold on
plot(t,err_o_1)
plot(t,err_o_2)
plot(t,beta_hat(t)+gamma_hat)
legend('$|\hat{e}_{1}|$','$|\hat{e}_{2}|$','$\hat{\beta}(x_b,t)+\hat{\gamma}(d_b)$','Interpreter','latex','FontSize',20)
xlim([0,6])
ylim([0,2])
xlabel('time','FontSize',15)
ylabel('estimate error','FontSize',15)



 
xobs = xobs';
x_obs = xobs(3:4,:);
 %% test for quantization


 

% S = rand(1)*2-1
% C = 0;
% L = 1;
% N  = 100;
% Se  = Quanti(S,C,L,N,0)
% Sd = Quanti(Se,C,L,N,1)

%% Initial quantization parameters

 
C_0 = xc;  
L_0 = [1;1]*xb; 

N = 4; %Quantization level
% N=20;
Tq = 0.1; % time between each transmission 

% built zonotope

% a zonotope can be defined by a center vector c and a projection matrix G
C_Z{1} = C_0;

G_Z{1} =eye(2)*xb;

%% start simulating
close all
% Ck(:,1) = C_0;
Cd(:,1) = C_0;  % docoded packet


% for j = 1:n
%     G_Zr{1} = reduceorder(G_Z{1},n);
%     Lk1(j,1)= C_Z{1}(j) + G_Zr{1}(j,j);    % initial quantization region's element-wise upper bound  
%     Lk2(j,1)= C_Z{1}(j) - G_Zr{1}(j,j);    % lower bound
% end

Cq(:,1) = C_0;
Lk(:,1) = L_0;

xk(:,1) = x0;       % plant state
xs(:,1) = x0;      %actual state
T_total = 0:Tq:Tend; 

% store
T_e = cell(0);
T_d = cell(0);
T_k = cell(0);

x = cell(0);
xp = cell(0);
xe = cell(0);
xr = cell(0);
xL = cell(0);
Tr_plot = [];
Xr_plot = [];


figure()
hold on
xlabel('$\hat{x}_1$','FontSize',15,'Interpreter','latex')
ylabel('$\hat{x}_2$','FontSize',15,'Interpreter','latex')
title('Dynamic Quantization Scheme using Hypercubes','FontSize',15)

Tsim = step;


maxplot = 20;
for k = 1:length(T_total)-1
% for k = 1:10


%     [To;9_k{k},x{k}] = ode45(@(t,x) nlinearFode(t,x,A,H,n), [(k-1)*Tq k*Tq], [xk(:,k);Cd(:,k)]);
    
    T_k{k} = (k-1)*Tq:Tsim:k*Tq;
    xp{k} = x_obs(:,round((k-1).*Tq./step)+1:round(k.*Tq./step));           %store plant trajc
    xr{k} = lsim(sys,0*T_k{k},0:Tsim:Tq,Cd(:,k));                        % store reconstucted state
%     xp{k} = xp{k}';
    xr{k} = xr{k}';
    
%  plot over-approximated quantization region
    if k<maxplot
       p(1) = plot([xp{k}(1,:), x_obs(1,round(k.*Tq./step)+1)] ,[xp{k}(2,:) x_obs(2,round(k.*Tq./step)+1)],'-g');
       p(2) =  scatter(xp{k}(1,1),xp{k}(2,1),40,[0 0 0],'filled');

        G_Z{k} = [Lk(1,k) 0; 0 Lk(2,k)];
        C_Z{k} = Cq(:,k);
        ppppp=PlotZonotope(G_Z{k},C_Z{k});
    end


    % update state for next ode's initial plant state
     xk(:,k+1) =  xp{k}(:,end);  
     xs(:,k+1) =  xp{k}(:,1); 
    
%     Lk1(:,k+1) = 1./N*(exp(Lx*Tq).*Lk1(:,k)) + xL{k}(1:n,end);
% 
%     Lk2(:,k+1) = 1./N*(exp(Lx*Tq).*Lk2(:,k)) + xL{k}(n+1:2*n,end);
        
    
    %begin quantization
    
    Se = Quanti(xp{k}(:,1),Cq(:,k),Lk(:,k),N,0);
    Sd = Quanti(Se,Cq(:,k),Lk(:,k),N,1);
    Sd = Sd';
    % Update decoded packet for updating the next quantization centroid
    Cd(:,k+1) =  Sd;
    
%     check for overflow
        Se = Se';
    if max(Se) >= N
        warning('overflooting') 
    end
    if min(Se) < 0
        warning('overflooting')
   
    end
    
%     store overflow index
    overflow(:,k) = Se ;

%         plot quantization subregion and decoded state

    if k==1||k == 2||k==4||k==6||k==8
        [x,y] = meshgrid(C_Z{k}(1)-G_Z{k}(1,1)+ 2*(0:N).*(G_Z{k}(1,1))./N, C_Z{k}(2)-G_Z{k}(2,2)+ 2*(0:N).*(G_Z{k}(2,2))./N);
        pp = plot(x(1),y(1),'--','Color','r');
        plot(x,y,'--','Color','r')
plot(x',y','--','Color','r')

    end
    
     if k<maxplot
  ppp = scatter(Cd(1,k+1),Cd(2,k+1),40,[0 0 1],'filled');  
     end

    % calculate the beta functions
    Beta = (exp(Tq*norm(A,'inf'))-1)./norm(A,'inf').*(U);
    Betad = beta_hat(k*Tq)+gamma_hat;
    
    Betae = (exp(Tq*norm(A,'inf'))-1)./norm(A,'inf').*(norm(K*C,'inf')*Betad);
    Lam = expm(A.*Tq);
    LamHAT = abs(Lam);
    
    
     % Zk+1 = Lam*Zk + Z_beta
    
 
     

      
     
     
     %overapproximating the terminal reachable set using zonotope
    % and update quantization parameters
    

    
    Cq(:,k+1) = Lam*Cd(:,k+1);
    for i = 1:n
        Lk(i,k+1) = Beta+Betae+ exp(norm(A,'inf').*Tq)*Lk(i,k)./N;
    end
        
    
    %      plot the terminal reachable set zonotope
    if k<maxplot-1
        Czz = Lam*Sd;
%         Gzz = Lam/N*G_Z{k};
        Gzz = [(Beta+Betae+ exp(norm(A,'inf').*Tq)*Lk(i,k)./N)*eye(2)];
%         pppppp=PlotZonotope2(Gzz,Czz);
    end

      
    



    
%     G_Z{i} = reduceorder(G_Z{i},8);

    %store reconstruced state
    
%     plot(T_k{k},xr{k},'Color','red')
    Tr_plot  = [Tr_plot,T_k{k}];
    Xr_plot = [Xr_plot,xr{k}];
end
% txt = 'k = 1';
% text(1,1.75,txt,'FontSize',15)
% 
% txt = 'k = 3';
% text(0.2-0.05,1,txt,'FontSize',15)
% 
% txt = 'k = 4';
% text(-0.1,0.9 ,txt,'FontSize',15)
% 
% txt = 'k = 5';
% text(-0.3,0.8 ,txt,'FontSize',15)

legend([p(1) p(2) ppppp pp ppp],'local state estimate $\hat{x}(t)$','local state estimate at the transmission time $\hat{x}(t_k)$','quantization region $S_Q^k$','quantization subregions $S_q^k$','decoded state $P_d^k$','FontSize',15,'Interpreter','latex')
%%
% plot quantization level
figure()
ne = 10;

subplot(2,1,1)
hold on
title('Quantization region for each state','FontSize',15,'Interpreter','latex')
ylabel('$\hat{x}_1$','FontSize',15,'Interpreter','latex')
p(1) = plot(0:0.01:10*ne*0.01,x_obs(1,1:ne*10+1),'Color',[.5 0 .5]);
for i = 1:ne+1
     p(2) = line([0.1*(i-1),0.1*(i-1)],[Cq(1,i) + Lk(1,i) , Cq(1,i)-Lk(1,i)],'Color','red');
    for j =0:N
       line( [0.1*(i-1)-0.01 ,0.1*(i-1)+0.01], [ Cq(1,i)-Lk(1,i)+2*j/N*Lk(1,i) Cq(1,i)-Lk(1,i)+2*j/N*Lk(1,i)]   ,'Color','red'    )
    end
    p(3) = scatter(0.1*(i-2),Cd(1,i),10,[0 0 1],'filled');
end
scatter(0.1*(i-1),Cd(1,i+1),10,[0 0 1],'filled')

legend([p(1) p(2) p(3)],'local state estiamte $\hat{x}_1$','quantization range $[C_1-L_1,C_1+L_1]$','deocoded state $P_{d,1}$','FontSize',11,'Interpreter','latex')

xlim([0.07,1.01])
% semilogy(xplot,(err_q(1,1+1:ne+1)),'-*r')
xlabel('time','FontSize',15)
hold off

subplot(2,1,2)
hold on
ylabel('$\hat{x}_2$','FontSize',15,'Interpreter','latex')
p(1) = plot(0:0.01:10*ne*0.01,x_obs(2,1:ne*10+1),'Color','green');
% semilogy(xplot,(err_q(2,1+1:ne+1)),'-*r')
xlabel('time','FontSize',15)
for i = 1:ne+1
     p(2) = line( [0.1*(i-1),0.1*(i-1)],[Cq(2,i) + Lk(2,i) , Cq(2,i)-Lk(2,i)],'Color','red' );
        for j =0:N
        line( [0.1*(i-1)-0.01 ,0.1*(i-1)+0.01], [ Cq(2,i)-Lk(2,i)+2*j/N*Lk(2,i) Cq(2,i)-Lk(2,i)+2*j/N*Lk(2,i)]   ,'Color','red'    )
        end
   p(3)=  scatter(0.1*(i-2),Cd(2,i),10,[0 0 1],'filled');
end
scatter(0.1*(i-1),Cd(2,i+1),10,[0 0 1],'filled')
legend([p(1) p(2) p(3)],'local state estimate $\hat{x}_2$','quantization range $[C_2-L_2,C_2+L_2]$','deocoded state $P_{d,2}$','FontSize',15,'Interpreter','latex')
xlim([0.07,1.01])
hold off


%%
err_max = Lk/N;
err_q = abs(xs-Cd);
L_q_sim = err_max(:,end);
% L_q_sim = vpa(L_q_sim,3);
figure()
ne = 101;
xplot = 0:ne-1;
subplot(2,1,1)
hold on
title('Quantization error at k-th tranmission instant','FontSize',20)
ylabel('Error','FontSize',20)
semilogy(xplot,(err_max(1,1:ne)),'-*b')
semilogy(xplot,(err_q(1,1+1:ne+1)),'-*r')
legend('$\bar{e}_{q,1}$','$e_{q,1}$','Interpreter','latex','FontSize',20)
hold off

text(max(xlim)-20, 0.1, sprintf('$L_{q,1}= 0.0684 $'), 'Horiz','left', 'Vert','bottom','FontSize',20,'Interpreter','latex')

subplot(2,1,2)
hold on
ylabel('Error','FontSize',20)
xlabel('Transmission instant, k','FontSize',20)
semilogy(xplot,(err_max(2,1:ne)),'-*b')
semilogy(xplot,(err_q(2,1+1:ne+1)),'-*r')
legend('$\bar{e}_{q,2}$','$e_{q,2}$','Interpreter','latex','FontSize',20)
hold off
text(max(xlim)-20, 0.1, sprintf('$L_{q,2}= 0.0684 $'), 'Horiz','left', 'Vert','bottom','FontSize',20,'Interpreter','latex')





%% testing for error convergance
(eye(n)-LamHAT./N)^(-1)*[1;1]*Beta /4
