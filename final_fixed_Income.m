close all
clear all
clc
format short
rng(123456789)
%**************************************************************************
%% Data manipulation
%**************************************************************************
Dates= readtable(('Euribor.xlsx'),Sheet="Sheet2"); 
df1 = table2array(readtable("data.txt"));
data = (flip(df1))/100;
%**************************************************************************
%% Vasicek 1
%**************************************************************************
a = -0.0013; mu = -0.0070;sigma= 1.551742635401969e-04 ; r0 =0.0390 ;T = 3; % VaR Horizon set to 3 for 3-years Horizon
nsteps = 12*21; dt = T/nsteps% days in a month * number of months in a year  
M = 1E3; 
[ri1 rall1] = VasicekModelCW(a , mu, sigma,r0,T,nsteps,M);

timesteps=[0:nsteps]*dt;
Theoretical_Vasiceck_mean = mu + exp(-a*timesteps)*(r0-mu);
Simulated_mean = mean(rall1')';
Theoretical_Vasicek_Variance = (sigma^2/(2*a)) *  (1-exp(-2*a*timesteps));
Simulated_Variance = var(rall1')';

figure
plot(rall1)
hold on 
plot(Theoretical_Vasiceck_mean)
xlim([0 nsteps])
xlabel('Time','Interpreter','latex')
ylabel('r(i)','Interpreter','latex')
legend('Simulated paths','Expected path','Interpreter','latex')
title('Simulated paths of the Vasicek model with a = -0.0013, mu = -0.0070 , sigma =0.000157','Interpreter','latex')

figure
hold on
plot(ri1,'o','MarkerFaceColor',[0 0.4470 0.7410]) 
plot(ri1,'-','LineWidth',1.5,'Color',[0 0.4470 0.7410]) 
xlabel('Time','Interpreter','latex')
ylabel('r(i)','Interpreter','latex')
title('Istantaneous Rate','Interpreter','latex')
grid minor

fixed_maturity = linspace(1,10,1000); 
P_t1 = zeros(nsteps+1, length(fixed_maturity));
R_t1=zeros(nsteps+1, length(fixed_maturity));
R_t_continuously1 = zeros(nsteps+1, length(fixed_maturity));
R_t_proxy1 = zeros(nsteps+1, length(fixed_maturity));
L_t1 = zeros(nsteps+1, length(fixed_maturity));
Vts1 = (sigma^2*(1-exp(-2*a*dt)))/(2*a);

for i = 1:nsteps+1 
    for k = 1:length(fixed_maturity)
        rij = ri1(i);
        maturity = fixed_maturity(k);
        B_t = (1-exp(-a*maturity))/(a);
        A_t = (B_t-maturity)* (a^2*mu-(sigma^2)/2)/(a^2) - ((sigma^2)/(a*4))* B_t^2;
        P_t1(i,k) = exp(-B_t*rij+A_t);
        R_t1(i,k) = -log(P_t1(i,k))/maturity;
        R_t_continuously1(i,k) = (B_t/maturity)*rij + A_t/maturity;
        R_t_Proxy1 = ((B_t/maturity) * rij + (A_t / maturity)) * exp(-a*dt);
        L_t_Proxy1 = ((B_t/maturity)*mu - A_t/maturity)*(1-exp(-a*dt));
        R_t_Proxy2 = sigma * (B_t / maturity) * (randn(1,1)*((sigma^2)/(2*a))+mu);
        R_t_proxy1(i,k) = R_t_Proxy1 +L_t_Proxy1 + R_t_Proxy2;
        accrual = (max(fixed_maturity))-maturity;
        L_t1(i,k) = 1/(accrual) * ((1-P_t1(i,k))/P_t1(i,k));
    end
end
%additional_quantitites
I_t1 = ri1(1:end-1)+ri1(2:end) *(dt./2);
MMA1 = (exp(I_t1(1:end)));

figure %Add subplots
plot(fixed_maturity,P_t1(1:10,:))
xlabel('Time', 'interpreter', 'latex'); 
ylabel('R(t)','interpreter','latex');
legend('Simulated DIscount Factors','interpreter','latex');
title('Spot Rate R(t) for different fixed maturities', 'interpreter','latex');

%**************************************************************************
%% Vasicek Estimated from real World Calibration
%**************************************************************************
X = [ones(size(data,1)-1,1) data(1:end-1)];Y = [(data(2:end,1))];

[coefficients, intervals, residuals] = regress(Y, X);alpha_hat = -log(coefficients(2,1))/dt;mu_hat = coefficients(1)/(1-coefficients(2));
sigma_hat = std(residuals)/ (sqrt((1-exp(-2*alpha_hat*dt))/(2*alpha_hat)));
r0=data(end,1);T =3;nsteps = 12*21;M = 1E3;
figure('Color',[1 1 1])
scatter(Y,X(:,2),'.')
lsline()
fitlm(X(:,2),Y,"linear")
xlabel('Euribor t','interpreter','latex')
ylabel('Euribor t+1','interpreter','latex')
title('OLS Estimation of VSCK parameters','interpreter','latex')
first_values = [ coefficients(2,1) ...
mu_hat alpha_hat sigma_hat  ];
titles = {'Beta Hat', 'Mu Hat', 'Alpha Hat', 'Sigma Hat'};
initial_param_estimates = table(first_values(1), first_values(2), first_values(3), first_values(4), ...
    'VariableNames', titles)
%%
[ri2 rall2]= VasicekModelCW(alpha_hat, mu_hat, sigma_hat,r0,T, nsteps,M );

P_t2 = zeros(nsteps+1, length(fixed_maturity));
R_t2=zeros(nsteps+1, length(fixed_maturity));
R_t_continuously2 = zeros(nsteps+1, length(fixed_maturity));
R_t_proxy23 = zeros(nsteps+1, length(fixed_maturity));
L_t2 = zeros(nsteps+1, length(fixed_maturity));
Vts2 = (sigma_hat^2*(1-exp(-2*alpha_hat*dt)))/(2*alpha_hat);

for i = 1:nsteps+1 
    for k = 1:length(fixed_maturity)
        rij = ri2(i);
        maturity = fixed_maturity(k);
        B_t = (1-exp(-alpha_hat*maturity))/(alpha_hat);
        A_t = (B_t-maturity)* (alpha_hat^2*mu_hat-(sigma_hat^2)/2)/(alpha_hat^2) - ((sigma_hat^2)/(alpha_hat*4))* B_t^2;
        P_t2(i,k) = exp(-B_t*rij+A_t);
        R_t2(i,k) = -log(P_t2(i,k))/maturity;
        R_t_continuously2(i,k) = (B_t/maturity)*rij + A_t/maturity;
        R_t_Proxy1 = ((B_t/maturity) * rij + (A_t / maturity)) * exp(-alpha_hat*dt);
        L_t_Proxy1 = ((B_t/maturity)*mu_hat - A_t/maturity)*(1-exp(-alpha_hat*dt));
        L_t_Proxyregr(k) = ((B_t/maturity)*mu_hat - A_t/maturity);
        R_t_Proxy2 = sigma_hat * (B_t / maturity) * (randn(1,1)*((sigma_hat^2)/(2*alpha_hat))+mu_hat);
        R_t_proxy23(i,k) = R_t_Proxy1 +L_t_Proxy1 + R_t_Proxy2;
        accrual = (max(fixed_maturity))-maturity;
        L_t2(i,k) = 1/(accrual) * ((1-P_t2(i,k))/P_t2(i,k));
    end
end
%additional_quantitites
I_t1 = ri2(1:end-1)+ri2(2:end) *(dt./2);
MMA1 = (exp(I_t1(1:end)));
delta = 1;

figure
plot(R_t_proxy23(:,2))
%%
Theoretical_Vasiceck_mean = mu_hat + exp(-alpha_hat*timesteps)*(r0-mu_hat);
Simulated_mean = mean(rall2')';
Theoretical_Vasicek_Variance = (sigma_hat^2/(2*alpha_hat)) *  (1-exp(-2*alpha_hat*timesteps));
Simulated_Variance = var(rall2')';

figure('Color',[1 1 1])
plot(rall2)
hold on 
plot(Theoretical_Vasiceck_mean)
xlim([0 nsteps])
xlabel('Time','Interpreter','latex')
legend('Simulated paths','Expected path','Interpreter','latex')
title('Simulated paths of the Vasicek model Calibrated on Market Data','Interpreter','latex')

figure('COlor',[1 1 1])
hold on
%plot(ri2,'o','MarkerFaceColor',[0 0.4470 0.7410]) 
plot(ri2,'-','LineWidth',1.5,'Color',[0 0.4470 0.7410]) 
xlabel('Time','Interpreter','latex')
ylabel('r(i)','Interpreter','latex')
title('Istantaneous Rate','Interpreter','latex')
xlim([0 nsteps])

%%


figure('COlor',[1 1 1])
plot(P_t2([2 4 10 100 150 200 250],1))
xlabel('Time','Interpreter','latex')
ylabel('P(t)','Interpreter','latex')
title('Discount Factor Calibrated VSCK','Interpreter','latex')
%xlim([0 nsteps])





%%
Y = [diff(R_t_proxy23(1,:),delta)]';
X = [ones(size(Y,1),1) (R_t_proxy23(1,delta:end-1))'];

[cal_coef, intervals, residuals_cal] = regress(Y, X);
alpha = cal_coef(1,1)*100;
b_real = exp(-alpha*delta); a_real = b_real * L_t_Proxyregr(delta,delta);
sigma_real = sqrt(Vts2)*(B_t/maturity)*sqrt(abs(1-b_real))*2*alpha;
%sigma_real = std(residuals_cal)/ (sqrt((1-exp(-2*cal_coef(1,1)*dt))/(2*cal_coef(1,1))));

%[b_real a_real sigma_real] = recoverParameters(cal_coef(1),1,B_t,A_t,1,mu_hat,sigma_hat)


[ri3 rall3]= VasicekModelCW(a_real, b_real, sigma_real,r0,T, nsteps,M );


P_t3 = zeros(nsteps+1, length(fixed_maturity));
R_t3=zeros(nsteps+1, length(fixed_maturity));
R_t_continuously3 = zeros(nsteps+1, length(fixed_maturity));
R_t_proxy3 = zeros(nsteps+1, length(fixed_maturity));
L_t3 = zeros(nsteps+1, length(fixed_maturity));
Vts3 = (sigma_real^2*(1-exp(-2*a_real*dt)))/(2*a_real);

for i = 1:nsteps+1 
    for k = 1:length(fixed_maturity)
        rij = ri3(i);
        maturity = fixed_maturity(k);
        B_t = (1-exp(-a_real*maturity))/(a_real);
        A_t = (B_t-maturity)* (a_real^2*b_real-(sigma_real^2)/2)/(a_real^2) - ((sigma_real^2)/(a_real*4))* B_t^2;
        P_t3(i,k) = exp(-B_t*rij+A_t);
        R_t3(i,k) = -log(P_t3(i,k))/maturity;
        R_t_continuously3(i,k) = (B_t/maturity)*rij + A_t/maturity;
        R_t_Proxy1 = ((B_t/maturity) * rij + (A_t / maturity)) * exp(-a_real*dt);
        L_t_Proxy1 = ((B_t/maturity)*b_real - A_t/maturity)*(1-exp(-a_real*dt));
        L_t_Proxyregr3(k) = ((B_t/maturity)*b_real - A_t/maturity);
        R_t_Proxy2 = sigma_real * (B_t / maturity) * (randn(1,1)*((sigma_real^2)/(2*a_real))+b_real);
        R_t_proxy3(i,k) = R_t_Proxy1 +L_t_Proxy1 + R_t_Proxy2;
        accrual = (max(fixed_maturity))-maturity;
        L_t3(i,k) = 1/(accrual) * ((1-P_t3(i,k))/P_t3(i,k));
    end
end

figure('Color',[1 1 1])
plot(R_t_proxy23(1,:))
hold on
plot((abs(R_t_proxy3(1,:))))
xlabel('Maturity','Interpreter','latex')
ylabel('R(t)','Interpreter','latex')
title('Proxy Values of r(t) comparison over T','Interpreter','latex')
legend('R(t)1 Proxy','Calibrated Proxy','interpreter','latex')


%**************************************************************************
%% Vasicek Estimated Risk Neutral Calibration to Real World
%**************************************************************************

[ri4 rall4] = VasicekModelCW(0.0625068816710969,3.0487211332473/100,0.000246489412412363,1.53436155851349/100,3,nsteps,M);

spot = (readmatrix("Spot.txt"))/100;
spot_mat = (readmatrix("right_SPOT.txt"))/100;

P_t4 = zeros(nsteps+1, length(fixed_maturity));
R_t4=zeros(nsteps+1, length(fixed_maturity));
R_t_continuously4 = zeros(nsteps+1, length(fixed_maturity));
for i = 1:nsteps+1 
    for k = 1:length(fixed_maturity)
        rij = ri4(i);
        maturity = fixed_maturity(k);
        B_t = (1-exp(-alpha_hat*maturity))/(alpha_hat);
        A_t = (B_t-maturity)* (alpha_hat^2*mu_hat-(sigma_hat^2)/2)/(alpha_hat^2) - ((sigma_hat^2)/(alpha_hat*4))* B_t^2;
        P_t4(i,k) = exp(-B_t*rij+A_t);
        R_t4(i,k) = -log(P_t4(i,k))/maturity;
        R_t_continuously4(i,k) = (B_t/maturity)*rij + A_t/maturity;

    end
end

figure
plot(P_t4(1,:))
xlabel('Time', 'interpreter', 'latex'); 
ylabel('R(t)','interpreter','latex');
legend('Simulated DIscount Factors','interpreter','latex');
title('Spot Rate R(t) for different fixed maturities', 'interpreter','latex');

%%
S = spot(end-size(ri4,1)+1:end,1);
for n = length(S) - 1;
Sx = sum(S(1:end-1));
Sy = sum(S(2:end));
Sxx = sum(S(1:end-1).^2);
Sxy = sum(S(1:end-1).*S(2:end));
Syy = sum(S(2:end).^2);
mu_hat3 = (Sy*Sxx - Sx*Sxy) / (n*(Sxx - Sxy) - (Sx^2 - Sx*Sy));
alpha_hat3 = ((Sxy - mu_hat3*Sx - mu_hat3*Sy + n*mu_hat3^2) / (Sxx - 2*mu_hat3*Sx + n*mu_hat3^2)) / dt;
a = 1 - alpha_hat3 * dt;
sigmah2 = (Syy - 2*a*Sxy + a^2*Sxx - 2*mu_hat3*(1-a)*(Sy - a*Sx) + n*mu_hat3^2*(1-a)^2) / n;
sigma_hat3 = sqrt(sigmah2 * 2 * alpha_hat3 / (1 - a^2));
end
%%
[parameters, expectedinfomatrix, lvalue] = recoverParameter(S,6);
parameters
%%
a = alpha_hat3; mu = mu_hat3;sigma= sigma_hat3 ; r0 = spot(1);T = 3;
MLE = [alpha_hat3 mu_hat3 sigma_hat3];
%a = parameters(1); mu = parameters(2);sigma= parameters(3) ; r0 = spot(1);T = 3; % VaR Horizon
nsteps = 21*12; % days in a year  
M = 1E3; dt = T/nsteps;
%price_function = vasicekPrice(a,mu,sigma,r0,3,5,1,0.2)
Vts = (sigma^2*(1-exp(-2*a*dt)))/(2*a);
rall = [];
for j = 1:M
	r = zeros(nsteps+1,1); r(1) = r0;
	Z = randn(1,nsteps);%To ask
	for i = 1:nsteps
        Its = (r(i)*exp(-a*dt)) + (mu*(1-exp(-a*dt)));
        r(i+1) = Its + (sqrt(Vts)*Z(i));
    end
	rall = [rall, r];
end
ri3 = [rall(:,end)] ;

P_t3 = zeros(nsteps+1, length(fixed_maturity));
R_t3=zeros(nsteps+1, length(fixed_maturity));
R_t_continuously3 = zeros(nsteps+1, length(fixed_maturity));
for i = 1:nsteps+1 
    for k = 1:length(fixed_maturity)
        rij = ri3(i);
        maturity = fixed_maturity(k);
        B_t = ((1-exp(-alpha_hat3*maturity))/(alpha_hat3));
        A_t = ((B_t-maturity)* (alpha_hat3^2*mu_hat3-(sigma_hat3^2)/2)/(alpha_hat3^2) - ((sigma_hat3^2)/(alpha_hat3*4))* B_t^2);
        P_t3(i,k) = exp(-B_t*rij+A_t);
        R_t3(i,k) = -log(P_t3(i,k))/maturity;
        R_t_continuously3(i,k) = (B_t/maturity)*rij + A_t/maturity;
    end
end

figure
plot(R_t3(1,:))
xlim([1 size(P_t3,2)])
xlabel('Time to maturity','Interpreter','latex')
ylabel('P(t)3','Interpreter','latex')
title('Simulated ZCB Price')
%%
figure
plot((ri3(2:end,2)))
hold on
yyaxis right
plot(spot(end-250:end,1))
%%
frenchtesco = 6:6:60;
figure('Color',[1 1 1])
plot(R_t3(1,frenchtesco))
hold on 
plot(spot_mat)
plot(R_t2(1,frenchtesco))
%plot(R_t4(1,:))
%plot(R_t3(end,:))
xlabel('Time to Maturities','Interpreter','latex')
ylabel('R(t)','Interpreter','latex')
legend('R(t)1 Initial','Observed R(t)','R(t)2 Best Fit','Interpreter','latex')
title('Comparison','Interpreter','latex')

%**************************************************************************
%% Bond Pricing
%**************************************************************************
dt = 3/size(frenchtesco,1)
MMA = ri3(1:end-1)+ri3(2:end) .*(dt./2);
MMA = (exp(MMA(1:end)));
frenchtesco = 6:6:36;
%Pts = rall(frenchtesco,:)
Pts = P_t3(frenchtesco,:);
c = (4/100 *1000)/2;
c_s = repmat(c,size(Pts,1),1);
discount_factor = mean(Pts)';
FV = 1000;
bond_exact = sum(c_s./MMA(frenchtesco))+ FV/MMA(36)
%**************************************************************************
%% Bond Pricing at the VaR horizon and P&L disttibution
%**************************************************************************
ci = 0.20;
for i = 1:size(P_t3,1)
    for j = 1:size(fixed_maturity,2)
        BP_t(i,j) = (ci * P_t3(i,j));
    end
end
BP_t_VaR1 = sum(BP_t(1:35,:)*0.5,2)+ P_t3(36,:)*FV
BP_t_VaR = sum(BP_t(1:35,:)*0.5)+ P_t3(36,:)*FV

VaR_RN = ((BP_t_VaR(1,1)*sqrt(31)*sigma_hat*norminv(0.01))+mu_hat*sqrt(31))
ES = ((-BP_t_VaR(1,1)*sqrt(31)*Vts * normpdf(norminv(0.01)))/(0.01))*100
VaR_exact_violations = BP_t_VaR(BP_t_VaR<VaR_RN);
[VaR_RN ES]
%%
%**************************************************************************
%% Bond Pricing at the VaR horizon and P&L disttibution with MMA
%**************************************************************************
%%

for i = 1:size(ri4,1)

    term_structure_from(i) = get_df_Vasicek([parameters  ri4(i)],3)
end
figure
plot(term_structure_from+0.02880)
hold on
plot(P_t4(:,1))
%%

accrual = (max(frenchtesco)- frenchtesco)/252;

MMAs = rall3(1:end-1,:) - rall3(2:end,:).*(dt./2);
MMAs = (exp(MMAs));
bond_exact_sim= sum(c./MMAs(frenchtesco,:)) + FV./MMAs(36,:);

bond_exact_sim_RN = (sum(c*P_t3(frenchtesco,2))+(FV*P_t3(36,1)));

scarto = (bond_exact - BP_t_VaR1);





%**************************************************************************
%% Bond Pricing at the VaR horizon and P&L disttibution with MMA
%**************************************************************************
MMAs = rall(1:end-1,:) - rall(2:end,:).*(dt./2);
MMAs = (exp(MMAs));
bond_exact_sim= sum(c./MMAs(frenchtesco,:)) + FV./MMAs(36,:);

bond_exact_sim_RN = (sum(c*P_t3,2)+(FV*P_t2(:,end)))

scarto = (bond_exact - bond_exact_sim)

figure('Color',[1 1 1])
histfit(scarto)
ylabel('Frequncies','Interpreter','latex')
title('P&L distribution','Interpreter' , 'latex')
[h,p,jbstat,critval] = jbtest(scarto)
jbtest(scarto,0.001)

%% NoN Parametric VaR NoN Parametric EXpected Shortfall
iter = sort(scarto,'ascend');
VaR1 = prctile(iter,(1-0.995)*100); 
Var_Non_Parametric_in_returns = VaR1;

figure('Color',[1 1 1])
plot(scarto)
hold on
line([0,length(scarto)], [Var_Non_Parametric_in_returns, Var_Non_Parametric_in_returns],'Color','r',  'LineWidth', 2)
line([0,length(scarto)], [VaR_RN, VaR_RN], 'LineWidth', 2,'Color','g')
legend('P & L','VaR Non Parametric','VaR Analytical','Interpreter','latex')
%% Analytical VaR Analytical EXpected Shortfall
lossess_exceeding = mean(scarto(scarto< prctile(scarto,(1-0.995)*100)));
Non_Parametric_ES = -(lossess_exceeding-prctile(scarto,(1-0.995)*100))
ES/100
[Non_Parametric_ES lossess_exceeding]
[Var_Non_Parametric_in_returns VaR_RN]
%%
VaR_RN = (scarto(1,1)*sqrt(11)*Vts*norminv(0.01))*100
ES = ((-scarto(1,1)*sqrt(31)*Vts * normpdf(norminv(0.01)))/(0.01))*100
VaR_exact_violations = scarto(scarto<VaR_RN);
[VaR_RN ES]

