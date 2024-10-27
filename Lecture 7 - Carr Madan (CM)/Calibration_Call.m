clear; close all;
%>European Call Options
Maturity=4/12;
Data=[25.30	200 %Prices, Strikes
10.00	225
18.20	210
41.85	180
15.20	215
7.85	230
12.36	220
46.05	175
70.70	150
50.85	170
60.40	160
29.40	195];
%Data=[Strike,Time to Maturity,Market Price]
Data1=[Data(:,2), Maturity*ones(size(Data(:,2))), Data(:,1)];
Maturity=4/252;
Data=[0.24	230
9.25	210
2.26	220
0.76	225
19.03	200
5.20	215
14.00	205
0.10	235];
Data=[Data(:,2), Maturity*ones(size(Data(:,2))), Data(:,1)];
Data=[Data1;Data];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r=0.02; spot=218.75;
% sigma=params(1);
% lambda=params(2);
% p=params(3);
% lambdap=params(4);
% lambdam=params(5);
x0=[0.5, 5, 0.8, 15, 15];
LB=[0.1, 0, 0, 1,1];
UB=[0.8,20,1 , 50,50];
[params,error]=lsqnonlin(@(params) ...
        fun(params,spot,Data(:,1),r,Data(:,2),Data(:,3)),...
        x0,LB,UB)
figure
plot(Data(:,1),Data(:,3),'+','markersize',16,'linewidth',5); hold on
Price=fun(params,spot,Data(:,1),r,Data(:,2),0);
plot(Data(:,1),Price,'s','markersize',16,'linewidth',5); 
legend('Market price', 'Model price')













