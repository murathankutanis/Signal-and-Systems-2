%%%%%%%%%%%%%%%%%%%%%%% EE305 HW Q2 %%%%%%%%%%%%%%%%%%%%%%%

clear 
%% Part a
t=0:0.01:10;
figure
hold on
for i=1:8
    theta=2*pi/8*(i-1);
    w=i;
    x=2*cos(w*t+theta);
    legendString = ['x(' num2str(i) ') = 2*cos(' num2str(w) 't + ' num2str(theta) ')'];
    plot(t,x,'DisplayName',legendString)
    xlabel('t')
    ylabel('x(t)')
    grid
    legend
end

%% Part b
figure;
set(gcf,'Position',[100 100 1000 1000]);
t1=0.8;
no_bins=100;
MCnum=10*10^3;
w=10*rand(1,MCnum);
theta_values=2*pi*rand(1,MCnum);
x_values=2*cos(w*t1+theta_values);
histogram(x_values,no_bins)
hold on
min_x=min(x_values);
 max_x=max(x_values);
A=max(abs(min_x),max_x);
x1=linspace(-A,A,1000);
plot(x1,1./(pi*sqrt(A^2-x1.^2))*MCnum/(no_bins/(2*A)))
title(['t_1=',num2str(t1)])
xlabel('x_1')
hold off

% The histogram matches the theoretical density, it validates that the 
% process behaves as expected for a random phase cosine, indicating that 
% the first-order density is accurately represented in the simulation.

%% Part c
figure;
set(gcf,'Position',[100 100 1000 1000]);
count=1;
MCnum=50*10^3;
w=10*rand(1,MCnum);
theta_values=2*pi*rand(1,MCnum);
for Delta=[0 0.95 1.97 10.84]
    t1=0.8+Delta; t2=1.23+Delta;
    no_bins=50;
    x1_values=2*cos(w*t1+theta_values);
    x2_values=2*cos(w*t2+theta_values);
    subplot(2,2,count)
    hist3([x1_values',x2_values'],'Nbins',[no_bins no_bins])
    title(['t_1=',num2str(t1),' t_2=',num2str(t2),'Delta=',num2str(Delta)],'Fontsize',14)
    xlabel('x_1')
    ylabel('y_1')
    count=count+1;
end

% The histograms look similar for all the different Δ values, it indicates 
% that the second-order density does not change significantly with 
% different time intervals, supporting the idea of second-order stationarity.