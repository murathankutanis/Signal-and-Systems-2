clear all
samples=10000;

%% Part A

% As we increase the number of the samples, the probability distribution
%regress around the mean

%Also, the expectaions are not calculated but we can say that, as we 
%increse the n, the mean of the sample become more closer to the 
%expectation of Xi.


VarianceCalculation=[];
i=1;
figure;
set(gcf,'Position', [100 100 1000 1000]);
for NumSums=[1 3 30]
    y=unifrnd(1,4,NumSums,samples);
    SumOfRandomNumbers=sum(y,1);
    subplot(3,2,2*i-1);
    histogram(SumOfRandomNumbers,100,'Normalization','PDF');
    title('PDF of Sn for n=' + string(NumSums));
    xlabel('Value');
    ylabel('f(x)');
    subplot(3,2,2*i);
    cdfplot(SumOfRandomNumbers);
    title('CDF of Sn for n='+ string(NumSums));
    xlabel('Value');
    ylabel('F(x)');
    VarianceCalculation(i,:)=SumOfRandomNumbers;
    i=i+1;

end

%% Part B

XiMeanPartB = (1 + 4) / 2;
disp(['Xi Mean is ' num2str(XiMeanPartB)]);

XiVarPartB = ((4 - 1)^2) / 12;
disp(['Xi Variance is ' num2str(XiVarPartB)]);

count = 1;
for NumSums = [1 3 30]
    SnMeanPartB = sum(VarianceCalculation(count, :)) / samples;
    disp(['Sn Mean for n = ' num2str(NumSums) ' is ' num2str(SnMeanPartB)]);
    
    SnSum = 0;
    for element = 1:length(SumOfRandomNumbers)
        SnSum = SnSum + (VarianceCalculation(count, element) - SnMeanPartB)^2;
    end
    SnVariancePartB = SnSum / length(SumOfRandomNumbers);
    disp(['Sn Variance for n = ' num2str(NumSums) ' is ' num2str(SnVariancePartB)]);
    
    count = count + 1;
end


%% Part C

% As we increase the number of samples, the PDF's become more 
%closer to the Gaussian distribution (From CLT)

figure;
set(gcf, 'Position', [100 100 1000 1000]);
i=1;
for NumSums=[1 3 30]
    y=unifrnd(1,4,NumSums,samples);
    SumOfRandomNumbers=sum(y,1);
    subplot(3,2,2*i-1);
    histogram(SumOfRandomNumbers,100,'Normalization','PDF');
    title('PDF of Sn for n=' + string(NumSums));
    hold on
    SnMeanPartB=mean(SumOfRandomNumbers);
    SnStd=std(SumOfRandomNumbers);
    fplot(@(x) normpdf(x, SnMeanPartB,SnStd), [min(SumOfRandomNumbers), max(SumOfRandomNumbers)], 'r')
    hold off
    xlabel('Value')
    ylabel('F(x)')
    legend('','Gaussian')
    subplot(3,2,2*i);
    cdfplot(SumOfRandomNumbers);
    title('CDF of Sn for n='+ string(NumSums));
    xlabel('Value');
    ylabel('F(x)');
    i=i+1;
end

%% Part D

% In this part, we empirically estimate the mean and variance of the sum of
% random numbers generated from a discrete distribution, and compare the
% distributions with Gaussian distributions.

% Define discrete distribution
Q = [1 2 3 4];
w = [1/6 2/6 1/6 2/6];
w = w / sum(w);


NumSums_values = [1 3 30];
NumSums_length = length(NumSums_values);
y2 = zeros(max(NumSums_values), samples);
figure;
set(gcf, 'Position', [100 100 1000 1000]);

for i = 1:NumSums_length
    NumSums = NumSums_values(i);
    cp = [0, cumsum(w)];
    
    for m = 1:NumSums
        for n = 1:samples
            r = rand;
            ind = find(r > cp, 1, 'last');
            y2(m, n) = Q(ind);
        end
    end
    
    SumOfRandomNumbers = sum(y2(1:NumSums, :), 1);
    
    % Plot histogram
    subplot(NumSums_length, 2, 2*i-1);
    histogram(SumOfRandomNumbers, 100, 'Normalization', 'probability');
    title(['PDF of Sn for n = ' num2str(NumSums)]);
    xlabel('Value');
    ylabel('f(x)');
    

    SnMeanPartD = mean(SumOfRandomNumbers);
    SnStd = std(SumOfRandomNumbers);
    hold on;
    x_values = linspace(min(SumOfRandomNumbers), max(SumOfRandomNumbers), 100);
    plot(x_values, normpdf(x_values, SnMeanPartD, SnStd), 'r', 'LineWidth', 2);
    legend('', 'Gaussian');
    hold off;
    
    % Plot CDF
    subplot(NumSums_length, 2, 2*i);
    cdfplot(SumOfRandomNumbers);
    title(['CDF of Sn for n = ' num2str(NumSums)]);
    xlabel('Value');
    ylabel('F(x)');
    
    % Calculate and display mean and variance
    if i == 1
        XiMeanPartD = sum(Q .* w);
        disp(['Mean of Xi is ' num2str(XiMeanPartD)]);
        XiVarPartD = sum((Q - XiMeanPartD).^2 .* w);
        disp(['Variance of Xi is ' num2str(XiVarPartD)]);
    end
    
    SnMeanPartD = mean(SumOfRandomNumbers);
    SnVariancePartD = var(SumOfRandomNumbers);
    disp(['Mean of Sn for n = ' num2str(NumSums) ' is ' num2str(SnMeanPartD)]);
    disp(['Variance of Sn for n = ' num2str(NumSums) ' is ' num2str(SnVariancePartD)]);
end

