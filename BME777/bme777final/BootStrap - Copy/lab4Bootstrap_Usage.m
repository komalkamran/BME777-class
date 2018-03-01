load DataLab2_1.mat
Data = DataLab2_1;
ClassSplit = 50;
DataSplitRate = 0.4;
InitialParameterSet = [0 0 1];
LearningRate = 0.01;
Theta = 0;
MaxNoOfIteration = 300;
NumberOfBootStraps=20;

for i=1:20
    acc(i)=lab4Bootstrap(Data,ClassSplit,DataSplitRate, InitialParameterSet,LearningRate,Theta,MaxNoOfIteration,NumberOfBootStraps);
end

% max=max(acc);
% min=min(acc);
% 
% scatter(i,max,'r');
% hold on;
% scatter(i,min,'k');
% hold on;
stem(acc);
% legend('Max Accuracy','Min Accuracy');
axis([0 20 0 100]);
xlabel('Iteration');
ylabel('Accuracy (%)');
title('Bootstrap Accuracies');