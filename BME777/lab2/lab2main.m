% load DataLab2_1.mat
Data = DataLab2_1;
ClassSplit = 50;
DataSplitRate = 0.5;
InitialParameterSet = [0 0 1];
LearningRate = 0.01;
Theta = 0.8;
MaxNoOfIteration = 300;
[OptimizedParameterSet,NoOfIteration] = lab2(Data,ClassSplit,DataSplitRate, InitialParameterSet,LearningRate,Theta,MaxNoOfIteration);
disp(OptimizedParameterSet);
disp(NoOfIteration);