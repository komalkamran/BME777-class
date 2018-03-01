load('DataLab2_1.mat');
Data = DataLab2_1;
ClassSplit = 50;
DataSplitRate = 0.4;
InitialParameterSet = [0 0 1];
LearningRate = 0.01;
Theta = 0;
MaxNoOfIteration = 300;

lab4JackKnife(Data,ClassSplit,DataSplitRate, InitialParameterSet,LearningRate,Theta,MaxNoOfIteration);
