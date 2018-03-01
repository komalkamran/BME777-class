Eta = 0.1;
Theta = 0.001;
MaxNoOfIteration = 300;
Problem = 0;

load DataLab3.mat;
Data= DataLab3;

[J,w] = Lab3Data(Eta,Theta,MaxNoOfIteration,Problem,Data);


 