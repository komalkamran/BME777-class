load('Diabetes.mat');
Data = Diabetes;
FeatureForClassification = 1;
LabelColumn = 3;
FeatureX = 100;

[PosProb, G] = lab1(FeatureX, Diabetes,FeatureForClassification, LabelColumn);

[difference, index_At_G_Equals_0] = min(abs(G));
x0 = FeatureX(index_At_G_Equals_0);


fprintf('The Optimal Threshold is approximately at: %.2f\n', x0);    
% fprintf('The G value is\n', G); 