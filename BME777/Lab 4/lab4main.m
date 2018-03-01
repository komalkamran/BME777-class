load DataLab4
Data=Breast_Tissue;
InitMean1=[300 0.7 5800]';
InitMean2=[350 0.3 600]';
InitMean3=[400 1.1 1082]';
MaxNoOfIteration=400;

[PrevMeant,Label]=lab4(Data,InitMean1,InitMean2,InitMean3, MaxNoOfIteration);
PrevMeant

figure;
scatter3(Data(:,1), Data(:, 2), Data(:, 3), '*');
hold on
scatter3(PrevMeant(:,1), PrevMeant(:,2), PrevMeant(:,3), 'k');
xlabel('X1'); ylabel('X2'); zlabel('X3');
title('Scatter plot of Data with final mu');