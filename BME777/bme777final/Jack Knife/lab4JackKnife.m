% clear;
% close all;
% load DataLab2_1.mat
% Data = DataLab2_1;
% ClassSplit = 50;
% DataSplitRate = 0.4;
% InitialParameterSet = [0 0 1];
% LearningRate = 0.01;
% Theta = 0;
% MaxNoOfIteration = 300;

function [TrainedParameterSet,NoOfIteration] =lab2(Data,ClassSplit,DataSplitRate,InitialParameterSet,LearningRate,Theta,MaxNoOfIteration)
close all;

[Len,~] = size(Data);
Class1 =[];
Class2 =[];

i=1; %count
k=1; %count for # of data points in class 2
while i<=Len
    if Data(i,3)==1 %class 1
        Class1(i,1)=Data(i,1); %Class 1 feat 1
        Class1(i,2)=Data(i,2); %Class 1 feat 2
    else %class 2
        Class2(k,1)=Data(i,1); %Class 2 feat 1
        Class2(k,2)=Data(i,2); %Class 2 feat 2
        k=k+1;
    end
    i=i+1;
end

% Calculate the number of training samples.
[ClassSize1,~]=size(Class1);
[ClassSize2,~]=size(Class2);
Train_Num1=ClassSize1;
Train_Num2=ClassSize2;

RowToTest=1; %row 1 <- "Test" is the row # we will use to jack knife/leave out

Accuracy=zeros(length(Data),1);

while (RowToTest<=length(Data)) %iterate through all rows of Data

    % Split the data in class 1 into training and testing sets. 
    Train_Class1=[];
    i=1;
    while i<=ClassSize1
        Train_Class1(i,:)=Class1(i,:);
        i=i+1;
    end
    
    % Split the data in class 2 into training and testing sets.
    Train_Class2=[];
    i=1;
    while i<=ClassSize2
        Train_Class2(i,:)=Class2(i,:);
        i=i+1;
    end

    % Append all training and test sets
    Train_Class1= [ones(Train_Num1,1) Train_Class1];
    Train_Class2= (-1).*[ones(Train_Num2,1) Train_Class2]; %Normalized Training Class 2

    % Prepare the training data including all training samples of classs 1 and 2.
    Tot_Train_Num=Train_Num1+Train_Num2;
    Train_Data=zeros(Tot_Train_Num,3);
    
    % column/Feature 1
    Train_Data(1:Train_Num1, 1) = Train_Class1(:,1); 
    Train_Data(Train_Num1+1:Tot_Train_Num, 1) = Train_Class2(:,1);  % rows 21 to 50
    
    % column/Feature 2
    Train_Data(1:Train_Num1, 2:3)=Train_Class1(:,2:3);
    Train_Data(Train_Num1+1:Tot_Train_Num, 2:3) = Train_Class2(:,2:3);

    %Row to be jack knife'd / the row that was left out
    Test_Data(1,:)=Train_Data(RowToTest,:);

    % Update Training data set to remove the jack knife'd row
    Train_Data(RowToTest,:)=[]; %nullified

    % Implement basic gradient algorithm.
    OptParams = InitialParameterSet;
    PerceptronFunction = zeros(MaxNoOfIteration,1);
    Criterion = 1;
    NoOfIteration = 1;
    
    while ((Criterion>Theta))
        % Update the PerceptionFunction and The GradientOfCost.
            PredictedValue=OptParams*Train_Data';
            k=find(PredictedValue<0,(Train_Num1+Train_Num2));
            GradientOfCost=sum(Train_Data(k,:));      
        % Update the optimized parameters.
            OptParams = OptParams + LearningRate*GradientOfCost;
            PerceptronFunction(NoOfIteration)=sum(OptParams.*GradientOfCost);
        % Update the value of the criterion to stop the algorithm.
            Criterion= sqrt(GradientOfCost(1).^2+GradientOfCost(2).^2+GradientOfCost(3).^2)*LearningRate;
        %Break the algorithm when the NoOfIteration = MaxNoOfIteration.
            if(NoOfIteration == MaxNoOfIteration)
                break;
            end
            NoOfIteration = NoOfIteration + 1;   
    end

    % Calculate the classification accuracy of the predictions on the test data.
    PredictedValueTest=OptParams*Test_Data';
    NoOfAccuracy=find(PredictedValueTest<0,(length(Test_Data)));

    ClassificationAccuracy =((1-length(NoOfAccuracy))/1*100);
    Accuracy(RowToTest)=ClassificationAccuracy;

    RowToTest=RowToTest+1;
end

FinalAccuracy=mean(Accuracy);
text=['Jack Knife Accuracy: ', num2str(FinalAccuracy),'%'];
disp(text);