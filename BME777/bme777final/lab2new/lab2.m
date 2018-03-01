
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BME777: LAB 2: Linear Discriminant Functions.
% Indian Liver Patient Dataset.
% Link: https://archive.ics.uci.edu/ml/datasets/ILPD+%28Indian+Liver+Patient+Dataset%29#
% Class1: Liver patient. Class2: non Liver patient.
% DataLab2_1: Features: TP Total Proteins and ALB Albumin with modification for problem simplification. 
% Features 8-9. 
% DataLab2_2: Features: TP Total Proteins and A/G Ratio	Albumin and
% Globulin Ratio. Features 8-10.
% 50 samples were extracted for each class.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% 1. Data: 100x3 dataset. The first column contains the feature x1, the second
% column contains the feature x2. The class labels are given in the third
% column.
% 2. ClassSplit: Threshold where classes are divided. See the third
% column of the Data to choose the correct threshold.
% 3. DataSplitRate: Threhold to split the data in each class into training and testing data.
% For e.g., DataSplitRate = 0.4 ==> 40% data of class 1,2 is for training.
% 60% of the data is for testing.
% 4. InitialParameterSet: Initial values of the set of parameters. For
% e.g., InitialParameterSet = [0 0 1].
% 5. LearningRate: Learning rate when updating the algorithm.
% 5. Theta: The expected cost that the optimized parameter set may give.
% 6. MaxNoOfIteration: the maximum number of iterations the algorithm may run.
%
% Output:
% 1: TrainedParameterSet: The set of optimized parameters.
% 2: NoOfIteration: The number of iteration when the algorithm converges.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example of use:
% load DataLab2_1.mat
% Data = DataLab2_1;
% ClassSplit = 50;
% DataSplitRate = 0.4;
% InitialParameterSet = [0 0 1];
% LearningRate = 0.01;
% Theta = 0;
% MaxNoOfIteration = 300;
% [OptimizedParameterSet,NoOfIteration] = ...
% lab2(Data,ClassSplit,DataSplitRate, ... 
% InitialParameterSet,LearningRate,Theta,MaxNoOfIteration);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [TrainedParameterSet,NoOfIteration] =lab2(Data,ClassSplit,DataSplitRate,InitialParameterSet,LearningRate,Theta,MaxNoOfIteration)

close all;

[Len,~] = size(Data);

% Split the data into two classes based on ClassSplit. 
Class1 = Data(1:ClassSplit, 1:2);
Class2 = Data((ClassSplit + 1):end, 1:2);

Class1_Labels = Data(1:ClassSplit, :);
Class2_Labels = Data((ClassSplit + 1):end, :);

% Calculate the number of training samples.
Train_Num1 = DataSplitRate * length(Class1)
Train_Num2 = DataSplitRate * length(Class2)

% Split the data in class 1 into training and testing sets. 
Train_Class1 = Class1(1:Train_Num1, :);
Test_Class1 = Class1(Train_Num1 +1  : end, :);
Test_Class1_Labels = Class1_Labels(Train_Num1 + 1:end, :);

% Split the data in class 2 into training and testing sets.
Train_class2 = Class2(1:Train_Num2, :); % Do not forget to normalize the training data of class 2;
Train_class2 = -1 .* Train_class2;
Test_Class2 =  Class2(Train_Num2 +1  : end, :);

Test_Class2_Labels = Class2_Labels(Train_Num2 + 1:end, :);

% Prepare the training data including all training samples of classs 1 and
% 2.
Train_Data = zeros(Train_Num1 + Train_Num2,3);
Train_Data(1:Train_Num1,1) = ones; 
Train_Data(Train_Num1+1:Train_Num1+Train_Num2,1) = -1 ;
Train_Data(1:Train_Num1,2:3) = Train_Class1;
Train_Data(Train_Num1+1:Train_Num1 + Train_Num2,2:3)= Train_class2;

% Prepare the test data including all test samples of classs 1 and
% 2.
if((Train_Num1+Train_Num2)~=(length(Class1)+length(Class2)))
	Test_Data = zeros(length(Test_Class1) + length(Test_Class2), 3);
	Test_Data(1:length(Test_Class1) + length(Test_Class2),1) = ones;
	Temp = [Test_Class1; Test_Class2];
	Test_Data(:,2:3) = Temp;
    
    
    
end

test_data_labels = [Test_Class1_Labels; Test_Class2_Labels];

% Implement basic gradient algorithm.
OptParams = InitialParameterSet;
PerceptionFunction = zeros(MaxNoOfIteration,1);
Criterion = 1;
NoOfIteration = 1;

while ((Criterion>Theta))

    GradientOfCost = zeros(1,3);
    
    % Update the PerceptionFunction and The GradientOfCost.
    for i=1:(Train_Num1 + Train_Num2)    
        % Use the current OptParams and the ith train data to predict the class.
        PredictedValue = OptParams * Train_Data(i,:)';  
        if (PredictedValue <= 0)
            PerceptionFunction(NoOfIteration) = PerceptionFunction(NoOfIteration)- PredictedValue ;
            GradientOfCost = GradientOfCost + (-1* (Train_Data(i,:))); 
        end
    end
    
    % Update the optimized parameters.
    OptParams = OptParams - LearningRate*GradientOfCost;
    
    % Update the value of the criterion to stop the algorithm.
    if abs(LearningRate * GradientOfCost) < Theta
        
        Criterion = Theta;
    end
    
    % Break the algorithm when the NoOfIteration = MaxNoOfIteration.
    if(NoOfIteration == MaxNoOfIteration)
        break;
    end
    NoOfIteration = NoOfIteration + 1;
    
end

% Plot data of class 1, class 2 and the estimated boundary.
scatter(Class1(:,1),Class1(:,2), 100); 
hold on;
scatter(Class2(:,1),Class2(:,2), 100);
hold off;

%estimatedboundary = (OptParams(1:1) + (OptParams(1:3)*Train_Data(i:);
% Plot the values of the perception function.
figure;
plot(PerceptionFunction);

% Calculate the classification accuracy of the predictions on the test data.
NoOfAccuracy = 0;

if((Train_Num1+Train_Num2)~=(length(Class1)+length(Class2)))
	for j=1:length(Test_Data)
        PredictedValue2 = OptParams * Test_Data(j,:)';
        
        if PredictedValue2 > 0 & test_data_labels(j,3) == 1
            NoOfAccuracy = NoOfAccuracy + 1;
        else if PredictedValue2 < 0 & test_data_labels(j,3) == 2
                NoOfAccuracy = NoOfAccuracy + 1;
             
        % Update the number of correct prediction here. 
            end
            
	end
end
	
ClassificationAccuracy = NoOfAccuracy/length(Test_Data)
disp('Classification Accuracy=')
disp(ClassificationAccuracy);
%NoOfIteration = 12;

TrainedParameterSet = OptParams
disp(TrainedParameterSet);
end