function [BootstrapAccuracy] =lab2(Data,ClassSplit,DataSplitRate,InitialParameterSet,LearningRate,Theta,MaxNoOfIteration,NumberOfBootStraps)
close all;

[Len,~] = size(Data);

i=1; j=1;
while i<=Len
    if Data(i,3)==1 % divide data into Class 1 portion
        Class1(i,1)=Data(i,1);
        Class1(i,2)=Data(i,2);
    else % divide the rest of the data into Class 2 portion
        Class2(j,1)=Data(i,1);
        Class2(j,2)=Data(i,2);
        j=j+1;
    end
    i=i+1;
end

% Calculate the number of training samples.
[ClassSize1,~]=size(Class1);
[ClassSize2,~]=size(Class2);
Train_Num1=ClassSize1;
Train_Num2=ClassSize2;

% Define a training set for Class 1
for i=1:ClassSize1
    Train_Class1(i,:)=Class1(i,:);
end

% Define a training set for Class 2
for i=1:ClassSize2
    Train_Class2(i,:)=Class2(i,:);
end

% Append all training sets
Train_Class1= [ones(Train_Num1,1) Train_Class1];
Train_Class2= (-1).*[ones(Train_Num2,1) Train_Class2]; %Normalized Training Class 2

% Prepare the training data including all training samples of classs 1 and 2.
Tot_TrainNum = Train_Num1 + Train_Num2;
CompiledData=zeros(Tot_TrainNum,3);

%column/Feature 1
CompiledData(1:Train_Num1, 1) = Train_Class1(:,1); %rows 1-20
CompiledData(Train_Num1+1:Tot_TrainNum, 1) = Train_Class2(:,1);  %rows 21-50

%column/Feature 2 & 3
CompiledData(1:Train_Num1, 2:3) = Train_Class1(:, 2:3);
CompiledData(Train_Num1+1:Tot_TrainNum, 2:3) = Train_Class2(:, 2:3);

for iteration=1:NumberOfBootStraps

    randnums=randi(100,1,100); %returns an 1-by-100 matrix of psuedorandom #s from 1-100
    BStrap_Train_num=randnums(1:ClassSize1);
    BStrap_Test_num=randnums(ClassSize1+1: ClassSize1+ClassSize2);

    %Update test/train sets for Class 1
    for i=1:ClassSize1
        Train_Data1(i,:)=CompiledData(BStrap_Train_num(i),:);
        Test_Data1(i,:)=CompiledData(BStrap_Test_num(i),:);
    end
    
    %Update test/train sets for Class 2
    for i=1:ClassSize2
        Train_Data2(i,:)=CompiledData(BStrap_Test_num(i),:);
        Test_Data2(i,:)=CompiledData(BStrap_Train_num(i),:);
    end

    % Implement basic gradient algorithm.
    OptParams = InitialParameterSet;
    PerceptronFunction = zeros(MaxNoOfIteration,1);
    Criterion = 1;
    NoOfIteration = 1;
    while ((Criterion>Theta))
        % Update the PerceptionFunction and The GradientOfCost.
            PredictedValue=OptParams*Train_Data1';
            intstep=find(PredictedValue<0,(Tot_TrainNum));
            GradientOfCost=sum(Train_Data1(intstep,:));
        % Update the optimized parameters.
            OptParams = OptParams + LearningRate*GradientOfCost;
            PerceptronFunction(NoOfIteration)=sum(OptParams.*GradientOfCost);
        % Update the value of the criterion to stop the algorithm.
            % Criterion= sqrt(GradientOfCost(1).^2+GradientOfCost(2).^2+GradientOfCost(3).^2)*LearningRate;
        for i=1:(Tot_TrainNum) 
            if (PredictedValue>0)
                Criterion =0;
            end
        end
        
        %Break the algorithm when the NoOfIteration = MaxNoOfIteration.
            if(NoOfIteration == MaxNoOfIteration)
                break;
            end
            NoOfIteration = NoOfIteration + 1;   
    end

    % % Calculate the classification accuracy of the predictions on the test data.
    PredictedValueTest=OptParams*Test_Data1';
    NoOfAccuracy=find(PredictedValueTest<0,(length(Test_Data1)));

    ClassificationAccuracy1 =((length(Test_Data1)-length(NoOfAccuracy))/length(Test_Data1)*100);


    % Repeat Implement basic gradient algorithm.
    OptParams = InitialParameterSet;
    PerceptronFunction = zeros(MaxNoOfIteration,1);
    Criterion = 1;
    NoOfIteration = 1;
    while ((Criterion>Theta))
        % Update the PerceptionFunction and The GradientOfCost.
            PredictedValue=OptParams*Train_Data2';
            intstep=find(PredictedValue<0,(Tot_TrainNum));
            GradientOfCost=sum(Train_Data2(intstep,:));
        % Update the optimized parameters.
            OptParams = OptParams + LearningRate*GradientOfCost;
            PerceptronFunction(NoOfIteration)=sum(OptParams.*GradientOfCost);
        % Update the value of the criterion to stop the algorithm.
            %Criterion= sqrt(GradientOfCost(1).^2+GradientOfCost(2).^2+GradientOfCost(3).^2)*LearningRate;
        for i=1:(Train_Num1 + Train_Num2) 
            if (PredictedValue>0)
                Criterion =0;
            end
        end
        %Break the algorithm when the NoOfIteration = MaxNoOfIteration.
            if(NoOfIteration == MaxNoOfIteration)
                break;
            end
            NoOfIteration = NoOfIteration + 1;   
    end

    % % Calculate the classification accuracy of the predictions on the test data.
    PredictedValueTest=OptParams*Test_Data2';
    NoOfAccuracy=find(PredictedValueTest<0,(length(Test_Data2)));

    ClassificationAccuracy2 =((length(Test_Data2)-length(NoOfAccuracy))/length(Test_Data2)*100);

    FinalClassificationAccuracy(iteration)=(ClassificationAccuracy1+ClassificationAccuracy2)/2;

end

BootstrapAccuracy=mean(FinalClassificationAccuracy);
% text=['Bootstrap Accuracy = ', num2str(BootstrapAccuracy), '%'];
% disp(text);