%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BME777: Lab 4:  Unsupervised Learning and Algorithm Independent Machine Learning
% Breast Tissue Dataset: https://archive.ics.uci.edu/ml/datasets/Breast+Tissue
% Feature 1: I0	Impedivity (ohm) at zero frequency.
% Feature 2: HFS high-frequency slope of phase angle.
% Feature 3: AREA area under spectrum.
% Only ther first 3 classes of the orginal dataset are used.
% 14 samples of each class were extracted for clustering. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% 1. Data: Data for clustering.
% 2. InitialMean1: Initial mean vector of class 1.
% 3. InitialMean2: Initial mean vector of class 2.
% 4. InitialMean3: Initial mean vector of class 3.
% 5. MaxNoOfIteration: Maximum number of iteration.
% Outputs:
% 1. FinalMean1: Final mean vector of class 1 given by k-means.
% 2. FinalMean2: Final mean vector of class 2 given by k-means.
% 3. FinalMean3: Final mean vector of class 3 given by k-means.
% 4. Label: Label for each sample in the original data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example of use:
% load DataLab4.mat;
% Data = Breast_Tissue;
% InitMean1 = [choose your mean1 ];
% InitMean2 = [choose your mean2 ];
% InitMean3 = [choose your mean3 ];
% MaxNoOfIteration = 400;
% [FinMean1, FinMean2, FinMean3,Label] = lab4(Data,InitMean1,InitMean2,InitMean3,MaxNoOfIteration);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [PrevMeant,Label]=lab4(Data,InitMean1,InitMean2,InitMean3, MaxNoOfIteration)
% load DataLab4
% Data=Breast_Tissue;
% InitMean1=[400 0.7 5800]';
% InitMean2=[250 0.3 400]';
% InitMean3=[300 1.1 1082]';
% MaxNoOfIteration=400;

%allowing for a 2 or 3 class problem
if InitMean3~=0
   PrevMean = [InitMean1 InitMean2 InitMean3];
else
   PrevMean = [InitMean1 InitMean2];
end
PrevMeant=PrevMean';
Label = zeros(length(Data),1);
Itr = 0;

%% plot
figure;
scatter3(Data(:,1), Data(:, 2), Data(:, 3), '*');
hold on
scatter3(PrevMeant(:,1), PrevMeant(:,2), PrevMeant(:,3), 'R');
xlabel('X1'); ylabel('X2'); zlabel('X3');
title('Scatter plot of Data with mu');

while(1)
    
    Itr = Itr + 1;
    %%%%%%%%%%%%Compute Euclidean distance from each sample to the given means%%%%%%%%%%%%
    for i=1:length(Data)
       %% squared eucledian distance
       %%d1(i)= (x1(i)-mu1(1))^2+(x2(i)-mu1(2))^2+(x3(i)-mu1(3))^2
       %%d2(i)= (x1(i)-mu2(1))^2+(x2(i)-mu2(2))^2+(x3(i)-mu2(3))^2
       
       D1 = (Data(i, 1)-PrevMeant(1, 1))^2+(Data(i, 2)-PrevMeant(1, 2))^2+ (Data(i, 3)-PrevMeant(1,3))^2;
      
       D2 = (Data(i, 1)-PrevMeant(2, 1))^2+(Data(i, 2)-PrevMeant(2,2))^2+ (Data(i, 3)-PrevMeant(3,3))^2;
       if InitMean3~=0
           D3 = (Data(i, 1)-PrevMeant(3, 1))^2+(Data(i, 2)- PrevMeant(3,2))^2+ (Data(i, 3)-PrevMeant(3, 3))^2;
       end
       
    %%%%%%%%%%%%%Identify the minimum distance from the sample to the means%%%%%%%%%%%%%%% 
       if InitMean3~=0
           Dim=[D1 D2 D3];
           d=(min(Dim));
           l=find(Dim==d);
           Index=l;
           %[~,Index] = [~,l];
       else
           Dim=[D1 D2];
           d= find(min(Dim));
           Index=d;
          % [~,Index] = d;
       end
    %%%%%%%%%%%%Label the data samples based on the minimum euclidean distance%%%%%%%%%%%%  
       Label(i) = Index;
    end
    %%^this loop coputes the mu for each sample
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%Compute the new means%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Indexdata(1: length(Data),1)= Label;
    Indexdata(1: length(Data), 2:4)= Data;
    
    %%% calculating the mean for each class
    
    ind1 = Indexdata(:,1) == 1;
    ind2= Indexdata(:, 1) == 2;
    if InitMean3~=0
        ind3= Indexdata(:, 1) ==3;
    end
    
    A1 = Indexdata(ind1,:); %%a1 stores all of class 1
    [LenA1,~]=size(A1);
    A2 = Indexdata(ind2,:);% a2 stores all of class 2
    [LenA2, ~]=size(A2);
    if InitMean3~=0
        A3=Indexdata(ind3, :);% a3 stores all of class 3
        [LenA3,~]=size(A3);
    end
    %calculates average of mu1 mu2 and mu3
    total1 = sum(A1);
    FinalMean1=total1./LenA1;
    total2 = sum(A2);
    FinalMean2=total2./LenA2;
    if InitMean3~=0
        total3 = sum(A3);
        FinalMean3=total3./LenA3;
    else
        FinalMean3= 0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%Check for criterion function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%If criteria not met repeate the above%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %below- breaks if the mean didn't change between cycles
   
   if FinalMean1(1,2:4)==PrevMeant(1,:)
       break;
   end
   if (FinalMean2(1,2:4)==PrevMeant(2,:))
           break; 
   end
   if (InitMean3~=0)
       if FinalMean3(1,2:4)==PrevMeant(3, :)
           break;
       end
   end
%     if InitialMean3~=0
%         CurrMean =
%     else
%         CurrMean =
%     end
    if Itr==MaxNoOfIteration % Check conditions to stop the algorithm.
        break;
    end
    PrevMeant = [FinalMean1(1,2:4); FinalMean2(1,2:4); FinalMean3(1,2:4)]; 
   if Itr<2
       scatter3(PrevMeant(:,1), PrevMeant(:,2), PrevMeant(:,3), 'g');
   end
   
end
end

%Leave-One-Out and Bootstrap Methods
% Use Lab2 solution with DataLab2_1.mat
