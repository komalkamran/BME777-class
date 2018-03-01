%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BME777: LAB 1: Bayesian Decision Theory.

% The data for this lab is extracted from Pima Indians Diabetes Data Set: 
% https://archive.ics.uci.edu/ml/datasets/Pima+Indians+Diabetes
% The first two columns contain the 2nd and 3rd features of the original dataset. 
% 1st feature: Plasma glucose concentration.
% 2nd feature: Diastolic blood pressure (mm Hg).
% The third colum contatins the labels (1: positive, 2: negative) for diabetes.
% 268 samples were extracted for each class.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% 1. FeatureX: Feature value to be tested (to identify which class it belongs to).
% 2. Data: Matrix containing the training feature samples and class labels.
% 3. FeatureForClassification: Select type of feature used for
% classification. (1 or 2)
% 4. LabelColumn: Specify the column containing the labels of the data.
% Outputs:
% 1. PosteriorProbabilities: Posterior probabilities of class 1 and 2 given FeatureX.
% 2. DiscriminantFunctionValue: Value of the discriminant function given FeatureX.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example of use:
% load Diabetes.mat;
% FeatureX = 5;
% LabelColumn = 3;
% FeatureForClassification = 1;
% [PosProb, G] = lab1(FeatureX, Diabetes,FeatureForClassification, LabelColumn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PosteriorProbabilities,DiscriminantFunctionValue]=lab1(FeatureX,Data,FeatureForClassification, LabelColumn)

% Get number of samples.
[ro,~] = size(Data);

% Select feature for classification (1 or 2).  
SelectedFeature=Data(:,FeatureForClassification);

% Get class labels.
Label=Data(:,LabelColumn); 

ind1 = Data(:,3) == 1;
ind2 = Data(:,3) == 2;
Data1 = Data(ind1, :);
Data2 = Data(ind2, :);

%%%%%%%%Plot the histogram and box plot of features related to two classes%%%%%%%%%%

% Plot hist.
histogram(Data1(:,FeatureForClassification));
hold on
histogram(Data2(:,FeatureForClassification));
xlabel('Feature value');ylabel('Frequency');title('Histogram for selected feature');legend('Class 1','Class 2');
hold off
% Plot boxplot.
figure;
group = [ones(size(Data1(:,FeatureForClassification))); 2 * ones(size(Data2(:,FeatureForClassification)));];
boxplot([Data1(:,FeatureForClassification);Data2(:,FeatureForClassification)], group);
xlabel('Class');ylabel('Concentration');title('Boxplot for selected feature');
    
%%%%%%%%%%%%%%%%%%%%%%%Compute Prior Probabilities%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate prior probability of class 1 and class 2.

Pr1= length(Data1(:,FeatureForClassification))/length(Data(:,FeatureForClassification));
Pr2= length(Data2(:,FeatureForClassification))/length(Data(:,FeatureForClassification));

%%%%%%%%%%%%%%%%%%%%Compute Class-conditional probabilities%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the mean and the standard deviation of the class conditional density p(x/w1).
m11= mean(Data1(:,FeatureForClassification));
std11= std(Data1(:,FeatureForClassification));

% Calculate the mean and the standard deviation of the class conditional density p(x/w2).
m12= mean(Data2(:,FeatureForClassification));
std12= std(Data2(:,FeatureForClassification));

% Calculate the class-conditional probability of class 1 and class 2.
cp11= 1/sqrt((2*pi*std11))*exp(-(FeatureX'-m11).^2/(2*std11^2));
cp12= 1/sqrt((2*pi*std12))*exp(-(FeatureX'-m12).^2/(2*std12^2));

%%%%%%%%%%%%%%%%%%%%Compute feature probabilities%%%%%%%%%%%%%%%%%%%%%%%%%

pdf = (cp11)*(Pr1)+(cp12)*(Pr2);

%%%%%%%%%%%%%%%%%%%%%%%Compute the posterior probabilities%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Posterior probabilities for the test feature');

pos11= (cp11.*Pr1)./pdf;
pos12= (cp12.*Pr2)./pdf;

PosteriorProbabilities = [pos11,pos12];
PPt = PosteriorProbabilities'; 

fpErr = 1;
tnErr = 1; 
loss = [0 fpErr; tnErr 0]; 
riskFunc = loss*PPt;
minError = ones(size(riskFunc)) - riskFunc;

%%%%%%%%%%%Compute Discriminant function Value for Min Error Rate Classifier%%%%%%%%
disp('Discriminant function value for the test feature');

DiscriminantFunctionValue = minError(1,:) - minError(2,:); 
