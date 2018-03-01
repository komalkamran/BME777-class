% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % BME777: LAB 3: Multilayer Neural Networks.
% % Statlog (Heart) Dataset: https://archive.ics.uci.edu/ml/datasets/Statlog+%28Heart%29
% % The first two features are contained in the first two columns.
% % 1st feature: Resting blood pressure.
% % 2nd feature: Oldpeak = ST depression induced by exercise relative to rest.
% % The third column contains the label information.
% % Class +1: Absence of heart disease.
% % Class -1: Presence of heart disease.
% % 50 samples were extracted for each class.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Inputs:
% % 1. Eta: Learning rate.
% % 2. Theta: Threhold for the cost function to escape the algorithm.
% % 3. MaxNoOfIteration: Maximum number of iteration.
% % 4. Problem: 1: XOR, otherwise: Classification problem with given dataset.
% % 5. Data: the dataset used for training NN when problem ~=1.
% % Outputs:
% % 1. J: an array of cost.
% % 2. w: trained weight matrix.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Example of use:
% 
% %load DataLab3; 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [J,w] = Lab3Data(Eta ,Theta, MaxNoOfIteration, Problem, Data)
% 
% 
% %Initialization.
% 
function [J,w] = lab3(Eta ,Theta, MaxNoOfIteration, Problem, Data)

%Initialization.
% clear all;
% Eta = 0.1;
% Theta = 0.001;
% MaxNoOfIteration = 300;
% Problem = 1;
% load ('DataLab3.mat');
% Data = DataLab3;

 if Problem == 1 %XOR
    wih1(1,1:3) =[0.69 0.39 0.41]; %weight vector input to hidden unit no.1. %small converge faster
    wih2(1,1:3) =[0.65 0.83 0.37]; %weight vector input to hidden unit no.2.
    who1(1,1:3) =[0.42 0.59 0.56]; %weight vector hidden to output unit. 
    
    % Add data to feature 1,2 and label vectors.
    x1 = [-1 -1 1 1];
    x2 = [-1 1 -1 1];
    t = [-1 1 1 -1];
 else
     wih1(1,1:3) =[1.69 1.39 1.41]; %weight vector input to hidden unit no.1. %small converge faster
    wih2(1,1:3) =[1.65 1.83 1.37]; %weight vector input to hidden unit no.2.
    who1(1,1:3) =[1.42 1.59 1.56];%weight vector hidden to output unit. 
    
    % Add data to feature 1,2 and label vectors.
     x1 = Data(:,1)';
     x2 = Data(:,2)';
     t = Data(:,3)';
     end

% Initialize number of iteration and cost.
J = zeros(MaxNoOfIteration,1);
a=1; b=1;
r = 0; y=0;
Xm=zeros(length(x1),3);
Xm(:,1) =1;
Xm(:,2) =x1;
Xm(:,3) =x2;
Z = zeros(1, length(x1)); 
Y=zeros(1,3);
Y(1,1)=1;%this is the bias
    
while(1)
    
    
    % Initialize gradients of the three weight vectors.
    DeltaWih1 = 0; % Inputs of bias, x1,x2 to hidden neuron 1.
    DeltaWih2 = 0; % Inputs of bias, x1,x2 to hidden neuron 2.
    DeltaWho1 = 0; % Inputs of bias, y1,y2 to output neuron.
    % Initialize training sample order and predicted output.
    m = 0;
    
    while(m<length(x1))
        
        m = m + 1; %m indexes the paticular pattern presentation
       
        netj_1=Xm(m,:)*wih1.';
        netj_2=Xm(m,:)*wih2.';
        Y(1,2)=tansig(netj_1);
        Y(1,3)=tansig(netj_2);
        netk_1=who1*Y.';
        Z(m) =tansig(netk_1);

        
        %Derivatives of Sigmoids
        d_sigmoidk1=1-(Z(m)^2);
        d_sigmoidj1=1-(Y(1,2)^2);
        d_sigmoidj2=1-(Y(1,3)^2);
       
        % Calculate the sensitivity value of each hidden neuron and the output neuron.
        DeltaO1 =(t(m)-Z(m))*d_sigmoidk1; % Sensitivity value of the output neuron.
        DeltaH1 =d_sigmoidj1*who1(2)*DeltaO1;% Sensitivity value of hidden neuron 1.
        DeltaH2 =d_sigmoidj2*who1(3)*DeltaO1; % Sensitivity value of hidden neuron 2.
        
        % Update the gradient.
         DeltaWih1=DeltaWih1+Eta*DeltaH1*Xm(m,:);
         DeltaWih2=DeltaWih2+Eta*DeltaH2*Xm(m,:);
         DeltaWho1=DeltaWho1+Eta*DeltaO1*Y;
                    
    end
    % Update the weight vectors.
        wih1=wih1+DeltaWih1; % Weight vector input to hidden unit no.1
        wih2=wih2+DeltaWih2; % Weight vector input to hidden unit no.2
        who1=who1+DeltaWho1; % Weight vector hidden to output unit.
              
    % Check the condition to stop.
      r=r+1;
      J(r)=1/2.*sum((Z-t).^2);
    
      if (J<Theta)
          break;
      end

    if(r==MaxNoOfIteration)
        break;
    end  
end

w = [wih1; wih2; who1];
%Learning Curve


plot (J);
ylabel('J(r) Cost Function'); 
xlabel('# of Epoch'); 
title('Learning Curve'); 

%Decision Boundary
x_d1 = -1:0.01:1;
x_d2= -1:0.01:1;
[X1, X2] = meshgrid(x_d1,x_d2); 

LengthOfData = length(x_d1)*length(x_d2); 
XS = zeros(LengthOfData,3); 
XS(:,1) = 1; %bias
n=1; 
for i = 1:length(x_d1) 
    for j = 1:length(x_d2) 
        XS(n,2) = X1(i,j);
        XS(n,3) = X2(i,j);
        n=n+1;
    end
end


y(1,1)=1; %bias 
tot_len = length(XS);

z_d = zeros(tot_len,1);
a_y = zeros(tot_len,1);
b_y = zeros(tot_len,1);
    for i = 1:tot_len
        netj_1=(wih1)*XS(i,:).';
        netj_2=(wih2)*XS(i,:).';
        y(1,2)=(exp(netj_1)-exp(-netj_1))/(exp(netj_1)+exp(-netj_1));
        y(1,3)=(exp(netj_2)-exp(-netj_2))/(exp(netj_2)+exp(-netj_2));
        netk=(who1)*y.';
        zk=tanh(netk); 
        
        if zk> 0 
            z_d(i) = 1;
        else 
            z_d(i) = -1;
        end
        
        a_y(i)=y(1,2);
        b_y(i)=y(1,3);
       
    end
   X11 =XS(:,2);
   X22 =XS(:,3); 
   
   Ind1 = z_d>0; 
   Ind2 = z_d<=0; 
   z_d(Ind1) = 1; 
   z_d(Ind2) = -1;
figure; 
scatter3(X11(Ind1),X22(Ind1), z_d(Ind1)); 
hold on; 
scatter3(X11(Ind2),X22(Ind2), z_d(Ind2)); 
(view(0,90));
title('X-Space'); 
ylabel('X2'); 
xlabel('X1'); 

figure;
scatter3(a_y(Ind1), b_y(Ind1), z_d(Ind1)); 
hold on; 
scatter3(a_y(Ind2), b_y(Ind2), z_d(Ind2)); 
(view(0,90)); 
title('Y-Space'); 
ylabel('Y2'); 
xlabel('Y1'); 

