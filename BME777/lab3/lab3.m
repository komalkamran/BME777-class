%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BME777: LAB 3: Multilayer Neural Networks.
% Statlog (Heart) Dataset: https://archive.ics.uci.edu/ml/datasets/Statlog+%28Heart%29
% The first two features are contained in the first two columns.
% 1st feature: Resting blood pressure.
% 2nd feature: Oldpeak = ST depression induced by exercise relative to rest.
% The third column contains the label information.
% Class +1: Absence of heart disease.
% Class -1: Presence of heart disease.
% 50 samples were extracted for each class.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% 1. Eta: Learning rate.
% 2. Theta: Threhold for the cost function to escape the algorithm.
% 3. MaxNoOfIteration: Maximum number of iteration.
% 4. Problem: 1: XOR, otherwise: Classification problem with given dataset.
% 5. Data: the dataset used for training NN when problem ~=1.
% Outputs:
% 1. J: an array of cost.
% 2. w: trained weight matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example of use:

%load DataLab3; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [J,w] = lab3(Eta ,Theta, MaxNoOfIteration, Problem, Data)


%Initialization.

if Problem == 1
    wih1 = [.69 .39 .41]; %weight vector input to hidden unit no.1.
    wih2 = [.65 .83 .37]; %weight vector input to hidden unit no.2.
    who1 = [.42 .59 .56]; %weight vector hidden to output unit. 
    
    % Add data to feature 1,2 and label vectors.
    x1 = [-1 -1 1 1]; %Input1
    x2 = [-1 1 -1 1]; %Input 2
    t = [-1 1 1 -1]; %target tk
else
    wih1 = [.69 .39 .41]; %weight vector input to hidden unit no.1.
    wih2 = [.65 .83 .37]; %weight vector input to hidden unit no.2.
    who1 = [.42 .59 .56]; %weight vector hidden to output unit.
    % Add data to feature 1,2 and label vectors.
    
    
    x1 = Data(:,1)';
    x2 = Data(:,2)';
    t = Data(:,3)';
end
    
% Initialize number of iteration and cost.
r = 0;
J = zeros(MaxNoOfIteration,1);

while(1);
    
    r=r+1;
    
    % Initialize gradients of the three weight vectors.
    DeltaWih1 = [0 0 0]; % Inputs of bias, x1,x2 to hidden neuron 1.
    DeltaWih2 = [0 0 0]; % Inputs of bias, x1,x2 to hidden neuron 2.
    DeltaWho1 = [0 0 0]; % Inputs of bias, y1,y2 to output neuron.
    
    % Initialize training sample order and predicted output.
    m = 0;
    Z = zeros(1,length (x1));
    y1vec= zeros(1, length (x1));
    y2vec= zeros(1, length (x1));
    
    while(m<length(x1))
        
        m = m + 1;
        
       
        Xm = [1 x1(1,m) x2(1,m)] ; %assume bias = 1 
        c1 = wih1 * Xm'; %net activation of h1
        netj_1 = sum(c1);
        y1 =  tanh(netj_1);
        c2= wih2 * Xm' ; %net activation of h2
        netj_2 = sum(c2);
        y2 = tanh(netj_2);
        
        %inputs of y
        Ym  =  [1 y1 y2];
        yk = who1*Ym';
        net_k = sum(yk);
        zk= tanh(net_k);
        Z(m)= zk;
        tk = t(:,m);
        
      
        
        % Calculate the sensitivity value of each hidden neuron and the output neuron.
        DeltaO1 = (tk - zk)* (1 - zk^2);% Sensitivity value of the output neuron.
        DeltaH1 = who1(2)*DeltaO1*(1-(y1^2)); % Sensitivity value of hidden neuron 1.
        DeltaH2 = who1(3)*DeltaO1*(1-(y2^2)); % Sensitivity value of hidden neuron 2.
        
    DeltaWih1_temp = DeltaH1.*Xm.*Eta; % 
    DeltaWih2_temp = DeltaH2.*Xm.*Eta;
    DeltaWho1_temp = DeltaO1.*Ym.*Eta;
        
        % Update the gradient.
        DeltaWih1 = DeltaWih1 + DeltaWih1_temp;
        DeltaWih2 = DeltaWih2 + DeltaWih2_temp;
        DeltaWho1 = DeltaWho1 + DeltaWho1_temp;          
        
    end
    
    % Update the weight vectors.
    wih1 = wih1 + DeltaWih1; % Weight vector input to hidden unit no.1
    wih2 = wih2 + DeltaWih2; % Weight vector input to hidden unit no.2
    who1 = who1 + DeltaWho1; % Weight vector hidden to output unit.
    
    % Check the condition to stop.
    J(r)= 0.5*(tk-zk)^2;
%  grad_J = -1*(DeltaO1)*Ym;
%  grad_JCheck = sqrt (grad_J(1,1)^2 +grad_J (1,2)^2 +grad_J(1,3)^2);
%  
%     if(grad_JCheck <Theta)
%         break;
%     end
        
    if(J(r)==Theta)
        break;
    end
    if(r == MaxNoOfIteration)
        break;
    end

end


%wih1

%wih2

%who1

 figure;
 plot (J);
 xlabel ('Number of Iterations');
 ylabel ('Error');
 


w = [wih1; wih2; who1];    
    
 

figure;
gscatter(x1, x2, t, 'rg', 'os', 7);
xlabel('x1');
ylabel('x2');
title('graph in xspace');
hold on;

x1_temp=(-length(x1): length(x1));
x2_temp=(-wih1(:, 1)-wih1(:, 2)*x1_temp)/wih1(:,3);
plot(x1_temp, x2_temp);
hold on;

x2_temp=(-wih2(:, 1)-wih2(:, 2)*x1_temp)/wih2(:, 3);
plot(x1_temp, x2_temp);
title('decision surface');

 hold on;
  y1_temp=(-length(y1):length(y1));
  y2_temp=(-who1(:,1)-who1(:,2)*y1_temp)/(who1(:,3));
  figure;
  %   %%%%%%
  n=1:length(y1);
  gscatter(y1(1,n),y2(1,n),t,'rg','os',7); %axis([-2 2 -2 6]);
  hold on;
  plot(y1_temp,y2_temp);title('Decision Surface y');

  
 % Classification Accuracy 
  if Problem == 2
      X_test = ones(3,length(x1)); Y_test = ones(3,length(x1));
      X_test(2,:) = x1;
      X_test(3,:) = x2;
      t_test = t';

      
      Y_test(1,:) = w(1,:)*Y_test + w(2,:)*Y_test;
      Y_test(2,:) = w(1,:)*X_test;
      Y_test(3,:) = w(2,:)*X_test; %go from input to hidden...y=wij^t*x

      
      Z_test = w(3,:)*Y_test; %hidden to output...z=wjk^t*y
      Z_test = ((w(1,:)+w(2,:)+w(3,:))*X_test);
      count = 0;

      
      q = length(Z_test);

      
      for i = 1:q
          targVout =  (Z_test(1,q))*t_test(q,1); %The computed ouput and target output be the same sign, give positive result
          if targVout > 0    
            count = count + 1;
          end
      end

  end    
      Accuracy = (count/length(x1))*100;
      disp('Accuracy of the Classifier is: ');
      disp(Accuracy);
  end


