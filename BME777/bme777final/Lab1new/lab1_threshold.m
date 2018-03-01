
DataAdmission = keyadmission_type_id;

DataProcedures = keynum_lab_procedures;
DataMedications = keynum_medications;
DataTime = keytime_in_hospital;

DataGender = (keygender); 
DataGendernew = ones(1,14358);
for j= 1: length(DataGender)
    d = DataGender{j};
    if ~strcmp(d,'Male')
        DataGendernew (1,j) = 2;
    else
        DataGendernew (1,j) = 1;
    end
end 
 DataGendernew = DataGendernew';               

 DataReadmission = (keyreadmitted);
 DataReadmissionNew = ones (1,14358);
 for k = 1:length (DataReadmission)
     n = DataReadmission{k};
    if  strcmp (n, 'NO')
        DataReadmissionNew (1,k) = 1;
    elseif  strcmp(n, '>30')            
        DataReadmissionNew (1,k) = 2;
    else
        DataReadmissionNew (1,k) = 3;
    end
 end

 DataReadmissionNew = DataReadmissionNew';
           
 
 
DataAge = (keyage);
DataAgeNew = ones (1, 14358);
for  y = 1:length (DataAge)
     h = DataAge {y};
     if strcmp (h, '[0-10)')
         DataAgeNew(1,y) = 1;
     elseif strcmp (h, '[10-20)')
         DataAgeNew(1,y) = 2;
        elseif strcmp (h, '[20-30)')
         DataAgeNew(1,y) = 3;
            elseif strcmp (h, '[30-40)')
            DataAgeNew(1,y) = 4;
            elseif strcmp (h, '[40-50)')
             DataAgeNew(1,y) = 5;
             elseif strcmp (h, '[50-60)')
         DataAgeNew(1,y) = 6;
         elseif strcmp (h, '[60-70)')
         DataAgeNew(1,y) = 7;
         elseif strcmp (h, '[70-80)')
         DataAgeNew(1,y) = 8;
         elseif strcmp (h, '[80-90)')
         DataAgeNew(1,y) = 9;
     else strcmp (h, '[90-100)')
         DataAgeNew(1,y) = 10;
     end
end
DataAgeNew = DataAgeNew';


DataResults = (keyA1Cresult);
DataResultsNew = ones (1, 14358);
for z = 1:length (DataResults)
    m = DataResults{z};
    if strcmp (m, 'None')
        DataResultsNew (1,z) = 1; 
    elseif strcmp(m, 'Norm')
            DataResultsNew (1,z) = 2;
    elseif strcmp (m, '>8')
                DataResultsNew (1,z) = 3;
            else strcmp (m, '>7')
                DataResultsNew (1,z) = 4;
    end
end

DataResultsNew = DataResultsNew';


 
FeatureForClassification = 1;
LabelColumn = 2;
FeatureX = 100;



label = ones(1,14358);
LabelDischarge = double(keydischarge_disposition_id');
for i=1:length(LabelDischarge)
    if LabelDischarge(i) ~= 1
        label(1,i) = 2;
    else
        label(1,i) = 1;
    end
end
    
label = label';
Data = [DataAgeNew,label];
[PosProb, G] = lab1(FeatureX, Data, FeatureForClassification, LabelColumn);

[difference, index_At_G_Equals_0] = min(abs(G));
x0 = FeatureX(index_At_G_Equals_0);


fprintf('The Optimal Threshold is approximately at: %.2f\n', x0);    
% fprintf('The G value is\n', G); 