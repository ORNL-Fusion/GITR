constants

%Ionization
file_inz = 'scd50_w.dat';
%Recombination
file_rcmb = 'acd50_w.dat';

[IonizationTemp, IonizationDensity, IonizationRateCoeff, IonizationChargeState] = ADF11(file_inz);
[RecombinationTemp, RecombinationDensity, RecombinationRateCoeff, RecombinationChargeState] = ADF11(file_rcmb);

IonizationData.Temp = 10.^IonizationTemp;
IonizationData.Density = 10.^IonizationDensity.*1e6;
IonizationData.RateCoeff = 10.^IonizationRateCoeff./1e6;
IonizationData.ChargeState = IonizationChargeState;


RecombinationData.Temp = 10.^RecombinationTemp;
RecombinationData.Density = 10.^RecombinationDensity.*1e6;
RecombinationData.RateCoeff = 10.^RecombinationRateCoeff./1e6;
RecombinationData.ChargeState = RecombinationChargeState;


save('Processed/W/IonizationData.mat','IonizationData');
save('Processed/W/RecombinationData.mat','RecombinationData');

