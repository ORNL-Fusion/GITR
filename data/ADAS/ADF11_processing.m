clear all
close all

addpath(genpath('../../matlab'));

constants

Z = 13;

%Ionization
file_inz = 'scd89_al.dat';
%Recombination
file_rcmb = 'acd89_al.dat';

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


save('Processed/Al/IonizationData.mat','IonizationData');
save('Processed/Al/RecombinationData.mat','RecombinationData');

%Open the file
ncid = netcdf.create(['ADAS_Rates_Al.nc'],'NC_WRITE')
 
%Define the dimensions
dimScalar = netcdf.defDim(ncid,'scalar',1);
dimPair = netcdf.defDim(ncid,'pair',2);
dimTemp_Ionize = netcdf.defDim(ncid,'n_Temperatures_Ionize',length(IonizationData.Temp));
dimDensity_Ionize = netcdf.defDim(ncid,'n_Densities_Ionize',length(IonizationData.Density));
dimTemp_Recombine = netcdf.defDim(ncid,'n_Temperatures_Recombine',length(RecombinationData.Temp));
dimDensity_Recombine = netcdf.defDim(ncid,'n_Densities_Recombine',length(RecombinationData.Density));

dimChargeState_Ionize = netcdf.defDim(ncid,'n_ChargeStates_Ionize',length(IonizationData.ChargeState));
dimChargeState_Recombine = netcdf.defDim(ncid,'n_ChargeStates_Recombine',length(RecombinationData.ChargeState)); 
%Define IDs for the dimension variables (pressure,time,latitude,...)
Z_ID = netcdf.defVar(ncid,'Atomic_Number','int',[dimScalar]);
TempGridIonization=netcdf.defVar(ncid,'gridTemperature_Ionization','double',[dimTemp_Ionize]);
DensityGridIonization=netcdf.defVar(ncid,'gridDensity_Ionization','double',[dimDensity_Ionize]);
TempGridRecombination=netcdf.defVar(ncid,'gridTemperature_Recombination','double',[dimTemp_Recombine]);
DensityGridRecombination=netcdf.defVar(ncid,'gridDensity_Recombination','double',[dimDensity_Recombine]);

ChargeStateGridIonization = netcdf.defVar(ncid,'gridChargeState_Ionization','double',[dimChargeState_Ionize dimPair]);
ChargeStateGridRecombination = netcdf.defVar(ncid,'gridChargeState_Recombination','double',[dimChargeState_Recombine dimPair]);

 
%Define the main variable ()
IonizeCoeff = netcdf.defVar(ncid,'IonizationRateCoeff','double',[dimDensity_Ionize dimTemp_Ionize dimChargeState_Ionize]);
RecombineCoeff = netcdf.defVar(ncid,'RecombinationRateCoeff','double',[dimDensity_Recombine dimTemp_Recombine dimChargeState_Recombine]);

 
%We are done defining the NetCdf
netcdf.endDef(ncid);
 
%Then store the dimension variables in
netcdf.putVar(ncid,Z_ID,Z);
netcdf.putVar(ncid,TempGridIonization,log10(IonizationData.Temp));
netcdf.putVar(ncid,DensityGridIonization,log10(IonizationData.Density));
netcdf.putVar(ncid,TempGridRecombination,log10(RecombinationData.Temp));
netcdf.putVar(ncid,DensityGridRecombination,log10(RecombinationData.Density));

netcdf.putVar(ncid,ChargeStateGridIonization,IonizationData.ChargeState);
netcdf.putVar(ncid,ChargeStateGridRecombination,RecombinationData.ChargeState);

 
%Then store my main variable
netcdf.putVar(ncid,IonizeCoeff,log10(IonizationData.RateCoeff));
netcdf.putVar(ncid,RecombineCoeff,log10(RecombinationData.RateCoeff));

 
 
%We're done, close the netcdf
netcdf.close(ncid)

Coeff = interpn(log10(IonizationData.Density),log10(IonizationData.Temp),IonizationData.RateCoeff(:,:,1),log10(2.5e19),log10(15),'linear',0)
tion = 1/(Coeff*2.5e19);