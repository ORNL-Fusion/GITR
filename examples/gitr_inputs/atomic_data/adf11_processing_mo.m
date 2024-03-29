close all
clear all

ME = 9.10938356e-31;
MI = 1.6737236e-27;
Q = 1.60217662e-19;
EPS0 = 8.854187e-12;
Z = 42;
%Ionization
file_inz = 'scd89_mo.dat';
%Recombination
file_rcmb = 'acd89_mo.dat';

[IonizationTemp, IonizationDensity, IonizationRateCoeff, IonizationChargeState] = ADF11s(file_inz);
[RecombinationTemp, RecombinationDensity, RecombinationRateCoeff, RecombinationChargeState] = ADF11a(file_rcmb);

IonizationData.Temp = IonizationTemp;
IonizationData.Density = log10(10.^IonizationDensity.*1e6);
IonizationData.RateCoeff = log10(10.^IonizationRateCoeff./1e6);
IonizationData.ChargeState = IonizationChargeState;


RecombinationData.Temp = RecombinationTemp;
RecombinationData.Density = log10(10.^RecombinationDensity.*1e6);
RecombinationData.RateCoeff = log10(10.^RecombinationRateCoeff./1e6);
RecombinationData.ChargeState = RecombinationChargeState;


% save('Processed/W/IonizationData.mat','IonizationData');
% save('Processed/W/RecombinationData.mat','RecombinationData');

%Open the file
ncid = netcdf.create(['ADAS_Rates_Mo.nc'],'NC_WRITE')
 
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
netcdf.putVar(ncid,TempGridIonization,IonizationData.Temp);
netcdf.putVar(ncid,DensityGridIonization,IonizationData.Density);
netcdf.putVar(ncid,TempGridRecombination,RecombinationData.Temp);
netcdf.putVar(ncid,DensityGridRecombination,RecombinationData.Density);

netcdf.putVar(ncid,ChargeStateGridIonization,IonizationData.ChargeState);
netcdf.putVar(ncid,ChargeStateGridRecombination,RecombinationData.ChargeState);

 
%Then store my main variable
netcdf.putVar(ncid,IonizeCoeff,IonizationData.RateCoeff);
netcdf.putVar(ncid,RecombineCoeff,RecombinationData.RateCoeff);

 
 
%We're done, close the netcdf
netcdf.close(ncid)

% Example calculation of mean free path of ionization
impurity_mass = 96; % amu
impurity_charge = 0;
impurity_kinetic_energy = 10; % eV

impurity_velocity = sqrt(2*impurity_kinetic_energy*Q/impurity_mass/MI);

electron_temperature = 20; %eV
electron_density = 10^16; % m-3

Coeff = interpn(IonizationData.Density,IonizationData.Temp,IonizationData.RateCoeff(:,:,impurity_charge+1),log10(electron_density),log10(electron_temperature),'linear',0);

ionization_time = 1/(10^Coeff*electron_density); % in seconds

mean_free_path = ionization_time*impurity_velocity; % in m