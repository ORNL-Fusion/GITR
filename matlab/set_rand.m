n_MC_process = 6;
random = zeros(nP,nT,n_MC_process);

for i=1:nP


[s1,s2,s3,s4,s5,s6] = RandStream.create('mrg32k3a','NumStreams',n_MC_process,'Seed','shuffle'); %Include ,'Seed','shuffle' to get different values each time
%Ionization
random(i,:,1) = rand(s1,nT,1);
%Recombination
random(i,:,2) = rand(s2,nT,1);
%Cross-field diffusion direction
random(i,:,3) = rand(s3,nT,1);
%Parallel velocity diffusion +/-
random(i,:,4) = rand(s4,nT,1);
%Perpendicular velocity #1 diffusion +/-
random(i,:,5) = rand(s5,nT,1);
%Perpendicular velocity #2 diffusion +/-
random(i,:,6) = rand(s6,nT,1);

end