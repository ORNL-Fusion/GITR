struct ionize { 
  ...//Accept inputs including Temp, Dens, particles, dt etc.
  CUDA_CALLABLE_MEMBER_DEVICE 
  void operator()(std::size_t indx) const { 
    float ionTemperature = getIonTemperature()
    float ionDensity = getIonDensity()
    float P = expf(-dt/ionTemperature);
    float P1 = 1.0-P;
    float r1 = rand()
    float chargeAdd = 0.0;
    if(r1 <= P1){
      chargeAdd = 1.0;}
    particles->charge[indx] = particles->charge[indx]+chargeAdd;
  }
