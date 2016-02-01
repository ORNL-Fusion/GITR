
t_boris = [0:dt:nT*dt];

if trackHistory
    
    aHistory.x = 0;
    aHistory.y = 0;
    aHistory.z = 0;
   
    aHistory.vx = 0;
    aHistory.vy = 0;
    aHistory.vz = 0;
   
    aHistory.Z = 0;
    
    history(nT,nP) = aHistory;
    
end