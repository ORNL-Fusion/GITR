
t_boris = [0:dt:nT*dt];
if trackHistory
xHistory = zeros(nT, nP);
yHistory = zeros(nT, nP);
zHistory = zeros(nT, nP);
vxHistory = zeros(nT, nP);
vyHistory = zeros(nT, nP);
vzHistory = zeros(nT, nP);
Z_History = zeros(nT, nP);

coll_hist = zeros(nT,7);
end

end_pos = zeros(nP,7);