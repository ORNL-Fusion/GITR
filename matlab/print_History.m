if trackHistory
if printHistory
dlmwrite('xHistory_out.txt',xHistory,'delimiter','\t','precision',4)
dlmwrite('yHistory_out.txt',yHistory,'delimiter','\t','precision',4)
dlmwrite('zHistory_out.txt',zHistory,'delimiter','\t','precision',4)
dlmwrite('vxHistory_out.txt',vxHistory,'delimiter','\t','precision',4)
dlmwrite('vyHistory_out.txt',vyHistory,'delimiter','\t','precision',4)
dlmwrite('vzHistory_out.txt',vzHistory,'delimiter','\t','precision',4)
dlmwrite('Z_History_out.txt',Z_History,'delimiter','\t','precision',1)
end
end