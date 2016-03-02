figure(1)
%i=27
for i=1:iCPart
        TrajP(:,:)=TrajAr{1,1,i};
        color=[rand,rand,rand]
        plot3(TrajP(:,1),TrajP(:,2),TrajP(:,3),'Color',color)
        hold on
        clear TrajP;
end