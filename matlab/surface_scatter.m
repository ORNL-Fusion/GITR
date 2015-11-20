figure(10)

for i=1:nP
   
    plot(particles(i).y,particles(i).z,'*')
    hold on
end

xlabel('y axis')
ylabel('z axis')

title('Particles on Surface')
hold off