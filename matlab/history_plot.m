
    if mod(n_steps,100) == 0
        for p=1:nP
            pause(0.1)
            figure (2)
            plot_p = p;
            plot3(xHistory(1:n_steps,plot_p),yHistory(1:n_steps,plot_p),zHistory(1:n_steps,plot_p))
            xlabel('x axis')
            ylabel('y axis')
            zlabel('z axis')
            title('Boris Method')
            legend(num2str(time))
            hold on
            
            % figure (3)
            % plot(t_boris(1:n_steps),coll_hist(1:n_steps,1),t_boris(1:n_steps),coll_hist(1:n_steps,3))
            % figure (4)
            % plot(t_boris(1:n_steps),coll_hist(1:n_steps,2))
            % figure (5)
            % plot(t_boris(1:n_steps),coll_hist(1:n_steps,4),t_boris(1:n_steps),coll_hist(1:n_steps,5))
            % figure (6)
            % h1=histogram(coll_hist(1:n_steps,6))%,t_boris(1:n_steps),coll_hist(1:n_steps,7)
            % hold on
            % h2 = histogram(coll_hist(1:n_steps,7))
            % hold off
            pause(0.1)
        end
    end
    hold off