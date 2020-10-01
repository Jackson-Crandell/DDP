function [ ] = invPend_animation(x,horizon)
save_video = false;
O=[0 0];
axis(gca,'equal');
axis([-1.5 1.5 -1.5 1.5]);
grid on
for i=1:horizon
    P=1*[sin(x(1,i)) -cos(x(1,i))];
    
    pend=line([O(1) P(1)],[O(2) P(2)], 'LineWidth', 4);
    if (save_video)
        F(i) = getframe(gcf);
    end
    pause(0.01);
    if i<horizon
        
        delete(pend);
        
    end
end

if (save_video)
    video = VideoWriter('InvPend.avi','Uncompressed AVI');
    open(video)
    writeVideo(video,F)
    close(video)
end

end

