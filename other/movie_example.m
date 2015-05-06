
fig = figure('Visible','off');
Z = peaks;
surf(Z)
axis tight manual
ax = gca;
ax.NextPlot = 'replaceChildren';
set(gcf,'Renderer','zbuffer');

loops = 40;
F(loops) = struct('cdata',[],'colormap',[]);
for j = 1:loops
    X = sin(j*pi/10)*Z;
    surf(X,Z)
    drawnow
    F(j) = getframe;
end

writerObj = VideoWriter('peaks.avi');
open(writerObj);
for j = 1:loops
    writeVideo(writerObj,F(j));
end
close(writerObj)

%movie(fig,F,2)