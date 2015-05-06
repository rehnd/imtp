function psi_n = implicit_euler_it(H,en,psi,nsteps,dt,plots)
%Explicit Euler in Imaginary time
%    Propagates solution according to
%     \psi(t+dt) = \psi(t) - dt*H\psi(t)

fprintf('\nRunning Implicit Euler in Imaginary time\n')

% Set psi old to current psi
psi_o  = psi;
psi_ex = psi; % Exact solution 
err = zeros(nsteps,4); % error vectors


% Movie setup crap
f = figure('Visible',plots.showm);
plot(psi(:,1:4));
axis tight manual;
ax = gca;
ax.NextPlot = 'replaceChildren';
f.Renderer = 'zbuffer';
Figs(nsteps) = struct('cdata',[],'colormap',[]);

for k = 1:nsteps
    A = speye(size(H)) + dt*H;
    X = meshgrid(diag(en));
    psi_n = A\psi_o;
    psi_n(:,1:4) = psi_n(:,1:4).*exp((1-1i)*X(:,1:4)*dt);
    psi_o = psi_n;
    
    psi_ex(:,1:4) = psi_ex(:,1:4).*exp(-1i*X(:,1:4)*dt);
    
    err(k,1:4) = sqrt(sum(abs((psi_ex(:,1:4) - psi_n(:,1:4))).^2,1));
    
    if strcmp(plots.showm, 'on') || plots.savem
        plot(real(psi_o(:,1:4)));
        drawnow;
        Figs(k) = getframe;
    end
end
fprintf('Done. Ran for %i time steps\n',nsteps)

fe = figure('Visible',plots.showp);
loglog(err);
grid on;
ylim([min(err(:,1)),10]); xlim([0,nsteps]);
xlabel('Iteration Number')
ylabel('|\psi_{num}(t) - \psi_{ex}(t)|')
title('Error in Implicit Euler Imaginary Time')
saveas(fe,['plots/',plots.fn{3},'.png'])

% Options for showing/saving movies
if strcmp(plots.showm,'on')
    movie(f,Figs,2,30)
if plots.savem
    fprintf('Saving movie data. Could take a while...\n')
    wObj = VideoWriter([plots.fn{3},'.avi']);
    open(wObj);
    for j = 1:nsteps
        writeVideo(wObj,Figs(j));
    end
    close(wObj);
    fprintf('Done writing %s\n',plots.fn{3})
end

end