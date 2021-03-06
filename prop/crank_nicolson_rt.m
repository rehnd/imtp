function psi_n = crank_nicolson_rt(H,en,psi,nsteps,dt,plots)
%Implicit (Crank Nicolson) in Real time
%    Propagates solution according to
%     (I + i*dt*H)\psi(t+dt) = (I-i*dt*H)\psi(t)
%      doing a linear solve

fprintf('\nRunning Crank Nicolson in Imaginary time\n')

% Set psi old to current psi
psi_o  = psi;
psi_ex = psi; % Exact solution 
err = zeros(nsteps,4); % error vectors
psi2 = zeros(nsteps,4);
X = meshgrid(diag(en));

% Movie setup crap
f = figure('Visible',plots.showm);
plot(psi(:,1:4));
axis tight manual;
ax = gca;
ax.NextPlot = 'replaceChildren';
f.Renderer = 'zbuffer';
Figs(nsteps) = struct('cdata',[],'colormap',[]);

for k = 1:nsteps
    A = speye(size(H)) + 0.5i*dt*H;
    B = speye(size(H)) - 0.5i*dt*H;
    b = B*psi_o;
    psi_n = A\b;

    psi_o = psi_n;

    psi2(k,1:4) = dot(psi_o(:,1:4),psi_o(:,1:4));
    psi_ex(:,1:4) = psi_ex(:,1:4).*exp(-1i*X(:,1:4)*dt);

% %    Should not need to normalize for CN method
%    if mod(k,1) == 0
%        psi_o = psi_o./meshgrid(sqrt(sum(abs(psi_o).^2,1)));
%    end

    err(k,1:4) = sqrt(sum(abs((psi_ex(:,1:4) - psi_n(:,1:4))).^2,1));

    if strcmp(plots.showm, 'on') || plots.savem
        plot(real(psi_o(:,1:4)));
        drawnow;
        Figs(k) = getframe;
    end
end
fprintf('Done. Ran for %i time steps\n',nsteps)

plot_err(err,psi2,6,plots,nsteps);

% Options for showing/saving movies
if strcmp(plots.showm,'on')
    movie(f,Figs,2,30)
end
if plots.savem
    fprintf('Saving movie data. Could take a while...\n')
    wObj = VideoWriter(['mov/',plots.fn{6},'.avi']);
    open(wObj);
    for j = 1:nsteps
        writeVideo(wObj,Figs(j));
    end
    close(wObj);
    fprintf('Done writing %s\n',plots.fn{6})
end

end