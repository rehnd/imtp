function psi_n = explicit_euler_it(H,en,psi,nsteps,dt,plots)
%Explicit Euler in Imaginary time
%    Propagates solution according to
%     \psi(t+dt) = \psi(t) - dt*H\psi(t)

fprintf('\nRunning Explicit Euler in Imaginary time\n')

% Set psi old to current psi
psi_o  = psi;
psi_ex = psi; % Exact solution 
err = zeros(nsteps,4); % error vectors
psi2 = zeros(nsteps,4);
orthm = zeros(nsteps, 1);

% Movie setup crap
f = figure('Visible',plots.showm);
plot(psi(:,1:4));
axis tight manual;
ax = gca;
ax.NextPlot = 'replaceChildren';
f.Renderer = 'zbuffer';
Figs(nsteps) = struct('cdata',[],'colormap',[]);

%psi_o = psi_o.*exp(1i*meshgrid(diag(en))*dt);
X = meshgrid(diag(en));

for k = 1:nsteps
    %psi_o(:,1:4) = psi_o(:,1:4).*exp((1-1i)*X(:,1:4)*dt);
    psi_n = psi_o - dt*H*psi_o;
    psi_n(:,1:4) = psi_n(:,1:4).*exp((1-1i)*X(:,1:4)*dt);
    
    psi_o = psi_n;
    
    psi2(k,1:4) = dot(psi_o(:,1:4),psi_o(:,1:4));
    orthm(k) = dot(psi(:,1),psi(:,2))+dot(psi(:,1),psi(:,3))+...
        dot(psi(:,1),psi(:,4))+dot(psi(:,2),psi(:,3))+...
        dot(psi(:,2),psi(:,4))+dot(psi(:,3),psi(:,4));
    
    if mod(k,1) == 0
        psi_o = psi_o./meshgrid(sqrt(sum(abs(psi_o).^2,1)));
    end
    psi_ex(:,1:4) = psi_ex(:,1:4).*exp(-1i*X(:,1:4)*dt);
    
    err(k,1:4) = sqrt(sum(abs((psi_ex(:,1:4) - psi_n(:,1:4))).^2,1));

    if strcmp(plots.showm, 'on') || plots.savem
        plot(real(psi_o(:,1:4)));
        drawnow;
        Figs(k) = getframe;
    end
end
fprintf('Done. Ran for %i time steps\n',nsteps)

plot_err(err,psi2,plots.fn{1},plots,nsteps);

% Options for showing/saving movies
if strcmp(plots.showm,'on')
    movie(f,Figs,2,30)
if plots.savem
    fprintf('Saving movie data. Could take a while...\n')
    wObj = VideoWriter([plots.fn{1},'.avi']);
    open(wObj);
    for j = 1:nsteps
        writeVideo(wObj,Figs(j));
    end
    close(wObj);
    fprintf('Done writing %s\n',plots.fn{1})
end

end