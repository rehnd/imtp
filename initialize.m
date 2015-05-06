function [psi, en, H] = initialize(vars, plots)
%Sets up H = T + V and gets psi, energies
%   Solves the time independent SE

% Get vars arguments
h = vars.h; n = vars.n; xr = vars.xr;

% Kinetic Energy operator
D = sparse(1:n,1:n,-2*ones(1,n),n,n);
E = sparse(2:n,1:n-1,ones(1,n-1),n,n);
T = -(1/h^2)*(E + D + E');

% Harmonic Oscillator potential
x = linspace(xr(1), xr(2),n);
V = sparse(1:n,1:n, x.^2,n,n);

H = T + V;
H_full = full(H);

[psi, en] = eig(H_full);

fprintf('Constructed and Diagonalized H\n')

if plots.eigsts
    nwfs = 5;
    envs = zeros(n,nwfs);
    for i=1:nwfs
        envs(:,i) = en(i,i);
    end
    % scaling factor for plotting
    scfactor = en(1,1)*sqrt(n)/2.2;
    
    f = figure('visible',plots.showp);
    plot(x,psi(:,1:nwfs)*scfactor+envs)
    hold on;
    plot(x,diag(V));
    ylim([0,max(max(envs)*1.1)])
    xlim(xr)
    xlabel('x')
    ylabel('Re[\psi(x,t=0)]')
    legend('n=1','n=2','n=3','n=4','n=5')
    
    if plots.savep
        fname = 'plots/eigsts.png';
        saveas(f,fname)
        fprintf('Plotted Eigenstates: saved to %s\n',fname);
    end
end
end

