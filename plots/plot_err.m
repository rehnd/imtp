function [] = plot_err(err,psi2,fname,plots,nsteps)
%Plots stuff for each propagation method

if strcmp(plots.showp,'on') || plots.savep
    fe = figure('Visible',plots.showp);
    semilogy(err);
    grid on;
    ylim([min(err(:,1)),10]); xlim([0,nsteps]);
    xlabel('Iteration Number')
    ylabel('|\psi_{num}(t) - \psi_{ex}(t)|')
    title('Error in Explicit Euler Imaginary Time')
    saveas(fe,['plots/',fname,'.png'])

    f2 = figure('Visible',plots.showp);
    semilogx(psi2);
    grid on;
    ylim([0.9,1.05])
    xlim([1,nsteps])
    ylabel('|\psi(t)|^2')
    xlabel('Iteration number')
    legend('n=1','n=2','n=3','n=4');
    saveas(f2,['plots/',fname,'_psi2.png']);

    % f3 = figure('Visible',plots.showp);
    % semilogx(orthm);
    % grid on;
    % ylim([0.9,1.05])
    % xlim([1,nsteps])
    % ylabel('|\psi(t)|^2')
    % xlabel('Iteration number')
end

end
