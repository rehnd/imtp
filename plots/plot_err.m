function [] = plot_err(err,psi2,pn,plots,nsteps)
%Plots stuff for each propagation method

if strcmp(plots.showp,'on') || plots.savep
    fe = figure('Visible',plots.showp);
    semilogy(err);
    grid on;
    ylim([min(err(:,1)),10]); xlim([0,nsteps]);
    xlabel('Iteration Number')
    ylabel('|\psi_{num}(t) - \psi_{ex}(t)|')
    title(['Error in ',plots.names{pn}])
    saveas(fe,['plots/',plots.fn{pn},'.png'])

    f2 = figure('Visible',plots.showp);
    plot(psi2);
    grid on;
    ylim([0.95,1.05])
    xlim([1,nsteps])
    ylabel('|\psi(t)|^2')
    xlabel('Iteration number')
    title([plots.names{pn}])
    legend('n=1','n=2','n=3','n=4');
    saveas(f2,['plots/',plots.fn{pn},'_psi2.png']);

    % f3 = figure('Visible',plots.showp);
    % semilogx(orthm);
    % grid on;
    % ylim([0.9,1.05])
    % xlim([1,nsteps])
    % ylabel('|\psi(t)|^2')
    % xlabel('Iteration number')
end

end
