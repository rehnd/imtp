function [] = main()

[vars, plots] = setup();

[psi, en, H] = initialize(vars,plots);

% Freq components are (E_2-E_1)
% Use 1/10 of that value for now:
dt1 = 0.01/abs(en(2,2)-en(1,1))
dt2 = dt1;

%dt1 = vars.h^2/100000

nsteps = 1500;
dt1/vars.h^2

psi_ex_it = explicit_euler_it(H,en,psi,nsteps,dt1,plots);
psi_ex_rt = explicit_euler_rt(H,en,psi,nsteps,dt1,plots);

%psi_im_it = implicit_euler_it(H,en,psi,nsteps,dt1,plots);
%psi_im_rt = implicit_euler_rt(H,en,psi,nsteps,dt1,plots);
% 
% psi_cn_it = crank_nicolson_it(H,en,psi,nsteps,dt1,plots);
% psi_cn_rt = crank_nicolson_rt(H,en,psi,nsteps,dt,figname);
% 
% psi_sd_it = sec_order_diff_it(H,en,psi,nsteps,dt,figname);
% psi_sd_rt = sec_order_diff_rt(H,en,psi,nsteps,dt,figname);
end