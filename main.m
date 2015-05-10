function [] = main()

[vars, plots] = setup();

[psi, en, H] = initialize(vars,plots);

% Freq components are (E_2-E_1)
% Use 1/10 of that value for now:
dt = 0.01/abs(en(2,2)-en(1,1));

nsteps = 1500;
lam = dt/vars.h^2;

fprintf('\nUsing dt = %g\nN steps = %i\ndt/h^2 = %g\n',dt,nsteps,lam);

%psi_ex_it = explicit_euler_it(H,en,psi,nsteps,dt,plots);
%psi_ex_rt = explicit_euler_rt(H,en,psi,nsteps,dt,plots);

%psi_ex_it_sup = explicit_euler_it_sup(H,en,psi,nsteps,dt,plots);
%psi_im_rt = implicit_euler_rt(H,en,psi,nsteps,dt,plots);
% 
% psi_cn_it = crank_nicolson_it(H,en,psi,nsteps,dt,plots);
%psi_cn_rt = crank_nicolson_rt(H,en,psi,nsteps,dt,plots);

psi_cn_rt_sup = crank_nicolson_rt_sup(H,en,psi,nsteps,dt,plots);
% 
%psi_sd_it = sec_order_diff_it(H,en,psi,nsteps,dt,plots);
%psi_sd_rt = sec_order_diff_rt(H,en,psi,nsteps,dt,plots);
end