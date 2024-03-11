clear; clc;
% close all;
%% system parameters
% alpha = 1.7;
% beta = -0.002;
% F = 1.22;

alpha = 0;
beta = -0.0125;
F = 1.71;


gam = 1;

%% numerical settings
%transverse domain - distance
%longitudinal domain - time
N = 512*2;                  %transverse discretization
L = sqrt(8*pi^2/abs(beta)); %normalized mode circumference
Theta = L;                  %transverse domain size
naxis = (-N/2:N/2-1).';
dth = Theta/N;
theta = dth*naxis;          %distance domain
dw = 2*pi/Theta;
k = fftshift(dw*naxis);     %wavenumber domain
T = 64*pi;                  %propagation distance in time
dt = 0.001;                 %time step size

%intial condition 1
eta = L/2*pi;
A0 = sqrt(abs(beta/gam));
u0 = A0*sech(eta*theta);

%intial condition 2
% tau = 0.01;
% u0 = sqrt(1.2)*exp(-(theta/tau).^2);

%solution and time storage arrays
usave = u0;
tsave = 0;

%adaptive step sie parameters
err_upper = 1e-7;
err_lower = 1e-9;
dt_max = 0.001;

tic;

uf = fft(u0);
A = -(1+1i*alpha);
% LOf = (A+1i*beta/2*k.^2); %Linear operator in frequency domain
LOf = (A-1i*k.^2); %Linear operator in frequency domain
ELOf = exp(LOf*dt/2); %Exponential of linear op in freq domain with split step (for step size h)
ELOf2 = exp(LOf*dt/4); %Exponential of linear op in freq domain with split step (for step size h/2)
Linv_Ff = fft(F*ones(N,1))./LOf; 
NL = @(ut,h) exp(1i*gam*abs(ut).^2*h);
if ~all(LOf)
    Linv_Ff = zeros(N,1);
end

t = 0;
save_int = T/20;
t_nextsave = save_int;
while (t < T)
    uf1 = solve_single_step(uf,ELOf,Linv_Ff,NL,dt); %Full step
    uf2 = solve_single_step(uf,ELOf2,Linv_Ff,NL,dt/2); %First half step
    uf2 = solve_single_step(uf2,ELOf2,Linv_Ff,NL,dt/2); %Second half step
    rel_err = norm(uf2-uf1,inf)/norm(uf2,inf); %Checking relative error. Should be less than err_lower is step size is good enough
    
    % adjust the step size
    if  rel_err > err_upper
        dt = dt*2/3;
    else
        uf = (4*uf2-uf1)/3; %Richardson extrapolation
        t = t+dt;
        if rel_err < err_lower
            dt = dt*2;
            if dt>dt_max
                dt = dt_max;
            end
        end
    end
    
    if t > t_nextsave
        %disp([t;dt;rel_err]);
        fprintf("time: %.15f\nstep size: %.15f\nrelative error: %.15f\n\n",t,dt,rel_err);
        t_nextsave = t_nextsave + save_int;
        ut = ifft(uf);
        usave = [usave ut];
        tsave = [tsave t];
    end    
end

toc

figure
mesh(tsave,theta,abs(usave));

figure
plot(theta,abs(ifft(uf)).^2);

save temp.mat;
