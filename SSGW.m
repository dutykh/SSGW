function [zs,ws,PP] = SSGW(kd,kH2,N,tol)
% SSGW: Steady Surface Gravity Waves.
%       Computation of irrotational 2D periodic surface pure gravity waves 
%       of arbitrary length in arbitrary depth. 
%
% MANDATORY INPUT PARAMETERS:
% kd  = k*d   : relative depth (wavenumber "k" times mean water depth "d").
% kH2 = k*H/2 : steepness (half the total wave height "H" times the wavenumber "k").
%
% OPTIONAL INPUT PARAMETERS:
% N   : number of positive Fourier modes (default, N=2048).
% tol : tolerance (default, tol=1e-14).
%
% OUTPUT PARAMETERS:
% zs  = complex abscissas at the free surface (at the computational nodes).
% ws  = complex velocity at the free surface (at the computational nodes).
% PP  = Physical Parameters: PP(1)=depth, PP(2)=wavenumber, PP(3)=wavelenght, 
%       PP(4)=celerity c_e, PP(5)=celerity c_s, PP(6)=Bernoulli constant, 
%       PP(7)=crest height, PP(8)=trough height, PP(9)=impulse, 
%       PP(10)=potential energy, pp(11)=kinetic energy, PP(12)=radiation stress,
%       PP(13)=momentum flux, PP(14)=energy flux, PP(16)=group velocity.
%
% NOTE: The output quantities are dimensionless with the following scaling.
% In deep water:   rho = g = k = 1.
% In finite depth: rho = g = d = 1.
%
% EXAMPLE 1. To compute a wave of steepness kH2=0.3 in infinite depth:
% [zs,ws,PP]=SSGW(inf,0.3);
%
% EXAMPLE 2. To compute a cnoidal wave with height-over-depth=0.5 and 
% length-over-depth=100:
% Hd=0.5; Ld=100; kd=2*pi/Ld; kH2=pi*Hd/Ld; [zs,ws,PP]=SSGW(kd,kH2);
% 
% EXAMPLE 3. For steep and long waves, the default number of Fourier modes
% must be increased. For instance, in order to compute a cnoidal wave with 
% height-over-depth=0.7 and length-over-depth=10000:
% Hd=0.7; Ld=10000; kd=2*pi/Ld; kH2=pi*Hd/Ld; [zs,ws,PP]=SSGW(kd,kH2,2^19);
%
% The program works for all but the (almost) highest waves.
% Edit the m-file for more details.

% For details of the algorithm and the notations, read:
% Clamond, D. & Dutykh, D. 2017. Accurate fast computation of steady 
% two-dimensional surface gravity waves in arbitrary depth. Preprint.
% https://hal.archives-ouvertes.fr/hal-01465813/
%
% This m-file was written with the purpose of clarity. The notations closely 
% match those of the paper above.

% Authors: D. Clamond & D. Dutykh.
% Version: 2017-02-08.

%--------------------------------------------------------------------------

% Check input parameters.
if nargin<2                                             
   error('Two dimensionless scalar parameters must be provided.');
end
if kd<0 || imag(kd)~=0 || kH2<0 || imag(kH2)~=0
   error('Input scalar parameters kd and kH2 must be real and positive.');
end
if nargin<3
   N = 2048;
end
if nargin<4
   tol =1e-14;
end

% Determine depth and choose parameters.
if 1-tanh(kd) < tol                                                        % Deep water case.
   d   = inf;                                                                % Depth.
   k   = 1;                                                                  % Wavenumber.
   g   = 1;                                                                  % Acceleration due to gravity.
   lam = 1/k;                                                                % Characteristic wavelength lambda.
else                                                                       % Finite depth case.
   d   = 1;                                                                  % Depth.
   k   = kd/d;                                                               % Wavenumber.
   g   = 1;                                                                  % Acceleration due to gravity.
   lam = tanh(kd)/k;                                                         % Characteristic wavelength lambda.   
end
c02 = g*lam;                                                               % Linear phase velocity squared.
H   = 2*kH2/k;                                                             % Total wwave height.
L   = pi/k;                                                                % Half-length of the computational domain (with c_r=c_e).
dal = L/N;                                                                 % Delta alpha.
dk  = pi/L;                                                                % Delta k.

% Vectors.
va  = (0:2*N-1)'*dal;                                                      % Vector of abscissas in the conformal space.
vk  = [ 0:N-1 -N:-1 ]'*dk;                                                 % Vector of wavenumbers.

% Initial guess for the solution:
Ups = (H/2)*(1+cos(k*va));                                                 % Airy solution for Upsilon.
sig = 1;                                                                   % Parameter sigma.

% Commence Petviashvili's iterations.
err  = inf;                                                                % Enforce loop entry.
iter = 0;                                                                  % Iterations counter.
tic;                                                                       % Start clocking.
while (err > tol)

  % Compute sigma and delta.      
  mUps = mean(Ups);                                                        % << Upsilon >>.
  Ys   = Ups - mUps;                                                       % Y_s.
  if d==inf                                                                % Deep water.
     sig = 1;                                                                % sigma.
     CYs = real(ifft(abs(vk).*fft(Ys)));                                     % C{ Y_s }.
     mys = -Ys'*CYs/N/2;                                                     % << y_s >>.
  else                                                                     % Finite depth.
     C_hat  =  vk.*coth((sig*d)*vk);      C_hat(1)  = 1/(sig*d);             % Operator C in Fourier space.
     S2_hat = (vk.*csch((sig*d)*vk)).^2;  S2_hat(1) = 1/(sig*d)^2;           % Operator S^2 in Fourier space.
     Ys_hat  = fft(Ys);
     E  = mean(Ys.*real(ifft(C_hat.*Ys_hat))) + (sig-1)*d;                   % Equation for sigma.
     dE = d - d*mean(Ys.*real(ifft(S2_hat.*Ys_hat)));                        % Its derivative.
     sig = sig - E/dE;                                                       % Newton new sigma.
     mys = (sig-1)*d;                                                        % << y_s >>.
  end
  del = mys - mUps;                                                        % Parameter delta.
  C_hat  =  vk.*coth((sig*d)*vk);  C_hat(1) = 1/(sig*d);                   % Updated operator C in Fourier space.
 
  % Compute Bernoulli constant B.
  Ups2  = Ups.*Ups;                                                        % Upsilon^2.
  mUps2 = mean(Ups2);                                                      % << Upsilon^2 >>.
  CUps  = real(ifft(C_hat.*fft(Ups)));                                     % C{ Upsilon }.
  CUps2 = real(ifft(C_hat.*fft(Ups2)));                                    % C{ Upsilon^2 }.
  DCU   = CUps(N+1) -  CUps(1);                                            % C{ Upsilon }_trough - C{ Upsilon }_crest.
  DCU2  = CUps2(N+1) - CUps2(1);                                           % C{ Upsilon^2 }_trough - C{ Upsilon^2 }_crest.
  Bg    = 2*del - H/sig*(1+del/d+sig*CUps(1))/DCU + DCU2/DCU/2;            % B/g.
  
  % Define linear operators in Fourier space.
  Cinf_hat = abs(vk);  Cinf_hat(1) = 0;                                    % Operator C_inf.  
  CIC_hat  = tanh((sig*d)*abs(vk));                                        % Operator C_inf o C^{-1}.      
  if d==inf, CIC_hat(1) = 1; end                                           % Regularisation.
  L_hat    = (Bg-2*del)*Cinf_hat - ((1+del/d)/sig)*CIC_hat;                % Operator L.
  IL_hat   = 1./L_hat;  IL_hat(1) = 1;                                     % Operator L^-1.
 
  % Petviashvili's iteration.
  Ups_hat = fft(Ups);                                                      % Fourier transform of Upsilon.
  CUps_hat = C_hat.*Ups_hat;
  LUps = real(ifft(L_hat.*Ups_hat));                                       % L{Upsilon}.
  Ups2_hat = fft(Ups.*Ups);                                                % Fourier transform of Upsilon^2.
  NUps_hat = CIC_hat.*fft(Ups.*real(ifft(CUps_hat)));
  NUps_hat = NUps_hat + Cinf_hat.*Ups2_hat/2;                              % Nonlinear term in Fourier space.
  NUps = real(ifft(NUps_hat));                                             % N{ Upsilon }.
  S = (Ups'*LUps)/(Ups'*NUps);                                             % Weight.
  U = S*S*real(ifft(NUps_hat.*IL_hat));                                    % New Upsilon.
  U = H * ( U - U(N+1) ) / ( U(1) - U(N+1) );                              % Enforce mean value.
  
  % Update values.
  err = norm(U-Ups,inf);                                                   % Error measured by the L_inf norm.
  Ups = U;                                                                 % New Upsilon.
  iter = iter+1;
  
end
time = toc;

% Post processing.
IH_hat = -1i*coth(sig*d*vk);  IH_hat(1) = 0;                               % Inverse Hilbert transform.
Ys  = Ups - mean(Ups);
Ys_hat  = fft(Ys);
CYs = real(ifft(C_hat.*Ys_hat));
Xs  = real(ifft(IH_hat.*Ys_hat));
mys = -Ys'*CYs/N/2;
Zs  = Xs + 1i*Ys;
dZs = ifft(1i*vk.*fft(Zs));
zs  = va + 1i*mys + Zs;
dzs = 1 + dZs;
B   = g*Bg;
ce  = sum( (1+CYs)./abs(dzs).^2 )/2/N;
ce  = sqrt(B/ce);
cs  = sig*ce;
ws  = -ce./dzs;
a   = max(imag(zs));
b   = -min(imag(zs));

xs = [ real(zs(N+1:end))-2*pi/k ; real(zs(1:N)) ];
ys = [ imag(zs(N+1:end)) ; imag(zs(1:N)) ];

if d==inf 
    Bce2d = 0;
    IC = 1./abs(vk); IC(1) = 0;
else
    Bce2d = (B-ce^2)*d;
    IC = tanh(vk*sig*d)./vk; IC(1) = sig*d;                                % Inverse C-operator.
end
ydx  = real(dzs).*imag(zs);
intI = -ce*mys;                                                            % Impulse.
intV = mean(ydx.*imag(zs))*g/2;                                            % Potential energy.
intK = intI*ce/2;                                                          % Kinetic energy.
intSxx = 2*ce*intI - 2*intV + Bce2d;                                       % Radiation stress.
intS = intSxx - intV + g*d^2/2;                                            % Momentum flux.
intF = Bce2d*ce/2 + (B+ce^2)*intI/2 + (intK-2*intV)*ce;                    % Energy flux.
cg   = intF/(intK+intV);                                                   % Group velocity.
K1   = imag(zs)'*(imag(zs)/2-Bg)/2/N;
ICydx  = real(ifft(IC.*fft(ydx))); 
errfun = norm( imag(zs) - Bg + sqrt(Bg^2+2*K1-2*ICydx) ,inf);              % Residual.

% Output data.
PP(1)=d; PP(2)=k; PP(3)=H; PP(4)=ce; PP(5)=cs; PP(6)=B; PP(7)=a; PP(8)=b; 
PP(9)=intI; PP(10)=intV; pp(11)=intK; PP(12)=intSxx; PP(13)=intS;
PP(14)=intF; PP(15)=cg;

% Display results.
subplot(211)
plot(xs,ys,xs,0*ys,'k--','LineWidth',1.5)
if d==inf
   xlim([-pi pi]);
   ylim([-b a]*1.0333);
   xlabel('$k\ x$', 'interpreter', 'latex','FontSize',18);
   ylabel('$k\ \eta$', 'interpreter', 'latex','FontSize',18);
else
   xlim([xs(1) xs(end)]);
   ylim([-b a]*1.0333);
   xlabel('$x / d$', 'interpreter', 'latex','FontSize',18);
   ylabel('$\eta / d$', 'interpreter', 'latex','FontSize',18);
end
title('Free Surface', 'interpreter', 'latex','FontSize',24);
legend('Surface Elevation','Still Water Level');
set(gcf,'color','w');

subplot(212)
semilogy(0:2*N-1,abs(fft(imag(zs))),'LineWidth',1.5)
xlim([0 N-1]);
ylim([1e-17 max(abs(fft(imag(zs))))]*1.0333);
xlabel('Fourier modes', 'interpreter', 'latex','FontSize',18);
ylabel('$\log_{10} | \mathcal{F}\{\tilde{y}\} |$', 'interpreter', 'latex','FontSize',18);
title('Spectrum', 'interpreter', 'latex','FontSize',24);
set(gcf,'color','w');

% Print physical parameters.
fprintf('              %15.14f\n', []                     );
fprintf('NUMERICAL PARAMETERS.%15.14f\n', []              );
fprintf('Number of positive Fourier modes:     N = %9i\n', N);
fprintf('Tolerance:                          tol = %15.14e\n', tol);
fprintf('Number of iterations:              iter = %9i\n', iter);
fprintf('Residual:                           res = %9i\n', errfun);
fprintf('Iterations time (s)                time = %15.14f\n', time);
fprintf('    %15.14f\n', []                               );
fprintf('PHYSICAL PARAMETERS.%15.14f\n', []               );
fprintf('Mean depth:                           d = %15.14e\n', d);
fprintf('Acceleration due to gravity:          g = %15.14e\n', g);
fprintf('Wavelength:                      2*pi/k = %15.14e\n', 2*pi/k);
fprintf('     %15.14f\n', []                              );
fprintf('WAVE CHARACTERISTICS.%15.14f\n', []              );
fprintf('Wave height:                          H = %15.14f\n', H);
fprintf('Crest height (amplitude):             a = %15.14f\n', a);
fprintf('Trough height:                        b = %15.14f\n', b);
fprintf('Stokes first phase celerity:        c_e = %15.14f\n', ce);
fprintf('Stokes second phase celerity:       c_s = %15.14f\n', cs);
fprintf('Linear phase celerity:              c_0 = %15.14f\n', sqrt(g*lam));
fprintf('Bernoulli constant:                   B = %15.14f\n', B);
fprintf('     %15.14f\n', []                              );
fprintf('INTERGRAL QUANTITIES (in the frame of reference with zero circulation).%15.14f\n', []              );
fprintf('Impulse:                              I = %15.14f\n', intI);
fprintf('Potential energy:                     V = %15.14f\n', intV);
fprintf('Kinetic energy:                       K = %15.14f\n', intK);
fprintf('Radiation stress:                   Sxx = %15.14f\n', intSxx);
if d~=inf
fprintf('Momentum flux:                        S = %15.14f\n', intS);
end
fprintf('Energy flux:                          F = %15.14f\n', intF);
fprintf('Group celerity:                     c_g = %15.14f\n', cg);
