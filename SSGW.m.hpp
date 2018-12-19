// Automatically translated using m2cpp 2.0 on 2018-12-18 14:49:01

#ifndef SSGW_M_HPP
#define SSGW_M_HPP

#include <iostream>
#include "matlab2cpp.h"
#include <cstdio>
#include <armadillo>
#include <cmath>
#include <vector>
#include <iomanip>
using namespace arma ;

/*

% cSSGW: Periodic Surface Gravity Waves (c++ version) (https://github.com/FiniteTsai/cSSGW)
% This c++ code is modified from the SSGW of Denys Dutykh by using matlab2cpp
% Some information from SSGW are provided in the following:

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
%--------------------------------------------------------------------------*/

void SSGW(double kd, double kH2, cx_vec& zs, cx_vec& ws, vec& PP, int N=2048, double tol=1.e-14)
{

  double k,g,lam,d,c02,H,L,dal,dk,sig,err,mUps,mys,E,dE,del,mUps2,DCU,Bg,S,a,b,ce,cg,cs,B,Bce2d,K1,errfun,intF,intI,intK,intS,intSxx,intV,DCU2;

  int iter;

  vec va,vk,Ups,Ys,CYs,C_hat,S2_hat,Ups2,CUps,CUps2,Cinf_hat,CIC_hat,L_hat,IL_hat,LUps,NUps,U,Xs,xs,ys,IC,ydx,ICydx;

  cx_vec Ys_hat,Ups_hat,CUps_hat,Ups2_hat,NUps_hat,IH_hat,Zs,dZs,dzs;

  double time;

  // % Determine depth and choose parameters.
  if (kd<0 || kH2<0)
  {
    std::cerr << "Input scalar parameters kd and kH2 must be real and positive." << std::endl ;
  }
  if (1.-tanh(kd)<tol) // % Deep water case.
  {
    d = datum::inf ; // % Depth.
    k = double(1) ; // % Wavenumber.
    g = double(1) ; // % Acceleration due to gravity.
    lam = 1/k ; // % Characteristic wavelength lambda.
  }
  else // % Finite depth case.
  {
    d = 1. ; // % Depth.
    k = kd/d ; // % Wavenumber.
    g = 1. ; // % Acceleration due to gravity.
    lam = tanh(kd)/k ; // % Characteristic wavelength lambda.
  }

  c02 = g*lam ; // % Linear phase velocity squared.
  H = 2*kH2/k ; // % Total wwave height.
  L = datum::pi/k ; // % Half-length of the computational domain (with c_r=c_e).
  dal = L/(double)(N); // % Delta alpha.
  dk = datum::pi/L ; //  % Delta k.

  // % Vectors.

  va = arma::trans((m2cpp::fspan(0, 1, 2*N-1)))*dal ; // % Vector of abscissas in the conformal space.
  vk = arma::trans(arma::join_rows(m2cpp::fspan(0, 1, N-1), m2cpp::fspan(-N, 1, -1)))*dk ; // % Vector of wavenumbers.

  // % Initial guess for the solution:

  Ups = (H/2.0)*(1.+arma::cos(k*va)) ; // % Airy solution for Upsilon.
  sig = 1. ; // % Parameter sigma.

  // % Commence Petviashvili's iterations.

  err = datum::inf ; // % Enforce loop entry.
  iter = 0 ; // % Iterations counter.
  m2cpp::tic() ;

  while (err>tol)
  {
  // % Compute sigma and delta.
    mUps = mean(Ups) ; // % << Upsilon >>.
    Ys = Ups-mUps ; // % Y_s.

    if (d==datum::inf) // % Deep water.
    {
     sig = 1. ; // % sigma.
     CYs = arma::real(arma::ifft<cx_vec>(abs(vk)%arma::fft(Ys))) ; // % C{ Y_s }.
     mys = -dot(Ys,CYs)/(double)(N)/2.0 ; // % << y_s >>.
    }
    else // % Finite depth.
    {
     C_hat = vk/tanh((sig*d)*vk) ; // % Operator C in Fourier space.
     C_hat(0) = 1/(sig*d) ;
     S2_hat = arma::square((vk/sinh((sig*d)*vk))) ; // % Operator S^2 in Fourier space.
     S2_hat(0) = 1/pow((sig*d), 2) ;
     Ys_hat = arma::fft(Ys) ;
     E = mean(Ys%arma::real(arma::ifft<cx_vec>(C_hat%Ys_hat)))+(sig-1)*d ; // % Equation for sigma.
     dE = d-d*mean(Ys%arma::real(arma::ifft<cx_vec>(S2_hat%Ys_hat))) ; // % Its derivative.
     sig = sig-E/dE ; // % Newton new sigma.
     mys = (sig-1)*d ; // % << y_s >>.
    }

    del = mys-mUps ; // % Parameter delta.
    C_hat = vk/tanh((sig*d)*vk) ; // % Updated operator C in Fourier space.
    C_hat(0) = 1./(sig*d) ;

    // % Compute Bernoulli constant B.
    Ups2 = Ups%Ups ; // % Upsilon^2.
    mUps2 = mean(Ups2) ; // % << Upsilon^2 >>.
    CUps = arma::real(arma::ifft<cx_vec>(C_hat%arma::fft(Ups))) ; // % C{ Upsilon }.
    CUps2 = arma::real(arma::ifft<cx_vec>(C_hat%arma::fft(Ups2))) ; // % C{ Upsilon^2 }.

    DCU = CUps(N)-CUps(0) ; // % C{ Upsilon }_trough - C{ Upsilon }_crest.
    DCU2 = CUps2(N)-CUps2(0) ; // % C{ Upsilon^2 }_trough - C{ Upsilon^2 }_crest.
    Bg = 2.*del-H/sig*(1.+del/d+sig*CUps(0))/DCU+DCU2/DCU/2.0 ; // % B/g.

    // % Define linear operators in Fourier space.
    Cinf_hat = abs(vk) ;
    Cinf_hat(0) = 0. ;
    CIC_hat = tanh((sig*d)*abs(vk)) ;

    if (d==datum::inf) //% Regularisation.
    {
      CIC_hat(0) = double(1) ;
    }

    L_hat = (Bg-2*del)*Cinf_hat-((1.+del/d)/sig)*CIC_hat ; // % Operator L.
    IL_hat = 1/L_hat ; // % Operator L^-1.
    IL_hat(0) = 1. ;

    // % Petviashvili's iteration.

    Ups_hat = arma::fft(Ups) ; // % Fourier transform of Upsilon.
    CUps_hat = C_hat%Ups_hat ;
    LUps = arma::real(arma::ifft<cx_vec>(L_hat%Ups_hat)) ; // % L{Upsilon}.
    Ups2_hat = arma::fft(Ups%Ups) ; // % Fourier transform of Upsilon^2.
    NUps_hat = CIC_hat%arma::fft(Ups%arma::real(arma::ifft<cx_vec>(CUps_hat))) ;
    NUps_hat = NUps_hat+Cinf_hat%Ups2_hat/2.0 ; // % Nonlinear term in Fourier space.
    NUps = arma::real(arma::ifft<cx_vec>(NUps_hat)) ; // % N{ Upsilon }.

    S = dot(Ups,LUps)/dot(Ups,NUps) ; // % Weight
    U = S*S*arma::real(arma::ifft<cx_vec>(NUps_hat%IL_hat)) ; // % New Upsilon.
    U = H*(U-U(N))/(U(0)-U(N)) ; // % Enforce mean value.

    // % Update values.
    err = norm(U-Ups, "inf") ; // % Error measured by the L_inf norm.
    Ups = U ; // % New Upsilon.
    iter = iter+1 ;
  }

  time = m2cpp::toc() ;

  // % Post processing.
  IH_hat = -cx_double(0, 1)/tanh(sig*d*vk) ; // % Inverse Hilbert transform.
  IH_hat(0) = double(0) ;
  Ys = Ups-mean(Ups) ;
  Ys_hat = arma::fft(Ys) ;
  CYs = arma::real(arma::ifft<cx_vec>(C_hat%Ys_hat)) ;

  Xs = arma::real(arma::ifft<cx_vec>(IH_hat%Ys_hat)) ;
  mys = -dot(Ys,CYs)/N/2.0 ;
  Zs = Xs+cx_double(0, 1.)* Ys ;
  dZs = arma::ifft<cx_vec>(cx_double(0, 1)*vk%arma::fft(Zs)) ;
  zs = va+cx_double(0, 1.)*mys+Zs;
  dzs = 1.+dZs ;

  B = g*Bg ;
  ce = arma::sum((1+CYs)/arma::square(arma::abs(dzs)))/2.0/N ;
  ce = sqrt(B/ce) ;
  cs = sig*ce ;
  ws = -ce/dzs ;
  a = arma::max(arma::imag(zs)) ;
  b = -arma::min(arma::imag(zs)) ;

  xs = join_cols(arma::real(zs(arma::span(N, zs.n_rows-1)))-2*datum::pi/k, arma::real(zs(arma::span(0, N-1)))) ;
  ys = join_cols(arma::imag(zs(arma::span(N, zs.n_rows-1))), arma::imag(zs(arma::span(0, N-1))));

  if (d==datum::inf)
  {
    Bce2d = 0. ;
    IC = 1./arma::abs(vk) ;
    IC(0) = 0. ;
  }
  else
  {
    Bce2d = (B-pow(ce, 2))*d ;
    IC = tanh(vk*sig*d)/vk ; // % Inverse C-operator.
    IC(0) = sig*d ;
  }

  ydx = arma::real(dzs)%arma::imag(zs) ;
  intI = -ce*mys ; // % Impulse.
  intV = mean(ydx%arma::imag(zs))*g/2.0 ; // % Potential energy.
  intK = intI*ce/2.0 ; // % Kinetic energy.
  intSxx = 2*ce*intI-2*intV+Bce2d ; // % Radiation stress.
  intS = intSxx-intV+g*pow(d, 2)/2.0 ; // % Momentum flux.
  intF = Bce2d*ce/2.0+(B+pow(ce, 2))*intI/2.0+(intK-2*intV)*ce ; // % Energy flux.
  cg = intF/(intK+intV) ; // % Group velocity.
  K1 = arma::dot(arma::imag(zs),arma::imag(zs)/2.0-Bg)/2.0/N ;

  ICydx = arma::real(arma::ifft<cx_vec>(IC%arma::fft(ydx))) ;
  errfun = norm(arma::imag(zs)- Bg + arma::sqrt( Bg*Bg+2.*K1-2.*ICydx ),"inf") ; // % Residual.

  PP=vec(15);
  PP(0) = d ;
  PP(1) = k ;
  PP(2) = H ;
  PP(3) = ce ;
  PP(4) = cs ;
  PP(5) = B ;
  PP(6) = a ;
  PP(7) = b ;
  PP(8) = intI ;
  PP(9) = intV ;
  PP(10) = intK ;
  PP(11) = intSxx ;
  PP(12) = intS ;
  PP(13) = intF ;
  PP(14) = cg ;

  // % Display results.

  /* to be done
%{
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
%} */

  // % Print physical parameters.

  std::printf("             \n") ;
  std::printf("NUMERICAL PARAMETERS.\n") ;
  std::printf("Number of positive Fourier modes:     N = %9i\n", N) ;
  std::printf("Tolerance:                          tol = %15.14e\n", tol) ;
  std::printf("Number of iterations:              iter = %9i\n", iter) ;
  std::printf("Residual:                           res = %15.14e\n", errfun) ;
  std::printf("Iterations time (s)                time = %15.14e\n", time) ;
  std::printf("    \n") ;
  std::printf("PHYSICAL PARAMETERS.\n") ;
  std::printf("Mean depth:                           d = %15.14e\n", d) ;
  std::printf("Acceleration due to gravity:          g = %15.14e\n", g) ;
  std::printf("Wavelength:                      2*pi/k = %15.14e\n", 2*datum::pi/k) ;
  std::printf("    \n") ;
  std::printf("WAVE CHARACTERISTICS.\n") ;
  std::printf("Wave height:                          H = %15.14f\n", H) ;
  std::printf("Crest height (amplitude):             a = %15.14f\n", a) ;
  std::printf("Trough height:                        b = %15.14f\n", b) ;
  std::printf("Stokes first phase celerity:        c_e = %15.14f\n", ce) ;
  std::printf("Stokes second phase celerity:       c_s = %15.14f\n", cs) ;
  std::printf("Linear phase celerity:              c_0 = %15.14f\n", sqrt(g*lam)) ;
  std::printf("Bernoulli constant:                   B = %15.14f\n", B) ;
  std::printf("    \n") ;
  std::printf("INTERGRAL QUANTITIES (in the frame of reference with zero circulation).\n") ;
  std::printf("Impulse:                              I = %15.14f\n", intI) ;
  std::printf("Potential energy:                     V = %15.14f\n", intV) ;
  std::printf("Kinetic energy:                       K = %15.14f\n", intK) ;
  std::printf("Radiation stress:                   Sxx = %15.14f\n", intSxx) ;
  if (d!=datum::inf)
  {
    std::printf("Momentum flux:                        S = %15.14f\n", intS) ;
  }
  std::printf("Energy flux:                          F = %15.14f\n", intF) ;
  std::printf("Group celerity:                     c_g = %15.14f\n", cg) ;
}
#endif
