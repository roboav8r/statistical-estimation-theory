function [A] = Afun(xbar,mu)
%  AFUN    Generate the linearized state matrix A for the orbit 
%          propagation problem.  
%
%  INPUTS
%  xbar         nx-by-1 state about which the nonlinear dynamics function
%               f will be linearized, where xdot = f(x,t).  The state
%               vector elements are :
%               x = [rs;vs]
%               with rs the 3-by-1 position of the SV expressed in ECI
%               coordinates (J2000), and vs is the 3-by-1 velocity of the SV
%               with respect to the ECI coordinate system and expressed in
%               ECI coordinates.
%                            
%  mu           Earth's gravitational parameter (G*Mearth) (m^3/s^2)
%
%  OUTPUTS
%  A            The nx-by-nx linearized state matrix A where 
%               A = dftilde/dxtilde_k(t)
%+------------------------------------------------------------------+
% References:
%
%
% Author:  Todd Humphreys
%+==================================================================+
  
  nx = size(xbar,1);  % total elements in state vector x
  A = zeros(nx,nx);
  rs = xbar(1:3,1);
  rsmag = norm(rs);
  rx = rs(1); ry = rs(2); rz = rs(3);
  r3 = rsmag^3;  % shorthand for magnitude^3
  r5 = rsmag^5;  % shorthand for magnitude^5
  A = [0 0 0 1 0 0;
      0 0 0 0 1 0;
      0 0 0 0 0 1;
      2*mu*rx/rsmag^4 0 0 0 0 0;
      0 2*mu*ry/rsmag^4 0 0 0 0;
      0 0 2*mu*rz/rsmag^4 0 0 0];