function qdd = eom(params, th, phi, dth, dphi, u)
  % This is the starter file for the week5 assignment

  % Provided params are
  % params.g: gravitational constant
  % params.mr: mass of the "rod"
  % params.ir: rotational inertia of the rod
  % params.d: distance of rod CoM from the wheel axis
  % params.r: wheel radius

  % Provided states are:
  % th: wheel angle (relative to body)
  % phi: body pitch
  % dth, dphi: time-derivatives of above
  % u: torque applied at the wheel

  %qdd = [0;0];
  % THE STUDENT WILL FILL THIS OUT
  r  = params.r;
  g  = params.g;
  m  = params.mr;
  ir = params.ir;
  l  = params.d;
 
  A = [ m*r^2, m*r*(r + l*cos(phi)); ...
        m*r^2 + m*r*l*cos(phi), m*r^2 + 2*m*r*l*cos(phi) + m*l^2 + ir];
  B = [ u + m*r*l*sin(phi)*dphi^2; m*g*l*sin(phi) + m*r*l*sin(phi)*dphi^2];
  
  qdd = linsolve(A, B);
end