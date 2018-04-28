
function u = controllerNoisyEnc(params, t, obs, th, dth)
  % This is the starter file for the week5 assignment
  % Now you only receive noisy measurements for theta, and must use your EKF from week 3 to filter the data and get an estimate of the state
  % obs = [ay; az; gx] (same as last week)
  % New for 6b: you also have access to params.traj(t)

  % Template code (same as last week)
  xhat = EKFupdate(params, t, obs);
  phi = xhat(1);
  phidot = xhat(2);

  r = params.r;
  x = r*(th + phi);
  dx = r*(dth + phidot);
  
  kp1 = 0.15;
  kd1 = 0.15;
  phides = kp1*(params.traj(t) - x) + kd1*(0 - dx);
  
  % fill this out
  kp = 0.1;
  kd = 0.02;
  u = kp * sin(phi - phides) + kd * phidot;
end

function xhatOut = EKFupdate(params, t, z)
  % z = [ay; az; gx] with a* in units of g's, and gx in units of rad/s
  % You can borrow most of your week 3 solution, but this must only implement a single predict-update step of the EKF
  % Recall (from assignment 5b) that you can use persistent variables to create/update any additional state that you need.

  % Student completes this
  persistent P tlast X;
  if isempty(P)
      P = 1.*eye(2);
      tlast = t;
      X = [0 0]';
  end
  
  Q = 0.1.*diag([1, 1]);
  R = (0.01).*diag([1, 1, 1]);
  % Student completes this
  %X = rad2deg(xhat(:,1));
  dt = t - tlast;
  tlast = t;
  A = [ 1 dt; 0 1];
  update = A*X;
  P = A*P*A' + Q;
  H = [ cos(update(1)) 0; -sin(update(1)) 0; 0 1];
  h = [ sin(update(1)), cos(update(1)), update(2) ]'; 
  K = P * H' * pinv(H*P*H' + R);
  X = update + K*(z - h);
  P = (eye(2) - K*H)*P;
  xhatOut = X;
end
