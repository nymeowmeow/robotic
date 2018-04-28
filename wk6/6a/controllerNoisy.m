
function u = controllerNoisy(params, t, obs)
  % This is the starter file for the week5 assignment
  % Now you only receive noisy measurements for phi, and must use your EKF from week 3 to filter the data and get an estimate of the state
  % obs = [ay; az; gx] with a* in units of g's, and gx in units of rad/s

  % This template code calls the function EKFupdate that you must complete below
  xhat = EKFupdate(params, t, obs);
  phi = xhat(1);
  phidot = xhat(2);

  % The rest of this function should ideally be identical to your solution in week 4
  % Student completes this
  persistent e mytime;
  if isempty(e)
      e = 0;
      mytime = t;
  end
  
  dt = t - mytime;
  mytime = t;
  
  kp = 100;
  kd = 2;
  ki = 100;
  
  %e = e + (0 - phi*dt);
  %u = kp*(0 - phi) + kd*(0 - phidot) + ki*e;
  u = kp*phi + kd*phidot;
end

function xhatOut = EKFupdate(params, t, z)
  % z = [ay; az; gx] with a* in units of g's, and gx in units of rad/s
  % You can borrow most of your week 3 solution, but this must only implement a single predict-update step of the EKF
  % Recall (from assignment 5b) that you can use persistent variables to create/update any additional state that you need.

  % Student completes this
  % In this exercise, you will batch-process this data: you are provided a vector of timestamps (of length T), and a 3xT matrix of observations, z.
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