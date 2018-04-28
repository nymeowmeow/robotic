
function xhat = EKFstudent(t, z)
  % In this exercise, you will batch-process this data: you are provided a vector of timestamps (of length T), and a 3xT matrix of observations, z.
  xhat = zeros(2,length(t));
  xhat(:,1) = [0, 0];
  
  P = 1.*eye(2);
  Q = 0.1.*diag([1, 1]);
  R = (0.01).*diag([1, 1, 1]);
  % Student completes this
  %X = rad2deg(xhat(:,1));
  X = xhat(:,1);
  for k = 2:length(t)
      myz = z(:,k);
      %myz(3) = deg2rad(myz(3));
      dt = t(k) - t(k-1);
      A = [ 1 dt; 0 1];
      %update = A*xhat(:,k-1);
      update = A*X;
      P = A*P*A' + Q;
      H = [ cosd(update(1))*pi/180 0; -sind(update(1))*pi/180 0; 0 1];
      h = [ sind(update(1)), cosd(update(1)), update(2) ]'; 
      K = P * H' * pinv(H*P*H' + R);
      %xhat(:,k) = update + K*(z(:,k) - h);
      X = update + K*(myz - h);
      P = (eye(2) - K*H)*P;
      %xhat(:, k) = rad2deg(X);
      xhat(:, k) = X;
  end 
end
