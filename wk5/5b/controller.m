
function u = controller(params, t, phi, phidot)
  % STUDENT FILLS THIS OUT
  % 
  % Initialize any added state like this:
  % 
  % persistent newstate
  % if isempty(newstate)
  %   % initialize
  %   newstate = 0;
  % end
  % 
  persistent e;
  if isempty(e)
      e = 0;
  end
  persistent mytime;
  if isempty(mytime)
      mytime = t;
  end
  dt = t - mytime;
  mytime = t;
  kp = 30;
  kd = 1;
  ki = 800;
  
  e = e + (0 - phi*dt);
  u = kp*(0-phi) + kd*(0 - phidot) + ki*e;
  u = -u;
end

