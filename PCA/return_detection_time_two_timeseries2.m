function [t_detect,t_detect_series] = return_detection_time_two_timeseries2(x_1,x_2,t,var_m)
%This function will use Betsy Weatherhead's approach to determine detection time to detection of change for a time-series
%This will contain adaptive slope estimation
%
% x_1: first time-series (e.g. predictions by model or 0 vector)
% x_2: second time-series (e.g. measurements or 0 vector)
% t : temporal coordinate
% var_m: measurement uncertainty (use zero assuming perfect data)
%
% Output:
% phi = ar(1) autocorrelation parameter
% sigma = white noise parameter
% t_detect_min: n* from Weatherhead (central estimate for amount of time
% required)
% t_detect_robust: n* at a point where parameters and t_detect are stable
% and record length is at least M
% t_detect_series: time series of n* as a function of record length, n*(M)

%check length
if length(x_1)~=length(x_2)
  error('time series length not equal');
end
if length(x_1)~=length(t)
  error('time series length not equal');
end

%Make sure that the time-series are oriented consistently (as column vectors)
a = size(x_1);
if a(1)<a(2)
  x_1 = x_1';
end
a = size(x_2);
if a(1)<a(2)
  x_2 = x_2';
end
a = size(t);
if a(1)<a(2)
  t = t';
end
t = t-t(1);

%find detection time by starting with an initial value and determining t_detect from data uncertainty

a = x_1;
b = x_2;
trend_vec = zeros(1,length(t));
phi_vec = zeros(1,length(t));
sigma_vec = zeros(1,length(t));
t_detect_series = zeros(1,length(t));
%start_index = 5;
start_index = 1;

c = zeros(size(b));
tp = zeros(size(b));

% vectorized version of the commented code that follows it:
anan = isnan(a);
bnan = isnan(b);
notNan = ~anan & ~bnan;
tp = t.*notNan;
c = (a - b).*notNan;
lastIndex = find(notNan, 1, 'last');
if min(size(lastIndex)) == 0 % if all numbers are NaNs, then lastIndex is a 0x1 matrix
    t_detect_min = nan;
    t_detect_robust = nan;
    t_detect_series = nan*ones(size(t));
    return
end
% c = c(1:lastIndex);
% tp = tp(1:lastIndex);

% index = 1;
% for ii=1:length(b)
%    if notNan(ii)
%       tp(index) = t(ii);
%       c(index) = a(ii)-b(ii);
%       cc(index) = a(ii)-b(ii);
%       index = index+1;
%    end
% end
% c = c(1:index-1);
% cc = cc(1:index-1);

%negligible time required for this step

for i=start_index:length(tp)
    X = [ones(i, 1) tp(1:i)];
    a = X\c(1:i);
    trend_vec(i) = a(2);
end

%[p,s] = polyfit(tp,cc,1);
%trend_vec(start_index:length(tp)) = p(1);
%  for i=start_index:length(tp)
%    [p,s] = polyfit(tp(1:i),cc(1:i),1); 
%     trend_vec(i) = p(1);
%  end

%This step takes about 2.9 seconds

t_detect_series = length(t);
index1 = 1;
index2 = 120; % look at 10 years of data (assuming monthly data)
xp = c(index1:index2)-mean(c(index1:index2));
xp = detrend(xp);
% dummy = size(xp);
% if (dummy(1)<dummy(2))
%   xp = xp';
% end

%[wp,phi_val,var_e,sbc,fpe,th]=arfit(xp,1,1);
% [~,phi_val,var_e,~,~,~]=arfit(xp,1,1);
[w,phi_val,var_e, sbc, fpe, th]=arfit(xp,1,1);
if abs(phi_val)<1
   t_vals = start_index:length(tp);
   B = 2*(4/3)./sqrt(t_vals)*sqrt((1+phi_val)/(1-phi_val)); % in order to take into account uncertainty on parameters
   t_mid = (3.96*sqrt(var_e)/(1-phi_val)./abs(trend_vec(t_vals))).^(2/3); % this is n* from Weatherhead
   t_detect_series = t_mid.*exp(B)*(1+var_m/var_e)^(1/3); % 3.96 not 3.3 since using 95% confidence not 90%
end

t_detect = t_detect_series(end);
end