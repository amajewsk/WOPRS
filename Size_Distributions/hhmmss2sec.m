function [seconds]=hhmmss2sec(hhmmss_time)
%return seconds from start of day
hours = floor(hhmmss_time/1d4);
minutes = floor(rem(hhmmss_time,1d4)/1d2);
seconds = rem(hhmmss_time,1d2);

seconds = hours*3.6d3 + minutes*6d1 + seconds;

end