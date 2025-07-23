function [hhmmss]=sec2hhmmss(seconds)
%return HHMMSS.MS formatted time from seconds since start of day
hours = floor(seconds/3.6d3);
minutes = floor((seconds-hours*3.6d3)/6d1);
seconds = seconds - (hours*3.6d3) - (minutes*6d1);

hhmmss = hours*1d4 + minutes*1d2 + seconds;

end