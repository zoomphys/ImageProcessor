function str = getTime()
% return formatted time string
t = round(clock);
str = [num2str(t(1),'%04d') '-' num2str(t(2),'%02d') '-' num2str(t(3),'%02d') ' ' num2str(t(4),'%02d') ':' num2str(t(5),'%02d') ':' num2str(t(6),'%02d')];
