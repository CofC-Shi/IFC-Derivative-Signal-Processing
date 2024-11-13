%Linear Detrend

function dts = Linear_Detrend(sig,blocks)

s_len = length(sig);

%calculate block length
bl_len = round(s_len/blocks);

%calculate amount of points to average on ends of block
av_len = round(bl_len/10);

%create empty trend data
grid = zeros(1,s_len);

for k = 1:blocks

  %index based on block
  start = ((k-1)*bl_len)+1;
  fin = (k*bl_len);

  %make sure it goes to the end
  if(k == blocks)
      fin = s_len;
  end

  index = start:fin;

  %calculate the start and end values based on averaging 
  av_start = mean(sig(start:start+(av_len-1)));
  av_end = mean(sig(fin - (av_len+1):fin));

  %create the linear trend for the block
  grid(index) = linspace(av_start,av_end,(fin - start)+1);

end
%figure(1)
%plot(grid);

%subtract the trend from the signal
dts = sig - grid;

end
