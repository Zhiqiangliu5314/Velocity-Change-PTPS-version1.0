function [y]=bpfilt(x,dt,lf,hf);
  nyq=0.5/dt;
  wn=[lf/nyq,hf/nyq];
  [b,a]=butter(2,wn);
  nx=min(size(x));
  for ix=1:nx
     y(ix,:)=filtfilt(b,a,x(ix,:));
  end