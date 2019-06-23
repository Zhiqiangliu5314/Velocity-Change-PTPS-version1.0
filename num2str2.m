%-------------------new num2str-----------------------
function str=num2str2(val,nb);
s0=num2str(val);
if length(s0)<nb
    str(1:nb-length(s0))='0';
    str(nb-length(s0)+1:nb)=s0;
else
    str(1:nb)=s0(1:nb);
end
