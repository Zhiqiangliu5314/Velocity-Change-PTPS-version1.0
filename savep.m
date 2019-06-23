function savep(yp,staa,stab,hd,out)




for i=length(yp):-1:1
    if yp(i,2)==0
        yp(i,:)=[];
    end
end
% for i=size(yp,1):-1:1
%     if yp(i,1)==0 || yp(i,2)>0.01
%         yp(i,:)=[]
%     end

%calculate the theta of two stations
% x1=sta(nb).lo; y1=sta(nb).la;
% x2=sta(na).lo; y2=sta(na).la;
% x3=sta(na).lo; y3=1;
% theta=acosd(dot([x2-x1,y2-y1],[x3-x2,y3-y2])/(norm([x2-x1,y2-y1])*norm([x1-x2,y3-y2])));
% if x1-x2>0 && y1-y2>0
%     theta=theta
% elseif x1-x2>0 && y1-y2<0
%     theta=theta
% elseif x1-x2<0 && y1-y2<0
%     theta=360-theta
% elseif x1-x2<0 && y1-y2>0
%     theta=360-theta
% end

%print results
yp(:,4)=(yp(:,3)/2+yp(:,2));
yp(:,3)=(-yp(:,3)/2+yp(:,2));
yp(:,2)=yp(:,2);
yp2=yp;
yp=yp';
if ~isempty(yp)

fid=fopen(out,'a');
%fprintf(fid,'%s%s%s%s %f\n','EQ-ST=-',a,'-',b,theta);
fprintf(fid,'%s\n','//');
fprintf(fid,'%s%s%s%s \n','EQ-ST=-',staa,'-',stab);
fprintf(fid,'%f %f %f \n',hd(32),hd(33),hd(34));
fprintf(fid,'%f %f %f \n',hd(36),hd(37),hd(38));
fprintf(fid,'%s\n','Mode#=0');
fprintf(fid,'%s\n','DataQuality=grade1');
fprintf(fid,'%s\n','TimWin: ');
fprintf(fid,'%s\n','10 25');
fprintf(fid,'%s\n',['NPRD: Tmin: Tmax: (NSEG=16 MPRD=475)']);
fprintf(fid,'%d %d %d\n',size(yp2,1),yp2(1,1) ,yp2(size(yp2,1)));
fprintf(fid,'%s\n','Peroid VGroup:');
fprintf(fid,'%f %f %f %f\n',yp);
fclose(fid);
end
end