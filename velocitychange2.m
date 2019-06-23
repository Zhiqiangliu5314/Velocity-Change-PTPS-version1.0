function varargout = velocitychange2(varargin)
%% Program for measuring velocity change from direct surface wave

set(0,'defaultfigurecolor','w') 
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @velocitychange2_OpeningFcn, ...
                   'gui_OutputFcn',  @velocitychange2_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before velocitychange2 is made visible.
function velocitychange2_OpeningFcn(hObject, eventdata, handles, varargin)
clear global par;
global par;
handles.output = hObject;
guidata(hObject, handles);


%% read parmeters in the file of velocitychange2.ini
fini=fopen('velocitychange2.ini','r');
    leq=[];
    while isempty(leq)
        str=fgets(fini);
        leq=strfind(str,'//directory');
    end
    str=fgets(fini);
    leq=strfind(str,'=');
    while leq
        str2=str(leq(1)+1:length(str)-1);
         if strcmp(str(1:leq(1)-1),'Datrefer')
            par.frefer=strcat(str2);
         elseif strcmp(str(1:leq(1)-1),'lf')
             par.lf=str2num(str2);
         elseif strcmp(str(1:leq(1)-1),'hf')
             par.hf=str2num(str2);
         elseif strcmp(str(1:leq(1)-1),'Datstk')
            par.fstk=strcat(str2);     
         elseif strcmp(str(1:leq(1)-1),'Key')
            key=strcat(str2);     
         elseif strcmp(str(1:leq(1)-1),'reg')
            par.reg=str2num(str2);
         elseif strcmp(str(1:leq(1)-1),'time')
            par.time=str2num(str2);            
         elseif strcmp(str(1:leq(1)-1),'Datout')
            par.out=strcat(str2);
         elseif strcmp(str(1:leq(1)-1),'daylm')
            par.daylm=str2num(str2);     
         elseif strcmp(str(1:leq(1)-1),'dtime')
            par.dtime=str2num(str2);                 
         end  
        str=fgets(fini);
        leq=strfind(str,'=');
    end
fclose(fini);
par.ksta=0;
list=dir([par.frefer '*' key '*']);
ista=0;
for ii=1:length(list)
    stalist(ii)={list(ii).name};
end
par.stalist=stalist;

% --- Outputs from this function are returned to the command line.
function varargout = velocitychange2_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;
function edit1_Callback(hObject, eventdata, handles)
function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in list1.
function list1_Callback(hObject, eventdata, handles)
global par;
%clear value
par.xi=0;
par.x=[];
par.med='x'
par.day=[];
cla(handles.plot2,'reset')
cla(handles.plot3,'reset')
cla(handles.plot4,'reset')
cla(handles.plot5,'reset')
cla(handles.plot6,'reset')
cla(handles.plot7,'reset')
ehf=get(handles.edit2,'string');
elf=get(handles.edit4,'string');
if ~strcmp(ehf,'hf(s)')&&~strcmp(ehf,'lf(s)')
    par.lf=1/str2num(elf);
    par.hf=1/str2num(ehf);
end


set(handles.list1,'String',par.stalist);
ksta=get(hObject,'Value');
fout0= get(hObject,'String');
sta=fout0{get(hObject,'Value')};
if ksta~=par.ksta
    par.xi=0;
    par.x =[];
end
par.ksta=ksta;
par.sta=sta;

%read reference CCF
 refer=[par.frefer cell2mat(par.stalist(ksta))];
 [hd,hd2,hd3,datrf]=readsac0(refer);
 dt=str2num(num2str(hd(1)));
 par.dt=dt;
 name=cell2mat(par.stalist(ksta));
 dtn=0;
%dat1=0.5*(dat1+flipud(dat1));

 [datr1]=bpfilt(datrf(1:(length(datrf)+1)/2)',par.dt,par.lf,par.hf);
 [datr2]=bpfilt(datrf((length(datrf)+1)/2:end)',par.dt,par.lf,par.hf);
dat1=[datr1(1:end) datr2(2:end)];
 dat1=dat1-mean(dat1);

 
par.hd=hd;
par.hd2=hd2;
par.hd3=hd3;
par.drf=dat1';
%obtain mouse x 
set(gcf,'WindowButtonDownFcn',@ButttonDownFcn);
st=dir([par.fstk name(1:end-6) '\'  name(1:end-6) '*']);
if isempty(st)
    return
end

%%obatin start day and end day
if par.daylm(1)==0
par.std=str2num(st(1).name(end-6:end-4));
par.edd=str2num(st(length(st)).name(end-6:end-4))
else
par.std=par.daylm(2);
par.edd=par.daylm(3);
end

%read daily CCFs
for i=par.std:par.edd
     fdt=[par.fstk name(1:end-6) '\'  name(1:end-6) '.' num2str2(i,3) '.stk'];
if exist(fdt,'file')
    dtn=dtn+1;
    [hd0,hd02,hd03,dat01(i).dat]=readsac0(fdt);
else
    dat01(i).dat=[];
end
end

%%plot daily CCFs
subplot (handles.plot2);
for i=par.std:par.edd
if ~isempty(dat01(i).dat)
%reverse
%dat01(i).dat=0.5*(dat01(i).dat+flipud(dat01(i).dat));


 [datd1]=bpfilt(dat01(i).dat(1:(length(dat01(i).dat)+1)/2)',par.dt,par.lf,par.hf);
 [datd2]=bpfilt(dat01(i).dat((length(dat01(i).dat)+1)/2:end)',par.dt,par.lf,par.hf);
dat01(i).dat=[datd1(1:end) datd2(2:end)];
%[dat01(i).dat]=bpfilt(dat01(i).dat',par.dt,par.lf,par.hf);
dat01(i).dat=dat01(i).dat-mean(dat01(i).dat);
par.datst(i).dat=dat01(i).dat';
plot(-par.dtime:par.dt:par.dtime,dat01(i).dat/max(dat01(i).dat));
hold on;
end
end
%%theorical arrive time of the surface wave
Wstart=2.5; 
Wend=1.5;
line([hd(51)/Wend,hd(51)/Wend],[-1,1],'Color','g');
line([-hd(51)/Wend,-hd(51)/Wend],[-1,1],'Color','g');
line([hd(51)/Wstart,hd(51)/Wstart],[-1,1],'Color','g');
line([-hd(51)/Wstart,-hd(51)/Wstart],[-1,1],'Color','g');
set(handles.plot2,'xlim',[par.time(1) par.time(2)],'ylim',[-1 1]);
set(gca,'FontSize',15);
set(gca,'linewidth',1.2);


subplot(handles.plot6)
plot(par.hd(33),par.hd(32),'b^');
hold on;
plot(par.hd(37),par.hd(36),'r^');
plot(-155.6,19.5,'p');
plot(-155.5,19.7,'p');
plot(-155.25,19.4,'p');
line(handles.plot6,[par.hd(33),par.hd(37)],[par.hd(32),par.hd(36)]);
dis=round(111*distance(par.hd(33),par.hd(32),par.hd(37),par.hd(36)));
text(par.reg(2)-1,par.reg(4)-0.2,['dis=' num2str(dis) 'km'],'fontsize',15);
set(gca,'FontSize',12);
set(gca,'linewidth',1.2);
set(handles.plot6,'xlim',par.reg(1:2),'ylim',par.reg(3:4));

%% plot energy
subplot (handles.plot7);
for i=par.std:par.edd   
if ~isempty(par.datst(i).dat)
t=-par.dtime:par.dt:par.dtime;
y=par.datst(i).dat;
N=length(t)-1; %样点个数
%plot(t,y);
fs=100;%采样频率
df=fs/(N-1);%分辨率
f=(0:N-1)*df;%其中每点的频率
f2=1./f;
Y=fft(y(1:N))/N*2;%真实的幅值
%Y=fftshift(Y);
plot(f2(1:N/2),abs(Y(1:N/2)));
hold on
end
end
set(gca,'FontSize',12);
set(gca,'yscale','log') 
set(gca,'xlim',[1 10],'XTick',[2:2:10],'ylim',[100 100000]);
xlabel('period(s)')
ylabel('energy')


%%Function to obtain the x axis
function ButttonDownFcn(src,event)
global par

if strcmp(par.med,'x')
% if par.xi>=2
%     par.xi=0;
% end
par.xi=par.xi+1;
pt = get(gca,'CurrentPoint');
x = pt(1,1);
par.x(par.xi)=x;
line([x,x],[-1,1],'Color','r','linewidth',1);
%line([-x,-x],[-1,1],'Color','r');
text(x,-0.8,[par.med num2str(par.xi)],'fontsize',12);

elseif strcmp(par.med,'day')
  if ~isempty(par.day)
        delete(par.ax1);
        delete(par.ax2);
  end
pt = get(gca,'CurrentPoint');
day = round(pt(1,1));
par.day=day;
par.ax1=line([day,day],[-10,10],'Color','b','linewidth',1);
%line([-x,-x],[-1,1],'Color','r');
par.ax2=text(day,1,num2str(day),'fontsize',12);
end





%% Executes on button press in pushbutton1_pick
%% kernel coda for measuring velocity change 
function pushbutton1_Callback(hObject, eventdata, handles)
global par
if par.xi==0
    fprintf('please pick the peaks!!!!!')
    return
end
cla(handles.plot3,'reset')
cla(handles.plot4,'reset')
cla(handles.plot5,'reset')
xi1=par.xi;
x1=par.x;
par.med='day';

%%search  the peaks and trough of  refer data between x1 and x2
[pkrf(:,1),pkrf(:,2)]=findpeaks(abs(par.drf));
pkrf(:,2)=(pkrf(:,2).*par.dt)-par.dtime;
%%update x and xi
ii=0;
for i=1:length(pkrf)
    if pkrf(i,2)<x1(2) && pkrf(i,2)>x1(1)
        ii=ii+1;
        x(1,ii)=pkrf(i,2);
    end  
end
xi=length(x);
par.xi=xi;
par.x=x;


for ipk=1:xi
        [vp,pos]=min(abs(pkrf(:,2)-x(ipk)));
        xpkrf(ipk)=pkrf(pos,2);
end
        %%measure the shift of peaks of everyday
for i=par.std:par.edd
        %%load data
       if isempty(par.datst(i).dat)
              continue;
       end
              pkst=[];
              [pkst(:,1),pkst(:,2)]=findpeaks(abs(par.datst(i).dat));
       if isempty(pkst)
              continue;
       end
              pkst(:,2)=(pkst(:,2).*par.dt)-par.dtime;
       %%find the peaks of  measure data
       for ipk=1:xi
              [vp,pos]=min(abs(pkst(:,2)-x(ipk)));
              xpkdt(i,ipk)=pkst(pos,2);
       end
end
  
for i=1:xi
        shift(par.std:par.edd,i)=xpkdt(par.std:par.edd,i)-xpkrf(1,i);
        %shift(par.std:par.edd,i)=xpkdt(par.std:par.edd,i)-xpkdt(488,i);
end       
        
for iday=par.std:par.edd
        yaa=@(a)sum((a.*(xpkdt(iday,:))-shift(iday,:)).^2);
        [a(iday,1),fya(iday,1)]=fminbnd(yaa,-1,1);
        err=(shift(iday,:)-a(iday).*(xpkdt(iday,:))).^2;
        err=sum(err)/xi;
        err=sqrt(err/sum((xpkdt(iday,:).^2))); 
        a(iday,1)=-a(iday,1);
        a(iday,2)=err;
end

par.shift=shift;
par.a=a;


%%plot time-delay curves of the marked peaks and troughs 
subplot (handles.plot4);
plot(par.std:par.edd,100.*shift(par.std:par.edd,1:xi)./xpkdt(par.std:par.edd,1:xi));
set(handles.plot4,'xlim',[par.std par.edd],'ylim',[-1.5 1.5],'YTick',[-1.5:0.5:1.5]);
grid on;
set(gca,'FontSize',15);
set(gca,'linewidth',1.2);
ylabel('dt[i]/t[i](%)','FontSize',15);
xlabel('day','FontSize',15);

%%plot velocity change curve
subplot (handles.plot3);
errorbar(par.std:par.edd,a(par.std:par.edd,1)*100,a(par.std:par.edd,2)*100);
hold on;
plot(par.std:par.edd,a(par.std:par.edd,1)*100,'linewidth',2);
legend('errbar','dv/v');
plot([1 par.edd],[0,0],'k');
set(handles.plot3,'xlim',[par.std par.edd],'ylim',[-2 2],'YTick',[-2:0.5:2]);
 set(gca,'YDir','reverse');
set(gca,'FontSize',15);
set(gca,'linewidth',1.2);
grid on;
ylabel('dv/v(%)','FontSize',15);
xlabel('day','FontSize',15);





% --- Executes on button press in pushbutton4_Look.
function pushbutton4_Callback(hObject, eventdata, handles)
global par;
cla(handles.plot5,'reset');
cla(handles.plot7,'reset');

hd=par.hd;
subplot (handles.plot5);
plot(-par.dtime:par.dt:par.dtime,par.datst(par.day).dat/max(par.datst(par.day).dat),'r','linewidth',1);
hold on;
plot(-par.dtime:par.dt:par.dtime,par.drf/max(par.drf),'b','linewidth',1);
legend(['day=' num2str(par.day)],'rf');
xlim([-hd(51)/0.5 hd(51)/0.5]);
ylim([-1 1]);
plot(par.x,par.shift(par.day,:),'.r','markersize',15);
ax=-200:0.1:200;
plot(ax,-par.a(par.day,1)*ax,'k','linewidth',1);
xlim([par.time(1) par.time(2)]);
ylim([-1 1]);
grid on;
set(gca,'FontSize',15);
set(gca,'linewidth',1.2);
text(0.5*par.time(1),1.2,['Velocity change=' num2str(par.a(par.day,1)*100) '%±' num2str(par.a(par.day,2)*100) '%'],'FontSize',15)

%% plot energy
subplot (handles.plot7);


t=-par.dtime:par.dt:par.dtime;
y=par.drf;
N=length(t)-1; %样点个数
%plot(t,y);
fs=100;%采样频率
df=fs/(N-1);%分辨率
f=(0:N-1)*df;%其中每点的频率
f2=1./f;
Y=fft(y(1:N))/N*2;%真实的幅值
%Y=fftshift(Y);
plot(f2(1:N/2),abs(Y(1:N/2)),'r');
hold on

 i=par.day
if ~isempty(par.datst(i).dat)
t=-par.dtime:par.dt:par.dtime;
y=par.datst(i).dat;
N=length(t)-1; %样点个数
%plot(t,y);
fs=100;%采样频率
df=fs/(N-1);%分辨率
f=(0:N-1)*df;%其中每点的频率
f2=1./f;
Y=fft(y(1:N))/N*2;%真实的幅值
%Y=fftshift(Y);
plot(f2(1:N/2),abs(Y(1:N/2)),'b');
hold on
end


legend(['day=' num2str(par.day)],'rf');
set(gca,'FontSize',12);
set(gca,'yscale','log') 
set(gca,'xlim',[1 10],'XTick',[2:2:10],'ylim',[100 100000]);
xlabel('period(s)')
ylabel('energy')





% --- Executes on button press in pushbutton4_Save.
function pushbutton5_Callback(hObject, eventdata, handles)
global par;
staa=par.sta(1:4);
stab=par.sta(8:11);
yp(:,2:3)=par.a(:,1:2);
yp(:,1)=1:par.edd;
%save data
savep(yp,staa,stab,par.hd,par.out);
%save picture
a=getframe(gcf);
name=cell2mat(par.stalist(par.ksta));
%imwrite(a.cdata,[name(1:14) '.bmp'],'bmp');
if exist(['jpg'],'dir')==0
    mkdir(['jpg']);
end
set(gcf,'color','white','paperpositionmode','auto');
print(gcf,'-r300','-dbitmap',['jpg/' name(1:14) '.bmp'])
aresaved=[name(1:14)]
%print(gcf,'-r600','-depsc',[name(1:14) '.eps']);
%saveas(gcf,'exprimentLightBundles.bmp');
 




% --- Executes during object creation, after setting all properties.
function list1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
