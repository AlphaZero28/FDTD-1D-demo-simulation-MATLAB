clc;clear;

%% constant
c0 = 3*10^8;
meu = 1.256637061 * 10^-6;
epsilon = 8.854187817 * 10^-12;
mat_eps = [12.9]; % material epsilon
mat_meu = [1];
n=1;
nmax=sqrt(mat_eps*mat_meu);
nmin =1;
navg = (nmin+nmax)/2;

%% grid resolution
lambda=c0/(10^10*nmax);
NRES = 20;
dz1 = lambda/NRES;
dmin = 0.05;
NDRES = 1;
dz2 = dmin/NDRES;
dz3 = min(dz1,dz2);

%% initilizing dimension increment
dc = 0.04;
Nz=ceil(dc/dz3);
dz4 = dc/Nz;
dz=min(dz3,dz4);
dt = (0.5*nmin*dz)/c0;

z=linspace(0,1,Nz);

%% initializing feilds
Hx = zeros(1,Nz);
Ey = zeros(1,Nz);
H1=0;H2=0;H3=0;
E1=0;E2=0;E3=0;

%% simulation constat
nzsrc = ceil(Nz/2) -20;
MAKE_MOVIE = 1;
i = 1;
SHOW_PLOT=1;


%% source & iteration count
fmax = 10^10; 
tau = 0.4/fmax;

tprop = nmax*Nz*dz/c0;
Ts = 12*tau;
T = 5*tprop + Ts;
steps = round(T/dt);
STEPS = linspace(1*10^(-10),11*10^(-10),500);
STEPS =STEPS';
t = STEPS;

Hsource = -0.4*source(STEPS,tau);
Esource = 0.4*source(STEPS,tau);

Tsrc = Esource + Hsource;

len = length(Hsource);


%% initializing fourier kernel
Nfreq = 500;
freq = linspace(0,fmax,Nfreq);
Kn = exp(-i*2*pi*freq*dt);
EyR = zeros(1,Nfreq);
EyT = zeros(1,Nfreq);
Esrc = zeros(1,Nfreq);


material_shape = {[60:75]};
if SHOW_PLOT
    fig1 = figure('Name','FDTD Analysis','NumberTitle','off','color','w');
    figure(fig1);
    h1 = plot([1],[1],'LineWidth',1.4);
    h4 = get(h1,'Parent');
    set(h4,'Fontsize',12);
    hold on
    h2 = plot([1],[1],'LineWidth',1.4);
    title('Ey Hx Plot');
    h3 = text(3,2.7,'');
    xline(nzsrc,'r','Label','Source');
    xline(85,'k','Label','Transmittance Detector');
    xline(nzsrc-10,'b','Label','Reflectance Detector');
    uistack(h1,'top');
    uistack(h2,'top');
end


%% mateirals adding 
ER = ones(1,Nz);
UR = ones(1,Nz);
for i = 1:length(material_shape)
    m_shape = cell2mat(material_shape(i));
    ER(m_shape) = mat_eps(i);
    UR(m_shape) = mat_meu(i);
    if SHOW_PLOT
        hold on;
        fill([min(m_shape),m_shape,max(m_shape)],[-10,ones(1,length(m_shape))*10,-10],[0.5 0.5 0.5]);
        alpha(0.5);
        annotation('textbox', [0.61, 0.8, 0.1, 0.1], 'String', "Material")
    end
end

%% grid dispersion
k0 = (2*pi*fmax)/c0;
f = c0*dt*sin((k0*navg*dz)/2)/(navg*dz*sin((c0*k0*dt)/2))
UR = f*UR;
ER = f*ER;

%% update coefficients
mEy = (c0*dt) ./ (ER*dz);
mHx = (c0*dt) ./ (UR*dz);
Eyr=[];
Eyz=[];

%% main FDTD loop

for T = 1:len
    Es = Esource(T);
    Hs = Hsource(T);
    %% inject soft source: one directional
    Ey(nzsrc) = Ey(nzsrc) + Es;
    Hx(nzsrc) = Hx(nzsrc) + Hs;
    
    %% updating H from E 
    for nz = 1:Nz-1
        if nz ~= nzsrc-1
            Hx(nz) = Hx(nz) + mHx(nz)*(Ey(nz+1) - Ey(nz));
        else
            Hx(nz) = Hx(nz) + mHx(nz)*(Ey(nz+1) - Ey(nz)) - mHx(nz)*Es;
        end
    end 
    Hx(Nz) = Hx(Nz) + mHx(Nz)*( E3 - Ey(Nz)); 
    H3=H2;H2=H1;H1=Hx(1);
    
    %% updating E from H
    Ey(1) = Ey(1) + mEy(1)*(Hx(1) - H3);
    for nz = 2 : Nz
        if nz ~= nzsrc 
            Ey(nz) = Ey(nz) + mEy(nz)*(Hx(nz) - Hx(nz-1));
        else
            Ey(nz) = Ey(nz) + mEy(nz)*(Hx(nz) - Hx(nz-1)) - mEy(nz)*Hs;
        end  
    end
    E3=E2;E2=E1;E1=Ey(Nz);
    
    %% recording time domain value from the monitors 
    Eyr = [Eyr Ey(nzsrc-10)];
    Eyz = [Eyz Ey(Nz-5)];
    
    %% plot of electric field
    disp(T)
    remain = rem(i,10)
    if  remain== 0
        if SHOW_PLOT 
        set(h3,'String',sprintf('Step Num: %04d/%d',T,len));
        hold on;
        grid on;
        set(h1,'XData',[1:Nz],'YData',Ey);
        set(h2,'XData',[1:Nz],'YData',Hx);
        legend({'Ey','Hx'},'Location','SouthEast');
        xlabel('z direction');
        ylabel('Ey and Hx Fields');
        axis([1 Nz -3 3]);
        end
    end
    m(i) = getframe(gcf);
    i = i + 1;
    
end

%% completing Fourier transform on the time values from the monitors
EyR = (fftshift(fft(Eyr)));
EyT = (fftshift(fft(Eyz)));
Esrc =(fftshift(fft(Esource)));

%% reflectance and transmittance computation
Esrc=Esrc';
REF = abs(EyR./Esrc).^2;

TRN = abs(EyT./(Esrc)).^2;

CON = REF + TRN;

%% plotting
fmx = 0.5/dt;
faxis = linspace(-fmx,+fmx,Nfreq);

fig2 = figure('Color','w');
subplot(221)
p = plot(faxis,REF,'r',faxis,TRN,'b',faxis,CON,'k','LineWidth',1);
legend('reflectance','transmittance', 'conserv');

subplot(222);
plot(faxis,abs(Esrc));
title('source FFT');grid on;

subplot(223);
plot(faxis,abs(EyR));
title('reflection FFT');grid on;

subplot(224);
plot(faxis,abs(EyT));
title('tranns FFT');grid on;

if SHOW_PLOT
    hold off
end

%% making movies
if MAKE_MOVIE
    vidWriter = VideoWriter('SF and TF.avi');
    vidWriter.FrameRate = 40;
    set(vidWriter,'Quality',60);
    open(vidWriter);
    writeVideo(vidWriter,m);
    close(vidWriter);
end

%% source function
function src = source(l,tau)

t0 = 6*tau;
t = l;
src = exp([-((t-t0)/tau).^2]);
end
