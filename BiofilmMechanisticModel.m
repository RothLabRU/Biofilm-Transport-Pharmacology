%Biofilm Mechanistic Model
function [Eaval] = BiofilmMechanisticModel(t,S,Dose)
global v0;
 
%Constant Parameters
umax = 7.5; %Maximum specific growth rate
b = 0.01; %Specific death rate
Yxs = 0.45; %Yield coefficient
Ks = 2.55; %Monod coefficient
Ec = 0.115; %Cell volume fraction
px = 30000; %Cell intrinsic density
LL = 105*10^-6; %Concentration boundary layer thickness
Csi = 2500; %Substrate influent concentration [8 for CDC; 2500 for Perfusion]
tb = 1; %Biocide dose duration
tDoseStart = 0;

%For CDC Bioreactor
% if t>=tDoseStart && t<=tDoseStart+tb
%     Cbi = Dose; %Biocide influent concentrartion
% else
%     Cbi = 0;
% end

%For Perfusion Bioreactor
thalf = 2; % Half-life [Variable]
kh = log(0.5)/(thalf/24);
if t>=tDoseStart && t<=tDoseStart+tb
    Cbi = Dose*exp(kh*(t-tDoseStart)); %Biocide influent concentrartion
else
    Cbi = 0;
end

%Parameters Continued
Ds = 5.7*10^-4; %Substrate diffusion coefficient
Db = 1.5*10^-4; %Biocide diffusion coefficient
tau = 0.8; %Biofilm/bulk diffusivity ratio
kd = 7200; %Detachment rate coefficient
AV = 32; %Surface-area-to-volume ratio [10 for CDC; 32 for Perfusion]
QV = 80; %Dilution rate [10 for CDC; 80 for Perfusion]
kb = 75; %Biocide disinfection rate coefficinet
kr = 0.0385*24; %Biocide reaction rate coefficient
kr50 = 5; %EC50 [Variable]
krMax = kr*30; %Emax [Variable] 
gamma = 1; %Gamma [Variable]
kt = 10; %kt [Variable]
kp0 = 0.024; %kp0 [Variable]
kz0 = 2; %kz0 [Variable]

%Acquiring Matrix Components
N = (length(S)+2)/5;
h = 1/(N-1);
Ea(2:N) = S(1:N-1);
Lfa = S(N);
Cs(2:N-1) = S(N+1:2*N-2);
Cb(2:N-1) = S(2*N-1:3*N-4);
v(1:N) = v0;
XaS = S(3*N-3);
XiS = S(3*N-2);
CsS = S(3*N-1);
CbS = S(3*N);
Et(2:N) = S(3*N+1:4*N-1);
Ep(2:N) = S(4*N:5*N-2);
Ea(1) = Ea(2);
Et(1) = Et(2);
Ep(1) = Ep(2);
Cs(1) = Cs(2);
Cb(1) = Cb(2);
Cs(N) = (CsS* h*Lfa + Cs(N-1) * LL*tau)/ ( LL*tau + h*Lfa);
Cb(N) = CbS;

%Define Ordinary & Partial Differential Equations
for q = 2:N-1
    kp(q) = kp0*(1 - (Cs(q)/(Ks + Cs(q))));
    kz(q) = kz0*(Cs(q)/(Ks + Cs(q)));
end
kp(1) = kp(2);
kp(N) = kp(N-1);
kz(1) = kz(2);
kz(N) = kz(N-1);

for j = 2:N
    v(j) = Lfa*h*umax/Ec*(Ea(j-1)*Cs(j-1)/(Ks+Cs(j-1))+Ea(j)*Cs(j)/(Ks+Cs(j)))/2 + v(j-1);
end

for l = 2:N-1
     dCbdt(l) = Db*tau*(((Cb(l+1) - 2*Cb(l) + Cb(l-1))/h^2)*(1/Lfa^2)) - kr*Cb(l)*Ea(l)*px ...
         + ((Cb(l) - Cb(l-1))/h)*(v(N) - kd*Lfa^2)*(((l-1)/(N-1))/Lfa);
end  

for i = 2:N-1
    dCsdt(i) = Ds*tau*(((Cs(i+1) - 2*Cs(i) + Cs(i-1))/h^2)*(1/Lfa^2)) - umax/Yxs...
       *(Cs(i)/(Ks+Cs(i)))*Ea(i)*px + ((Cs(i) - Cs(i-1))/h)*(v(N)...
       - kd*Lfa^2)*(((i-1)/(N-1))/Lfa);
end  

for k = 2:N
    dEadt(k) = ((umax*Cs(k))/(Ks+Cs(k)) - b - ((krMax*(Cb(k))^gamma)/(kr50^gamma + (Cb(k))^gamma)) - kp(k))*Ea(k) + kz(k)*Ep(k)...
        + (v(k)/(Lfa))*((((k-1)/(N-1))/v(k))*(v(N)...
        - kd*Lfa^2) - 1)*((Ea(k) - Ea(k-1))/h)...
        - ((umax*Cs(k))/(Ks+Cs(k)))*(Ea(k-1))^2/Ec;
end

for m = 2:N
    dEtdt(m) = ((krMax*(Cb(m))^gamma)/(kr50^gamma + (Cb(m))^gamma))*Ea(m) - kt*Et(m)...
        + (v(m)/(Lfa))*((((m-1)/(N-1))/v(m))*(v(N)...
        - kd*Lfa^2) - 1)*((Et(m) - Et(m-1))/h)...
        - ((umax*Cs(m))/(Ks+Cs(m)))*(Et(m-1))^2/Ec;
end

for p = 2:N
    dEpdt(p) = kp(p)*Ea(p) - kz(p)*Ep(p)...
        + (v(p)/(Lfa))*((((p-1)/(N-1))/v(p))*(v(N)...
        - kd*Lfa^2) - 1)*((Ep(p) - Ep(p-1))/h)...
        - ((umax*Cs(p))/(Ks+Cs(p)))*(Ep(p-1))^2/Ec;
end

dLfdt = v(N) - kd*Lfa^2;
dXaSdt = (umax*CsS)/(Ks+CsS)*XaS - b*XaS - kb*XaS*CbS ...
     + kd*Ea(N)*px*Lfa^2*AV - QV*XaS;
dXiSdt = b*XaS - kb*XaS*CbS + kd*(Ec-Ea(N))*px*Lfa^2*AV - QV*XiS;
dCsSdt = QV*(Csi-CsS) - (umax*CsS)/(Ks+CsS)*XaS/Yxs - Ds*tau*((Cs(N)-Cs(N-1))/h/Lfa)*AV;
dCbSdt = QV*(Cbi-CbS) - kr*Cb(N)*(XaS) - Db*tau*((Cb(N)-Cb(N-1))/h/Lfa)*AV;

%Extract Differential Equation Outputs & Assemble Matrix
Eaval = [dEadt(2:N)'; dLfdt; dCsdt(2:N-1)'; dCbdt(2:N-1)'; dXaSdt; dXiSdt; dCsSdt; dCbSdt; dEtdt(2:N)'; dEpdt(2:N)'];
v0 = v;
 end