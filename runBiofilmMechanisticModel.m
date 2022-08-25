%run Biofilm Mechanistic Model

clear all;
close all;

global v0;

%Biocide Dose
DoseList = [1];
dLength = length(DoseList);

%Constant Parameters
Ec = 0.115; %Cell volume fraction
Lf = 5*10^-6; %Intial biofilm thickness
Csi = 2500; %Substrate influent concentration [8 for CDC; 2500 for Perfusion]
timeDays = 1; %Biocide dose duration
treatStart = 0;

N = 1000; %Number of Increments in Window

%Intializing
Ea0(1:N,1) = Ec;
Et0(1:N,1) = 0;
Ep0(1:N,1) = 0;
Cs0(1:N,1) = Csi;
Cb0(1:N,1) = 0; 
v0(1:N,1) = 0;
XaS = 0;
XiS = 0;
CsS = Csi;
CbS = 0;
 
%Solving
 for count = 1:dLength
 Dose = DoseList(count);
 tSpan = [0 timeDays];
 S0 = [Ea0(2:N);Lf;Cs0(2:N-1);Cb0(2:N-1);XaS;XiS;CsS;CbS;Et0(2:N);Ep0(2:N)];
 options = odeset('RelTol',1e-6,'AbsTol',1e-8);
 [tSol,SSol] = ode15s(@(t,S) BiofilmMechanisticModel(t,S,Dose), tSpan, S0, options);
     
%Plotting
if Dose == 0

figure
box on

subplot(2,2,1);
hold on
box on
plot(tSol,SSol(:,1:floor((N-1)/4):N-1));
xlabel('TIME (DAYS)','FontWeight','bold','FontAngle','italic');
ylabel('ACTIVE CELL VOLUME FRACTION','FontWeight','bold','FontAngle','italic');
title('Active Cell Volume Fraction');
legend('1-Bottom','2','3-Middle','4','5-Top');
xlim([0 timeDays]);

subplot(2,2,3);
box on
plot(tSol,SSol(:,N)*10^6);
xlabel('TIME (DAYS)','FontWeight','bold','FontAngle','italic');
ylabel('BIOFILM THICKNESS (\mum)','FontWeight','bold','FontAngle','italic');
title('Biofilm Thickness');
xlim([0 timeDays]);

subplot(2,2,4);
box on
plot(tSol,SSol(:,N+1:floor((N-1)/4)-1:2*N-2));
xlabel('TIME (DAYS)','FontWeight','bold','FontAngle','italic');
ylabel('SUBSTRATE CONCENTRATION (\mug/mL)','FontWeight','bold','FontAngle','italic');
title('Substrate Concentration');
legend('1-Bottom','2','3-Middle','4','5-Top')
xlim([0 timeDays]);

subplot(2,2,2);
box on
plot(tSol,SSol(:,4*N:floor((N-1)/4):5*N-2));
xlabel('TIME (DAYS)','FontWeight','bold','FontAngle','italic');
ylabel('PERSISTER CELL VOLUME FRACTION','FontWeight','bold','FontAngle','italic');
title('Persister Cell Volume Fraction');
legend('1-Bottom','2','3-Middle','4','5-Top')
xlim([0 timeDays]);

else
    
figure
box on

subplot(2,3,1);
box on
hold on
plot(tSol,SSol(:,1:floor((N-1)/4):N-1));
xlabel('TIME (DAYS)','FontWeight','bold','FontAngle','italic');
ylabel('ACTIVE CELL VOLUME FRACTION','FontWeight','bold','FontAngle','italic');
title('Active Cell Volume Fraction');
legend('1-Bottom','2','3-Middle','4','5-Top')
xlim([0 timeDays]);

subplot(2,3,4);
box on
plot(tSol,SSol(:,N)*10^6);
xlabel('TIME (DAYS)','FontWeight','bold','FontAngle','italic');
ylabel('BIOFILM THICKNESS (\mum)','FontWeight','bold','FontAngle','italic');
title('Biofilm Thickness');
xlim([0 timeDays]);

subplot(2,3,5);
box on
plot(tSol,SSol(:,N+1:floor((N-1)/4)-1:2*N-2));
xlabel('TIME (DAYS)','FontWeight','bold','FontAngle','italic');
ylabel('SUBSTRATE CONCENTRATION (\mug/mL)','FontWeight','bold','FontAngle','italic');
title('Substrate Concentration');
legend('1-Bottom','2','3-Middle','4','5-Top')
xlim([0 timeDays]);

subplot(2,3,3);
box on
plot(tSol,SSol(:,4*N:floor((N-1)/4):5*N-2));
xlabel('TIME (DAYS)','FontWeight','bold','FontAngle','italic');
ylabel('PERSISTER CELL VOLUME FRACTION','FontWeight','bold','FontAngle','italic');
title('Persister Cell Volume Fraction');
legend('1-Bottom','2','3-Middle','4','5-Top')
xlim([0 timeDays]);

subplot(2,3,2);
box on
hold on
plot(tSol,SSol(:,3*N+1:floor((N-1)/4):4*N-1));
xlabel('TIME (DAYS)','FontWeight','bold','FontAngle','italic');
ylabel('DAMAGED CELL VOLUME FRACTION','FontWeight','bold','FontAngle','italic');
title('Damaged Cell Volume Fraction');
legend('1-Bottom','2','3-Middle','4','5-Top')
xlim([0 timeDays]);

subplot(2,3,6);
box on
hold on
plot(tSol,SSol(:,2*N-1:floor((N-1)/4)-1:3*N-4));
xlabel('TIME (DAYS)','FontWeight','bold','FontAngle','italic');
ylabel('BIOCIDE CONCENTRATION (\mug/mL)','FontWeight','bold','FontAngle','italic');
title('Biocide Concentration');
legend('1-Bottom','2','3-Middle','4','5-Top')
xlim([0 timeDays]);
end

tLength = length(tSol);
for m = 1:tLength
    Xvb(m) = ((sum(SSol(m,1:N-1)))*SSol(m,N)/(Lf*Ec))/(N-1);
    Xvt(m) = ((sum(SSol(m,3*N+1:4*N-1)))*SSol(m,N)/(Lf*Ec))/(N-1);
    Xvp(m) = ((sum(SSol(m,4*N:5*N-2)))*SSol(m,N)/(Lf*Ec))/(N-1);
    Xv(m) = Xvb(m) + Xvt(m) + Xvp(m);
end
XvMax = 6/5*max(Xv);

figure
box on
patch([0 treatStart treatStart 0], [0 0 XvMax XvMax], 'w')
hold on
p2 = patch([treatStart timeDays timeDays treatStart], [0 0 XvMax XvMax], 'b');
plot(tSol,Xv,'-k')
hold off
set(p2,'FaceAlpha',0.3)
xlabel('TIME (DAYS)','FontWeight','bold','FontAngle','italic');
ylabel('FRACTION OF ORIGINAL','FontWeight','bold','FontAngle','italic');
title('Total Viable Cell Density');
xlim([0 timeDays]);
ylim([0 XvMax]);

plotValue1(count) = Xvb(tLength);
plotValue2(count) = Xvt(tLength);
plotValue3(count) = Xvp(tLength);
plotValue4(count) = Xv(tLength);
end

finalActive = plotValue1';
finalDamaged = plotValue2';
finalPersister = plotValue3';
finalTotal = plotValue4';