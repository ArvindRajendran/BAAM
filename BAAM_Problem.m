clear all
clc

load('Adsorbents_info.mat');

[rowf, colf]=size(den_DSL);


for i=2  %1:4        %2- Zeolite-13X
 
   
    
[qA1, qA2, qA3, qA4, qA5, qA6,qB1, qB2, qB3, qB4, qB5, qB6, pspan,YALL1,qAfLPP,qBfLPP,qAffeed,qBffeed,YOUT]=BAAM(den_DSL(i,:),Adsorbents(i,1));

%%Isotherm plot

semilogx(pspan,qA1,'linewidth',2); 
xlabel('Total Pressure [bar]');
ylabel(' CO_2 Loading [mol/kg]');
xlim([0 1]);
set(gca,'linewidth',1,'FontSize',14,'FontWeight','Bold');
caption = sprintf('CO_{2} Isotherm');
title(caption);

hold on

semilogx(pspan,qA2,'linewidth',2); hold on
semilogx(pspan,qA3,'linewidth',2); hold on
semilogx(pspan,qA4,'linewidth',2); hold on
semilogx(pspan,qA5,'linewidth',2); hold on
semilogx(pspan,qA6,'linewidth',2); hold on

semilogx(YALL1(:,1),YALL1(:,3),'linewidth',3); hold on

semilogx([YALL1(end,1) YALL1(1,1)],[YALL1(end,3) qAfLPP],'linewidth',3); hold on
semilogx([YALL1(1,1) YALL1(1,1)],[qAfLPP YALL1(1,3)],'linewidth',3); hold on

hold on

end

