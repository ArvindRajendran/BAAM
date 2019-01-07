%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% University of Alberta, Edmonton, Canada
% Laboratory of Advanced Separation Processes (Prof. Dr. Arvind Rajendran)
%
% Project:  ........
% Year:     2019
% MATLAB:   R2018b, Windows 64bit
% Authors:  Vishal Subramanian Balashankar (VS)
%           Ashwin Kumar Rajagopalan (AK)
%          
%
% Purpose:
% ...
%
% Last modified:
% - 2019-01-06, AK: Introduced header for this file
%
% Input arguments:
% - 
%
% Output arguments:
% - 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [qA1, qA2, qA3, qA4, qA5, qA6,qB1, qB2, qB3, qB4, qB5, qB6, pspand,YALL1,qAfLPPd,qBfLPPd,qAffeed,qBffeed,YOUT]=BAAM(Yinput,Adsorbents)

tic

%----allocating input parameters----
% 
rho=Yinput(1,1); % Particle Density [kg/m3]
qsbCO2=Yinput(1,2);                                % [mol/kg]
b0CO2=Yinput(1,3);                                 % [m3/mol]
dubCO2=Yinput(1,4)*1000;                           % [J/mol]
qsdCO2=Yinput(1,5);                                % [mol/kg]
d0CO2=Yinput(1,6);                                 % [m3/mol]
dudCO2=Yinput(1,7)*1000;                           % [J/mol]
qsbN2=Yinput(1,8);                                 % [mol/kg]
b0N2=Yinput(1,9);                                  % [m3/mol]
dubN2=Yinput(1,10)*1000;                           % [J/mol]
qsdN2=Yinput(1,11);                                % [mol/kg]
d0N2=Yinput(1,12);                                 % [m3/mol]
dudN2=Yinput(1,13)*1000;                           % [J/mol]

%------------------------------


%----adsorbent and column properties----

w=1;                                               % Adsorbent weight [kg]
e=0.37;                                            % Void fraction [-]
V=(w/(rho*(1-e)));                                 % Column volume [m3]
R= 8.314*10^(-5);                                  % Value of gas constant [(m3 bar)/( K mol)]  [for all calculations]
RE= 8.314;                                         % Value of gas constant [J/ mol K]
T=298.15;                                          % Feed temperature [K]
gamma=1.4;                                         % Adiabatic constant
neff=0.72;                                         % Vacuum pump efficiency

%----------------------------------------

%--------------DSL Parameters------------

bCO2=b0CO2*exp(-dubCO2/(RE*T));
bN2=b0N2*exp(-dubN2/(RE*T));
dCO2=d0CO2*exp(-dudCO2/(RE*T));
dN2=d0N2*exp(-dudN2/(RE*T));

%-----------------------------------------

%------------operating conditions---------

Phigh=1;                                           % High Pressure [bar]
Plow=0.03;                                         % Low  Pressure [bar]

delp=0.0001;                                       % For ODE solver
pspan=(Phigh:-delp:Plow)';                         % For ODE solver

yfeed=0.15;                                        % Feed concentration

YOUT(1,1:9)=0;                                     % For Printing the simulation output
kinc=1;

%---------Pint and Plow Span-----------------

X=Plow:0.001:0.1;                                  % Plow Range
Y=(Plow + 0.01):0.01:(Phigh - 0.01);               % Pint Range 

[Plgrid, Pigrid]=meshgrid(X,Y);                    %grid of Plow and Pint

[r, c]=size(Plgrid);


for i= 1:r
    for j=1:c
        if Pigrid(i,j)<Plgrid(i,j) || eq(Pigrid(i,j),Plgrid(i,j))
            Pigrid(i,j)=NaN;
        end      
    end
end



[YALL1]=BLOWEVAC();



for i= 1:size(Plgrid,2)
    
    lowp=find(abs(YALL1(:,1)-Plgrid(1,i))<10^(-6));
    
    [yAfLPP,~,qAfLPP, qBfLPP,NLPP]=mass_LPP_VSB(Plgrid(1,i),YALL1(lowp,2));
    [Nfeed,Nwaste,~,qAffeed,qBffeed]=mass_feed_VSB(yAfLPP);
    
    for j=1:size(Plgrid,1)
        
        if ~(isnan(Pigrid(j,i)) ||  isnan(Plgrid(1,i)))
            
            if (Plgrid(1,i)==Plow)
                qAfLPPd=qAfLPP;
                qBfLPPd=qBfLPP;
            end
            
            
            inter=find(abs(YALL1(:,1)-Pigrid(j,i))<10^(-6));
            nAbd=YALL1(inter,6);
            nBbd=YALL1(inter,7);
            nAevc=YALL1(lowp,6)-YALL1(inter,6);
            nBevc=YALL1(lowp,7)-YALL1(inter,7);
            
            purity=(nAevc/(nAevc+nBevc))*100;
            recovery=(nAevc/(Nfeed*yfeed))*100;
            enerbd=(YALL1(inter,8)*2.77778e-7)/(nAevc*44*1e-6);
            enerevc=((YALL1(lowp,8)-YALL1(inter,8))*2.77778e-7)/(nAevc*44*1e-6);
            enerT=enerbd+enerevc;
            wc=(nAevc)/(V*(1-e));
            PurP(j,i)=purity;
            RecP(j,i)=recovery;
            EnerP(j,i)=enerT;
            WcP(j,i)=wc;
            
            ncyclein=Nfeed;
            ncycleout=(Nwaste-NLPP) + (nAbd+nBbd) + (nAevc+nBevc);
            nTotalerr=((ncyclein-ncycleout)/ncycleout)*100;
            YOUT(kinc,:)=[Plgrid(1,i) Pigrid(j,i) purity recovery enerbd enerevc enerT wc nTotalerr];
            kinc=kinc+1;
            
            %fprintf(fid,'%f  %s  %s  %s  %f  %s  %f  %s  %f  %s  %f  %s  %f  %s  %f  %s  %f  %s  %f  %s  %f\n',adsindex,',',Adsorbents,',',Plgrid(1,i),',',Pigrid(j,i),',',purity,',',recovery,',',enerbd,',',enerevc,',',enerT,',',wc,',',nTotalerr);
            
        end
    end
    
end

PurP(PurP==0)=NaN;
PurPmax=PurP(:,1);
RecP(RecP==0)=NaN;
RecPmax=RecP(:,1);
EnerP(EnerP==0)=NaN;
WcP(WcP==0)=NaN;


pspand=0:0.001:1;
pspand=pspand';

[qA1, qB1]=isotherm(pspand,0);
[qA2, qB2]=isotherm(pspand,0.15);
[qA3, qB3]=isotherm(pspand,0.30);
[qA4, qB4]=isotherm(pspand,0.45);
[qA5, qB5]=isotherm(pspand,0.75);
[qA6, qB6]=isotherm(pspand,1);



%fclose(fid);


YOUT(:,10)=sqrt(YOUT(:,3).^2 + YOUT(:,4).^2);

inddis=find(YOUT(:,10)<110.3 & YOUT(:,10)>110.2);
maxdis=max(YOUT(:,10));

if ~isempty(inddis)
    YOUTdis(1:size(inddis,1),:)=YOUT(inddis(1:end,1),:);
    [minenergy, indmine]=min(YOUTdis(:,7));
    Fullenergy=(minenergy*1.1446) + 66.528;
    wcminenergy=YOUTdis(indmine,8);
    
    fid = fopen('BAAM_Results.txt', 'at');
    fprintf(fid,'%s%s%f%s%f%s%f%s%f%s%f\n',Adsorbents,',',T,',',maxdis,',',minenergy,',',Fullenergy,',',wcminenergy);
    fclose(fid);
else
    fid = fopen('BAAM_error.txt', 'at');
    fprintf(fid,'%s%s%f%s%f%s%s%s%s\n',Adsorbents,',',T,',',maxdis,',','-',',','-');
    fclose(fid);
end










    function [YALLBLOWEVAC]=BLOWEVAC()
        

        options=odeset('RelTol',1e-6,'AbsTol',1e-6);
        
        [PT, YT]= ode23s(@odeBLOWEVAC,pspan,yfeed,options);
        
        
        [qABLOWEVAC,qBBLOWEVAC]=isotherm(PT,YT);
        
        Ndes=0;
        NAdes=0;
        NBdes=0;
        energydes=0;
        
        
        for ides=1:(size(PT,1)-1)
            
            %Initial
            fAides=(PT(ides)*YT(ides)*V*e)/(R*T);
            fBides=(PT(ides)*(1-YT(ides))*V*e)/(R*T);
            NAides=fAides+(qABLOWEVAC(ides)*w);
            NBides=fBides+(qBBLOWEVAC(ides)*w);
            
            %Final
            fAfdes=(PT(ides+1)*YT(ides+1)*V*e)/(R*T);
            fBfdes=(PT(ides+1)*(1-YT(ides+1))*V*e)/(R*T);
            NAfdes=fAfdes+(qABLOWEVAC(ides+1)*w);
            NBfdes=fBfdes+(qBBLOWEVAC(ides+1)*w);
            
            Nbd=((NAides+NBides)-(NAfdes+NBfdes));
            Ndes(ides+1,1)=Ndes(ides,1)+Nbd;
            
            NAdes(ides+1,1)=NAdes(ides,1) + (NAides - NAfdes);
            
            NBdes(ides+1,1)=NBdes(ides,1) + (NBides - NBfdes);
            
            energydes(ides+1,1)= energydes(ides,1) + ( (gamma/(gamma-1))*((RE*T)/neff)* Nbd* ((1/PT(ides+1))^( (gamma-1)/gamma ) - 1));
            
            
            
        end
        
        YALLBLOWEVAC=[PT YT qABLOWEVAC qBBLOWEVAC Ndes NAdes NBdes energydes];
        
        
        function yprime=odeBLOWEVAC(P,y)
            
            yprime(1,1)=(R*T*w*y(1,1)*(y(1,1) - 1)*(R^2*T^2*bCO2*qsbCO2 - R^2*T^2*bN2*qsbN2 + R^2*T^2*dCO2*qsdCO2 - R^2*T^2*dN2*qsdN2 + P^2*bCO2*dN2^2*qsbCO2 + P^2*bN2^2*dCO2*qsdCO2 - P^2*bN2*dN2^2*qsbN2 - P^2*bN2^2*dN2*qsdN2 - 2*P^2*bCO2*dN2^2*qsbCO2*y(1,1) - 2*P^2*bN2^2*dCO2*qsdCO2*y(1,1) + 2*P^2*bN2*dN2^2*qsbN2*y(1,1) + 2*P^2*bN2^2*dN2*qsdN2*y(1,1) + P^2*bCO2*dCO2^2*qsbCO2*y(1,1)^2 + P^2*bCO2*dN2^2*qsbCO2*y(1,1)^2 + P^2*bCO2^2*dCO2*qsdCO2*y(1,1)^2 + P^2*bN2^2*dCO2*qsdCO2*y(1,1)^2 - P^2*bN2*dCO2^2*qsbN2*y(1,1)^2 - P^2*bN2*dN2^2*qsbN2*y(1,1)^2 - P^2*bCO2^2*dN2*qsdN2*y(1,1)^2 - P^2*bN2^2*dN2*qsdN2*y(1,1)^2 + 2*P^2*bCO2*bN2*dCO2*qsdCO2*y(1,1) - 2*P^2*bCO2*bN2*dN2*qsdN2*y(1,1) + 2*P^2*bCO2*dCO2*dN2*qsbCO2*y(1,1) - 2*P^2*bN2*dCO2*dN2*qsbN2*y(1,1) - 2*P^2*bCO2*bN2*dCO2*qsdCO2*y(1,1)^2 + 2*P^2*bCO2*bN2*dN2*qsdN2*y(1,1)^2 - 2*P^2*bCO2*dCO2*dN2*qsbCO2*y(1,1)^2 + 2*P^2*bN2*dCO2*dN2*qsbN2*y(1,1)^2 + 2*P*R*T*bCO2*dN2*qsbCO2 + 2*P*R*T*bN2*dCO2*qsdCO2 - 2*P*R*T*bN2*dN2*qsbN2 - 2*P*R*T*bN2*dN2*qsdN2 + 2*P*R*T*bCO2*dCO2*qsbCO2*y(1,1) - 2*P*R*T*bCO2*dN2*qsbCO2*y(1,1) + 2*P*R*T*bCO2*dCO2*qsdCO2*y(1,1) - 2*P*R*T*bN2*dCO2*qsdCO2*y(1,1) - 2*P*R*T*bN2*dCO2*qsbN2*y(1,1) + 2*P*R*T*bN2*dN2*qsbN2*y(1,1) - 2*P*R*T*bCO2*dN2*qsdN2*y(1,1) + 2*P*R*T*bN2*dN2*qsdN2*y(1,1)))/((w*((P*bCO2*qsbCO2)/(R*T + P*bN2 + P*bCO2*y(1,1) - P*bN2*y(1,1)) + (P*dCO2*qsdCO2)/(R*T + P*dN2 + P*dCO2*y(1,1) - P*dN2*y(1,1)) - (P^2*bCO2*qsbCO2*y(1,1)*(bCO2 - bN2))/(R*T + P*bN2 + P*bCO2*y(1,1) - P*bN2*y(1,1))^2 - (P^2*dCO2*qsdCO2*y(1,1)*(dCO2 - dN2))/(R*T + P*dN2 + P*dCO2*y(1,1) - P*dN2*y(1,1))^2) - w*y(1,1)*((P*bCO2*qsbCO2)/(R*T + P*bN2 + P*bCO2*y(1,1) - P*bN2*y(1,1)) - (P*bN2*qsbN2)/(R*T + P*bN2 + P*bCO2*y(1,1) - P*bN2*y(1,1)) + (P*dCO2*qsdCO2)/(R*T + P*dN2 + P*dCO2*y(1,1) - P*dN2*y(1,1)) - (P*dN2*qsdN2)/(R*T + P*dN2 + P*dCO2*y(1,1) - P*dN2*y(1,1)) + (P^2*bN2*qsbN2*(bCO2 - bN2)*(y(1,1) - 1))/(R*T + P*bN2 + P*bCO2*y(1,1) - P*bN2*y(1,1))^2 + (P^2*dN2*qsdN2*(dCO2 - dN2)*(y(1,1) - 1))/(R*T + P*dN2 + P*dCO2*y(1,1) - P*dN2*y(1,1))^2 - (P^2*bCO2*qsbCO2*y(1,1)*(bCO2 - bN2))/(R*T + P*bN2 + P*bCO2*y(1,1) - P*bN2*y(1,1))^2 - (P^2*dCO2*qsdCO2*y(1,1)*(dCO2 - dN2))/(R*T + P*dN2 + P*dCO2*y(1,1) - P*dN2*y(1,1))^2) + (P*V*e)/(R*T))*(R*T + P*bN2 + P*bCO2*y(1,1) - P*bN2*y(1,1))^2*(R*T + P*dN2 + P*dCO2*y(1,1) - P*dN2*y(1,1))^2);
            
        end
        
        
    end

    function [qA,qB]=isotherm(p,yA)
        
        
        c1=(p.*yA)./(R*T);
        c2=(p.*(1-yA))./(R*T);
        
        
        %qA=((qASat*b1*p*yA)/(1 + b1*p*yA + b2*p*(1-yA)))
        
        %qA=((qASat*b1*c1)/(1 + b1*c1 + b2*c2)) + ((qAdSat*d1*c1)/(1 + d1*c1 + d2*c2));
        
        qA=((qsbCO2*bCO2*c1)./(1 + bCO2.*c1 + bN2.*c2)) + ((qsdCO2*dCO2*c1)./(1 + dCO2.*c1 + dN2.*c2));
        
        %qB=((qBSat*b2*p*(1-yA))/(1 + b1*p*yA + b2*p*(1-yA) ));
        
        
        %qB=((qBSat*b2*c2)/(1 + b1*c1 + b2*c2 ))+ ((qBdSat*d2*c2)/(1 + d1*c1 + d2*c2 ));
        
        qB=((qsbN2*bN2*c2)./(1 + bCO2.*c1 + bN2.*c2 ))+ ((qsdN2*dN2*c2)./(1 + dCO2.*c1 + dN2.*c2 ));
        
        
    end



    function [yAfLPP,mbalLPPerror,qAfLPP, qBfLPP,NLPP]=mass_LPP_VSB(Pl,yAfevc)
        
        
        
        %options=optimset('Algorithm','Levenberg-Marquardt');
        options = optimoptions('fsolve','Display','off');
        
        zLPP=fsolve( @(x) LPP_VSB(x,Pl,yAfevc),[0 0],options);
        
        %Initial
        
        fAiLPP=(Pl*yAfevc*V*e)/(R*T);
        fBiLPP=(Pl*(1-yAfevc)*V*e)/(R*T);
        
        [qAiLPP, qBiLPP]=isotherm(Pl,yAfevc);
        
        NAiLPP=fAiLPP+(qAiLPP*w);
        NBiLPP=fBiLPP+(qBiLPP*w);
        
        %Final
        fAfLPP=(Phigh*zLPP(1)*V*e)/(R*T);
        fBfLPP=(Phigh*(1-zLPP(1))*V*e)/(R*T);
        [qAfLPP, qBfLPP]=isotherm(Phigh,zLPP(1));
        
        NAfLPP=fAfLPP+(qAfLPP*w);
        NBfLPP=fBfLPP+(qBfLPP*w);
        
        %fprintf('\n***PRESSURISATION STEP***');
        yAfLPP=zLPP(1);
        NLPP=zLPP(2);
        
        NLPPinitial=NAiLPP+(NLPP*yAfLPP)+NBiLPP+(NLPP*(1-yAfLPP));
        NLPPfinal=NAfLPP+NBfLPP;
        mbalLPPerror=abs(((NLPPinitial-NLPPfinal)/NLPPinitial)*100);
        
    end

    function zLPP=LPP_VSB(xLPP,Pl,yAfec)
        
        
        %xLPP(1) -> yAfLPP
        %xLPP(2) -> NLPP
        
        %Initial
        
        fAiLPP=(Pl*yAfec*V*e)/(R*T);
        fBiLPP=(Pl*(1-yAfec)*V*e)/(R*T);
        [qAiLPP, qBiLPP]=isotherm(Pl,yAfec);
        
        NAiLPP=fAiLPP+(qAiLPP*w);
        NBiLPP=fBiLPP+(qBiLPP*w);
        
        %Final
        fAfLPP=(Phigh*xLPP(1)*V*e)/(R*T);
        fBfLPP=(Phigh*(1-xLPP(1))*V*e)/(R*T);
        [qAfLPP, qBfLPP]=isotherm(Phigh,xLPP(1));
        
        NAfLPP=fAfLPP+(qAfLPP*w);
        NBfLPP=fBfLPP+(qBfLPP*w);
        
        %Overall Balance
        
        zLPP(1,1)=(NAiLPP+NBiLPP)-(NAfLPP+NBfLPP)+xLPP(2);
        zLPP(2,1)=NAiLPP-NAfLPP+(xLPP(2)*xLPP(1));
        
    end


    function [Nfeed,Nwaste,mbalfeederror,qAffeed,qBffeed]=mass_feed_VSB(yAfLPP)
        
        
        
        %options=optimset('Algorithm','Levenberg-Marquardt');
        options = optimoptions('fsolve','Display','off');
        
        zfeed=fsolve(@(x) feed_VSB(x,yAfLPP),[500 500],options);
        
        fAifeed=(Phigh*yAfLPP*V*e)/(R*T);
        fBifeed=(Phigh*(1-yAfLPP)*V*e)/(R*T);
        [qAifeed, qBifeed]=isotherm(Phigh,yAfLPP);
        
        NAifeed=fAifeed+(qAifeed*w);
        NBifeed=fBifeed+(qBifeed*w);
        
        %Final
        
        fAffeed=(Phigh*yfeed*V*e)/(R*T);
        fBffeed=(Phigh*(1-yfeed)*V*e)/(R*T);
        [qAffeed, qBffeed]=isotherm(Phigh,yfeed);
        
        NAffeed=fAffeed+(qAffeed*w);
        NBffeed=fBffeed+(qBffeed*w);
        
        
        %fprintf('\n\n***Feed STEP***');
        Nfeed=zfeed(1);
        Nwaste=zfeed(2);
        
        Nfeedinitial=NAifeed + (Nfeed*yfeed) -(Nwaste*yAfLPP) + NBifeed + Nfeed*(1-yfeed)-Nwaste*(1-yAfLPP);
        Nfeedfinal=NAffeed+NBffeed;
        
        mbalfeederror=abs(((Nfeedinitial-Nfeedfinal)/Nfeedinitial)*100);
        
        
    end

    function zfeed=feed_VSB(xfeed,yAfLPP)
        
        
        
        %xfeed(1) -> Nfeed
        %xfeed(2) -> NWaste
        
        %Initial
        
        fAifeed=(Phigh*yAfLPP*V*e)/(R*T);
        fBifeed=(Phigh*(1-yAfLPP)*V*e)/(R*T);
        [qAifeed, qBifeed]=isotherm(Phigh,yAfLPP);
        
        NAiLPP=fAifeed+(qAifeed*w);
        NBiLPP=fBifeed+(qBifeed*w);
        
        %Final
        fAffeed=(Phigh*yfeed*V*e)/(R*T);
        fBfLPP=(Phigh*(1-yfeed)*V*e)/(R*T);
        [qAffeed, qBffeed]=isotherm(Phigh,yfeed);
        
        NAffeed=fAffeed+(qAffeed*w);
        NBffeed=fBfLPP+(qBffeed*w);
        
        %Overall Balance
        
        zfeed(1,1)=(NAiLPP+NBiLPP)-(NAffeed+NBffeed) +xfeed(1) -xfeed(2);
        zfeed(2,1)=NAiLPP-NAffeed+(xfeed(1)*yfeed)-(xfeed(2)*yAfLPP);
        
    end



toc
end