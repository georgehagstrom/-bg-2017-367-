function CP = CP(I,T,P)

% This function calculates the C:P ratio of the maximally growing plankton
% type in the environment specified by the function arguments


% How to use: I is light level in microeinsteins per m^2 per second
% T is temperature in Celsius
% P is phosphate concentration in mol/L




rho=1.0;    % This is the density of water
pDryWeight=.47; % This is the percent of cell mass that is dry mass.

alphaS = .12;
gamma = .2;

    function mDryCell = mDryCell(pDryWeight, rho, r)
        mDryCell = 4*pi/3.0*pDryWeight*rho*r^3*1e6*1e-18; % This is mass in grams
    end



% This function determines the percent dry mass devoted to structure as a
% function of the cell radius

function S =  S(r)
        S= alphaS/r+gamma ;
end



% Here are costs of synthesis and maintenance of different cellular components
PhiS = .67;
PhiMT0 = .001;
T0= 25.0;
QPhi = 2.0;

    function PhiM =  PhiM(T)
        PhiM= PhiMT0*QPhi^((T-T0)/10);
    end

% Now for the Geider Version. 

% These are units of gC/gC/hour

QF=0.17;   % Nitrogen/Carbon ratio
k1 = 6.1e-4*QF*3600;
k2 = 1.4e-3*QF*3600;

aph = 11.6*QF; % m^2 /gC
phiM = 1.0e-6 ; % gC/mu mol photons

% Here F1 and F2 have as units grams of carbon. 

    function Pm = Pm(F1,F2)
        Pm= min([k1*F1,k2*F2]);
    end

    function Pa = Pa(F1,F2,I)
        Pmax = Pm(F1,F2);
        Pa= Pmax*(1-exp(-aph*phiM*F2*I*3600/Pmax));
    end

cc = aph*phiM*3600;

    function optF1alt = optF1alt(I,L)
        if L==0
            optF1alt = 0;
            return;
        end
            
         
            
        function dPaL =  dPaL(F1)
            if F1==0
                dPaL = k1;
                return;
            end
            
            dPaL = k1*(1-exp(-cc*I*(L-F1)/(k1*F1)))-I*cc*L/F1*exp(-cc*I*(L-F1)/(k1*F1));
        end     
        if dPaL(k2*L/(k1+k2))>0
            optF1alt = k2*L/(k1+k2);
            return
        end
        dPaLf = @(F) dPaL(F);
        optF1alt = fzero(dPaLf,[0,k2*L/(k1+k2)]);
    end
    
    

    function truePA = truePA(L,I)
        if L==0
            truePA = 0;
            return
        end
        F1 = optF1alt(I,L);
        truePA = Pa(F1,L-F1,I);
    end




% Now for the biosynthesis terms
kST0 = .168; % Biosynthesis efficiency at temperature T=T0
Q10k = 2.0;
    function fBio = fBio(E,T)
        fBio = E*kST0*Q10k^((T-T0)/10.0);
    end


% Now for uptake rate terms:


% Allometric Predictors:

    function KPAllometry = KPAllometry(r)
        Volume = 4.0/3.0*pi*r^3;
        KPAllometry = 10^(-1.4)*(Volume^.41)*1e-6;
    end

    function VMaxPAllometry = VMaxPAllometry(r)
        Volume = 4.0/3.0*pi*r^3;
        VMaxPAllometry = 10^(-9.1)*Volume*1e-6/24;
    end

    function uptakeRatePAllometry = uptakeRatePAllometry(r,P)
        uptakeRatePAllometry = P*VMaxPAllometry(r)/(KPAllometry(r)+P);
    end
    

% In order to calculate the growth rate of a strategy, we need to calculate the quota of said strategy


alphaE=.55;
PPhosphoLipid = .042*0;
PDNA = .095;
PRibosomeEu = .047;

NProtein = .16;
NRibosomeEu = .156;
NDNA = .16;
NPhospholipid = .016;

CProtein = .53;
CDNA = .36;
CPhospholipid = .65;
CRibosomeEu = .419;

CCarb = .4;
CLipid = .76;

% From this ratio we can calculate the stoichiometry (N:P) of the biosynthesis component. Essentially we use the rule that 
% N is 16% of dry mass, and that P is 4.2% by mass of the phospholipid pool, and that P is 4.7% by mass of the Ribosome pool.

    function totalP = totalP(E,r)
        totalP = ((.3*alphaS/r+.05)*PPhosphoLipid+alphaE*E*PRibosomeEu+.01*PDNA)/31;
    end

    function totalPM = totalPM(E,r)
        totalPM = (.1*PPhosphoLipid+alphaE*E*PRibosomeEu+.01*PDNA)/31;
    end


    function totalN = totalN(E,r)
        totalN = ((.3*alphaS/r+.05)*NPhospholipid+.7*alphaS*NProtein/r+(1-alphaE)*E*NProtein+alphaE*E*NRibosomeEu+(1-E-alphaS/r-gamma)*NProtein+.01*NDNA)/14.0;
    end

    function totalC= totalC(E,r)
        totalC = ((.3*alphaS/r+.05)*CPhospholipid+.7*alphaS*CProtein/r+.10*CLipid+.04*CCarb+(1-alphaE)*E*CProtein+alphaE*E*CRibosomeEu+(1-E-alphaS/r-gamma)*CProtein+.01*CDNA)/12.0;
    end





    function QP = QP(r,L)
        QP = mDryCell(pDryWeight,rho,r)*totalP(1-S(r)-L,r);
    end

    function muP = muP(r,L,P)
        muP = uptakeRatePAllometry(r,P)/QP(r,L);
    end



% We are in position to define the growth rate as a function of environmental conditions and strategy





    function muGeider = muGeider(r,E,I,T,P)
        L = 1-S(r)-E;
    
        BioSynthesisRate = min([(truePA(1-S(r)-E,I)-PhiM(T))/(1+PhiS), fBio(E,T), muP(r,L,P)]);
        if BioSynthesisRate>0
            muGeider = BioSynthesisRate;
            return;
        end
        muGeider = min([0,truePA(1-S(r)-E,100)-PhiM(T)]);
    end




% In this function we will calculate the optimal strategy 


function [rr,EE] = optimalStrategyGeider(I,T,P)
    % The first step will be to find E as a function of r
    CE = alphaE*PRibosomeEu/31.0;
    KST = kST0*Q10k^((T-T0)/10.0);
    a=CE*KST  ;
    CP = .3*alphaS*PPhosphoLipid/31.0;
    CS = (.05*PPhosphoLipid + .01*PDNA)/31.0;
    function b = b(r)
        b= KST/1.0*(CP/r+CS);
    end
    function c =  c(r)
        c =  -VMaxPAllometry(r)*P/(P+KPAllometry(r))*1/(mDryCell(pDryWeight,rho,r));
    end
    
    
    function E = E(r)
        E = (-b(r)+sqrt(b(r)^2-4*a*c(r)))/(2*a);
    end

    
    
    zeroFunction = @(r)fBio(E(r),T)-(truePA(1-S(r)-E(r),I)-PhiM(T))/(1+PhiS);
    
    
    rMinFunction = @(r) muP(r,0,P)-fBio(1-S(r)-E(r),10);
    

    if rMinFunction(alphaS/(1-gamma))*rMinFunction(1000)>0
        rr = nan;
        EE = nan;
        return 
    end
    
    rMin = fzero(rMinFunction,[alphaS/(1-gamma),1000]);
    
    if zeroFunction(alphaS/(1-gamma))*zeroFunction(1000)>0
        rr = nan;
        EE = nan;
        return
    end
    rOpt=fzero(zeroFunction,[alphaS/(1-gamma), 1000]);
    
    
    rr = rOpt;
    EE = E(rOpt);
    
end
    
    
            

    function optimalCPGeider = optimalCPGeider(I,T,P)
        [rr, EE]=optimalStrategyGeider(I,T,P);
    

        optimalCPGeider = totalC(EE,rr)/(totalPM(EE,rr));
        
   
    end
    CPStructure = optimalCPGeider(I,T,P);
    f = 2500.0;
    CP = 1/(1/CPStructure+f*P);
    
    
end


