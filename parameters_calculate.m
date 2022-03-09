%This function provides the semi-empiric parameters used to the criteria found in literature, assossiated to a
%specific high entropy alloy (HEA) given by the user. 
%It is expressed as:
%[parameters]=parameters_calculate(a_fraction,n_element,properties,enthalpy,imenthalpy,Tinput)
%
%Input: a_fraction -> Normalized atomic fraction of elements in the alloy
%       n_element -> Number of candidate constituent elements of the computed alloys
%       properties -> Physical and chemical properties of the elements in alloy
%       enthalpy -> Enthalpy of mixing of binary systems of the elements in alloy
%       imenthalpy -> Enthalpy of intermetallic compounds of the elements in alloy
%       Tinput -> Entered temperature
%Output:parameters -> Calculated parameter values
%
%The references of the data and equations are indicated at the end of the code
%
function [parameters]=parameters_calculate(a_fraction,n_element,properties,enthalpy,imenthalpy,Tinput)
    %variable that contains the average of the properties of the elements
    %in the alloy, according to their atomic fraction%
    average_properties=[];
    n_data=length(properties(:,1));
    for i=1:1:n_data
        average=0;
        for j=1:1:n_element
            average=average+a_fraction(j)*properties(i,j);
        end
        average_properties(i)=average;
    end
    R=8.314472; %Ideal gases constant (J/molK)
    
    %Extracts the averages necessary for the following calculations and
    %those that will be shown in results%
    r=average_properties(2); 
    tm=average_properties(3); 
    xp=average_properties(5); 
    xa=average_properties(6); 
    vec=average_properties(7); 
    ea=average_properties(8);
    AW=average_properties(9);
    rho=average_properties(10);
    Cp=average_properties(11);
    Thc=average_properties(12);
    
    D=properties(2,:)*2; %Saves each element atomic diameter
    
    %Thermodynamics and difference parameters calculation%
    D_r=0; %difference atomic radio variable
    D_xp=0; %difference of Pauling's electronegativity variable
    D_xa=0; %difference of Allen's electronegativity variable
    D_vec=0; %difference of VEC variable
    smix=0; %Entropy of mixing variable
    hmix=0; %Variable define as entalphy of mixing of AB%
    him=0; %Variable define as entalphy of AB for intermetallic%
    him_min=1000; %variable storing the lower entalphy of AB for intermetallic value
    
    %The following loop calculates the difference parameters before square
    %root and thermodynamics parameters
    for i=1:1:n_element
        D_r=D_r+a_fraction(i)*(r-properties(2,i))^2; %[1]
        D_xp=D_xp+a_fraction(i)*(xp-properties(5,i))^2; %[2]
        D_xa=D_xa+a_fraction(i)*(xa-properties(6,i))^2; %[3]
        D_vec=D_vec+a_fraction(i)*(vec-properties(7,i))^2;
        
        if a_fraction(i)~=0 %If atomic fration is different of cero%
            smix=smix+a_fraction(i)*log(a_fraction(i)); %Entropy of mixing calculation[4]
        end
        for j=1:1:n_element 
            %One quarter of intermetallic compounds enthalpy of mixing calculation%
            hmix_ij=a_fraction(i)*a_fraction(j)*enthalpy(i,j);
            hmix=hmix+hmix_ij; %[1]
            him_=imenthalpy(i,j);
            him_ij=a_fraction(i)*a_fraction(j)*imenthalpy(i,j);
            him=him+him_ij; %[5]
            if him_<him_min
                him_min=him_;
            end
        end
    end
    
    hmix=4*hmix; %kJ/mol
    smix=-R*smix; %J/molK
    tsmix=Tinput*smix/1000; %kJ/mol
    gmix=round(hmix-tsmix,2); %kJ/mol
    D_r=D_r^0.5; %-
    d_r=D_r/r*100; %%
    D_xp=D_xp^0.5; %-
    d_xp=D_xp/xp*100; %%
    D_xa=D_xa^0.5; %-
    d_xa=D_xa/xa*100; %%
    D_vec=D_vec^0.5; %-
    him=round(4*him,2); %kJ/mol
     
    %Omega parameter calculation - Studies the solid solution stability [6]%
    omega=round(tm*smix/abs(1000*hmix),2); %output
    if omega>100
        omega=100;
    end
    
    %Eta paraneter calculation - Studies the competition between phases Im
    %and SS[7]
    eta=round(tsmix/abs(him_min),3); %output
    
    %Gamma parameter calculation - Lattice distortion due to the atomic
    %radii difference[8]
    rs=1000;
    rl=0;
    for i=1:1:n_element
        if a_fraction(i)~=0
            if rs>properties(2,i)
                rs=properties(2,i);%smaller radii variable%
            end
            if rl<properties(2,i)
                rl=properties(2,i);%larger radii variable%
            end
        end
    end
    
    ws=1-(((rs+r)^2-r^2)/((rs+r)^2))^0.5; %Gamma parameter numerator%
    wl=1-(((rl+r)^2-r^2)/((rl+r)^2))^0.5; %Gamma parameter denominator%
    gamma=round(ws/wl,5); %output
    
    %Lambda parameter calculation - Geometric origin parameter to determine
    %distorted solid solution (DCC)[9]
    lambda=round(smix/(d_r)^2,3);
    
    %kcr parameter calculation - Critical value to predict the absence of
    %intermetallic phase [10]
    k1=him/hmix;
    k2=0.6;
    kcr=round(-tsmix/(hmix)*(1-k2)+1,2);
    D_k=round(kcr-k1,2); %output
    
    sim=k2*smix; %J/molK
    tsim=Tinput*sim/1000; %kJ/mol
    gim=round(him-tsim,2); %kJ/mol gibbs free energy of intermetallic compounds
    
    %Excess entropy calculation [11][12]
    xit=0; %Xi total variable
    for i=1:1:n_element
        xi(i)=(D(i)^3)*a_fraction(i); %Matrix containing each Xi value for the alloy elements%
        xit=xit+xi(i); %Total Xi value actualization - adition of every Xi value for the alloy entered%
    end
    y1=0; %y1 variable%
    y2=0; %y2 variable%
    y3=0; %y3 variable%
    for i=1:1:n_element
        j=i+1;
        while j<=n_element
            deltaij=((D(i)-D(j))^2)*((xi(i)*xi(j)*a_fraction(i)*a_fraction(j))^0.5)/(xit*D(i)*D(j)); %Delta ij variable
            y1=y1+deltaij*(D(i)+D(j))*(D(i)*D(j))^-0.5; %function Y1 calculation%
            for k=1:1:n_element
                y2=y2+deltaij*(xi(k)*(D(i)*D(j))^0.5)/(xit*D(k)); %function Y2 calculation%
            end
            j=j+1;
        end
        y3=y3+((xi(i)/xit)^(2/3))*a_fraction(i)^(1/3); %Cube root calculation for function Y3%   
    end
    y3=y3^3; %Y3 function output used to calculate the excess entropy%
    
    xibcc=0.68; %BCC phases packing factor [13] 
    Z_bcc=((1+xibcc+xibcc*xibcc)-3*xibcc*(y1+y2*xibcc)-y3*xibcc^3)*(1-xibcc)^-3; %F2 function calculation%
    t1_bcc=-1.5*(1-y1+y2+y3)+(3*y2+2*y3)*(1-xibcc)^-1+1.5*(1-y1-y2-y3/3)*(1-xibcc)^-2+(y3-1)*log(1-xibcc); %F1 function calculation%
    se_bcc=R*(t1_bcc-log(Z_bcc)-(3-2*xibcc)*((1-xibcc)^-2)+3+log((1+xibcc+xibcc^2-xibcc^3)*(1-xibcc)^-3)); %Excess entropy output [14]
    gmix_bcc=round(hmix-Tinput*(smix+se_bcc)/1000,2);
    
    xifcc=0.74; %FCC phases packing factor [13] 
    Z_fcc=((1+xifcc+xifcc*xifcc)-3*xifcc*(y1+y2*xifcc)-y3*xifcc^3)*(1-xifcc)^-3; %F2 function calculation%
    t1_fcc=-1.5*(1-y1+y2+y3)+(3*y2+2*y3)*(1-xifcc)^-1+1.5*(1-y1-y2-y3/3)*(1-xifcc)^-2+(y3-1)*log(1-xifcc); %F1 function calculation%
    se_fcc=R*(t1_fcc-log(Z_fcc)-(3-2*xifcc)*((1-xifcc)^-2)+3+log((1+xifcc+xifcc^2-xifcc^3)*(1-xifcc)^-3)); %Excess entropy output [14]
    gmix_fcc=round(hmix-Tinput*(smix+se_fcc)/1000,2);
    
    %Phi parameter calculation - differentiation of single-phase or complex microstructures [15]%
    phi_bcc=round((smix-abs(1000*hmix)/tm)/abs(se_bcc),3); %output
    phi_fcc=round((smix-abs(1000*hmix)/tm)/abs(se_fcc),3); %output
    
    %PSFE parameter calculation - atomic percent of pair sigma-forming
    %element[16]
    %9x16 - Matrix of atomic number for components of binary AB sigma phases [17]
    %the first column is the atomic number of component A: [Al; V; Cr;
    %Mn; Mn; Fe; Nb; Mo; Ta; W]
    %the next column is the atomic number for compoenents B, if it corresponds :
    %[component A, Mn, Fe, Co, Ni, Nb, Tc, Ru, Rh, Pd, Ta, Re, Os, Ir]
    sfe=[13	0	0	0	0	41	0	0	0	0	73	0	0	0	0	0
        23	25	26	27	28	0	0	0	0	0	0	75	0	0	0	0
        24	25	26	27	0	0	0	0	0	0	0	0	0	0	0	0
        25	0	0	0	0	0	43	0	0	0	0	75	0	0	0	0
        26	0	0	0	0	0	43	0	0	0	0	75	0	0	0	0
        41	0	26	0	0	0	0	0	45	46	0	75	76	0	78	0
        42	25	26	27	0	0	43	44	0	0	0	75	76	77	0	0
        73	0	0	0	0	0	0	0	45	0	0	75	76	77	78	79
        74	0	0	0	0	0	43	44	0	0	0	75	76	77	0	0];
    psfe=0;
    a_fracc=a_fraction;
    for i=1:1:n_element
        A=0;
        for j=1:1:length(sfe(:,1))
            if sfe(j,1)==properties(1,i)
                A=a_fracc(i);
                for ii=1:1:n_element
                    B=0;
                    for jj=2:1:length(sfe(1,:))
                        if sfe(j,jj)==properties(1,ii)
                            B=a_fracc(ii);
                            psfe=psfe+min([A B]);
                            a_fracc(i)=a_fracc(i)-min([A B]);
                            a_fracc(ii)=a_fracc(ii)-min([A B]);
                            A=a_fracc(i);
                            break
                        end
                    end
                end
                break
            end
        end
    end
    
    psfe=round(psfe*2*100,2); %output
    if psfe>100
        psfe=100;
    end
    
    %general output - parameters matrix showed in UITable when calculation button is pressed%
    parameters=[Tinput,round(tm,0),round(AW,2),round(rho,2),round(Cp,2),round(Thc,2),... %1-6
        round(r,2),round(xp,2),round(xa,2),round(vec,2),round(ea,2),round(d_r,2),... %7-12
        round(d_xp,3),round(D_xp,3),round(d_xa,3),round(D_xa,3),round(D_vec,2),... %13-17
        round(hmix,2),round(smix,2),round(gmix,2),round(him,2),round(se_bcc,5),round(se_fcc,5),... %18-23
        gamma,omega,lambda,phi_bcc,phi_fcc,eta,D_k,psfe,him_]; %24-32

end
%Reference
%[1] Y. Zhang, Y. Zhou, Solid Solution Formation Criteria for High Entropy Alloys, Mater. Sci. Forum. 561–565 (2007) 1337–1339. https://doi.org/10.4028/www.scientific.net/MSF.561-565.1337.
%[2] M.G. Poletti, L. Battezzati, Electronic and thermodynamic criteria for the occurrence of high entropy alloys in metallic systems, Acta Mater. 75 (2014) 297–306. https://doi.org/10.1016/j.actamat.2014.04.033. 
%[3] S. Fang, X. Xiao, L. Xia, W. Li, Y. Dong, Relationship between the widths of supercooled liquid regions and bond parameters of Mg-based bulk metallic glasses, J. Non. Cryst. Solids. 321 (2003) 120–125. https://doi.org/10.1016/S0022-3093(03)00155-8.
%[4] J.W. Yeh, S.K. Chen, S.J. Lin, J.Y. Gan, T.S. Chin, T.T. Shun, C.H. Tsau, S.Y. Chang, Nanostructured high-entropy alloys with multiple principal elements: Novel alloy design concepts and outcomes, Adv. Eng. Mater. 6 (2004) 299-303+274. https://doi.org/10.1002/adem.200300567.
%[5] O.N. Senkov, D.B. Miracle, A new thermodynamic parameter to predict formation of solid solution or intermetallic phases in high entropy alloys, J. Alloys Compd. 658 (2016) 603–607. https://doi.org/10.1016/j.jallcom.2015.10.279.
%[6] X. Yang, Y. Zhang, Prediction of high-entropy stabilized solid-solution in multi-component alloys, Mater. Chem. Phys. 132 (2012) 233–238. https://doi.org/10.1016/j.matchemphys.2011.11.021.
%[7] M.C. Troparevsky, J.R. Morris, P.R.C. Kent, A.R. Lupini, G.M. Stocks, Criteria for predicting the formation of single-phase high-entropy alloys, Phys. Rev. X. 5 (2015) 1–6. https://doi.org/10.1103/PhysRevX.5.011041.
%[8] Z. Wang, Y. Huang, Y. Yang, J. Wang, C.T. Liu, Atomic-size effect and solid solubility of multicomponent alloys, Scr. Mater. 94 (2015) 28–31. https://doi.org/10.1016/j.scriptamat.2014.09.010.
%[9] A.K. Singh, N. Kumar, A. Dwivedi, A. Subramaniam, A geometrical parameter for the formation of disordered solid solutions in multi-component alloys, Intermetallics. 53 (2014) 112–119. https://doi.org/10.1016/j.intermet.2014.04.019.
%[10] O.N. Senkov, D.B. Miracle, A new thermodynamic parameter to predict formation of solid solution or intermetallic phases in high entropy alloys, J. Alloys Compd. 658 (2016) 603–607. https://doi.org/10.1016/j.jallcom.2015.10.279.
%[11] A. Takeuchi, K. Amiya, T. Wada, K. Yubuta, W. Zhang, A. Makino, Entropies in alloy design for high-entropy and bulk glassy alloys, Entropy. 15 (2013) 3810–3821. https://doi.org/10.3390/e15093810.
%[12] Y.F. Ye, Q. Wang, J. Lu, C.T. Liu, Y. Yang, The generalized thermodynamic rule for phase selection in multicomponent alloys, Intermetallics. 59 (2015) 75–80. https://doi.org/10.1016/j.intermet.2014.12.011.
%[13] G.A. Mansoori, N.F. Carnahan, K.E. Starling, T.W. Leland, Equilibrium Thermodynamic Properties of the Mixture of Hard Spheres, J. Chem. Phys. 54 (1971) 1523–1525. https://doi.org/10.1063/1.1675048.
%[14] M.H. Tsai, K.C. Chang, J.H. Li, R.C. Tsai, A.H. Cheng, A second criterion for sigma phase formation in high-entropy alloys, Mater. Res. Lett. 4 (2015) 90–95. https://doi.org/10.1080/21663831.2015.1121168.
%[15] Y.F. Ye, Q. Wang, J. Lu, C.T. Liu, Y. Yang, Design of high entropy alloys: A single-parameter thermodynamic rule, Scr. Mater. 104 (2015) 53–55. https://doi.org/10.1016/j.scriptamat.2015.03.023.
%[16] M.H. Tsai, K.C. Chang, J.H. Li, R.C. Tsai, A.H. Cheng, A second criterion for sigma phase formation in high-entropy alloys, Mater. Res. Lett. 4 (2015) 90–95. https://doi.org/10.1080/21663831.2015.1121168.
%[17] E.O. Hall, H.S. Algie, The Sigma Phase, Metall. Rev. 11 (1966) 61–88.
