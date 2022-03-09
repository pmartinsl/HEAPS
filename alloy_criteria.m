%This function provides the semi-empiric parameters used to the criteria found in literature, assossiated to a
%specific high entropy alloy (HEA) given by the user. The criteria are: MC1.ΔHᵐ-δr, 
%MC2.Ω-δr, MC3.δr-ΔHᵐ, MC4.γ, MC5.δχᴬ-δr, MC6.Λ, MC7.T/Tm-δr-ΔHᵐ, MC8.φ, MC9.η-ΔHᶦᵐ, MC10.Δk,
%LSS1.VEC, LSS2. VEC-e/a, LSS3.T/Tm-VEC, LSS4.φ-VEC, LSS5.T-VEC-r-Δχᴬ-ΔVEC, ImF1.Δχᴾ,
%ImF2.PSFE-VEC, ImF3.δχᴬ-δr, FMP1.VEC
%It is expressed as:
%[criteria]=alloy_criteria(parameters)
%
%Input: parameters -> Calculated parameter values
%Output:criteria -> Evaluation of implemented criteria
%
%The references of the data and equations are indicated at the end of the code
%
function [criteria]=alloy_criteria(parameters)

    Tinput=parameters(1);
    tm=parameters(2);
    r=parameters(7);
    vec=parameters(10);
    ea=parameters(11);
    d_r=parameters(12);
    d_xp=parameters(13);
    D_xp=parameters(14);
    d_xa=parameters(15);
    D_xa=parameters(16);
    D_vec=parameters(17);
    hmix=parameters(18);
    smix=parameters(19);
    him=parameters(21);
    gamma=parameters(24);
    omega=parameters(25);
    lambda=parameters(26);
    phi_bcc=parameters(27);
    phi_fcc=parameters(28);
    eta=parameters(29);
    D_k=parameters(30);
    psfe=parameters(31);
    him_=parameters(32);
    
    %MC1 [1][2]
    if 0.5<d_r && d_r<6.5 && -17.5<hmix && hmix<5
        criteria(2)="SS";
    else
        criteria(2)="IM/BMG";
    end
    
    %MC2 [2]
    if 1.1<=omega && d_r<=6.6
        criteria(3)="SS";
    else
        criteria(3)="IM";
    end
    
    %MC3 [3]
    if -11.6<hmix && hmix<3.2 && d_r<=6.6
        criteria(4)="SS";
    elseif hmix<-12 && 6.4<d_r
        criteria(4)="BMG";
    else
        criteria(4)="Uncertain";
    end
    
    %MC4 [4]
    if gamma<=1.175
        criteria(5)="SS";
    else
        criteria(5)="IM/BMG";
    end
    
    %MC5 [5]
    if 1<d_r && d_r<6 && 3<d_xa && d_xa<6
        criteria(6)="SS";
    else
        criteria(6)="IM";
    end
    
    %MC6 [6]
    if 0.96<=lambda
        criteria(7)="SPSS";
    elseif 0.24<lambda && lambda<0.96
        criteria(7)="MPSS";
    else
        criteria(7)="IM";
    end
    
    %MC7 [7]
    value=Tinput/tm;
    if 0.9<=value
        if -15<=hmix && hmix<=5 && d_r<=6.6
            criteria(8)="SS(0.9<T/Tm)";
        elseif hmix<-15 || 6.6<d_r
            criteria(8)="IM(0.9<T/Tm)";
        else
            criteria(8)="Uncertain";
        end
    elseif 0.5<=value && value<0.9
        if -7.5<=hmix && d_r<=3.3
            criteria(8)="SS(0.5<T/Tm<0.9)";
        elseif -7.5<=hmix || 3.3<d_r
            criteria(8)="IM(0.5<T/Tm<0.9)";
        else
            criteria(8)="Uncertain";
        end
    else
        criteria(8)="Uncertain";
    end
    
    %MC8 [8]
    if phi_bcc>=20 || phi_fcc>=20
        criteria(9)="SPSS";
    else
        criteria(9)="MPSS/IM";
    end
    
    %MC9 [9]
    if 1<eta && him_<3.7
        criteria(10)="SPSS";
    elseif eta<1
        criteria(10)="IM";
    else
        criteria(10)="Uncertain";
    end
    
    %MC10 [10]
    if 0<=D_k
        criteria(11)="SS";
    else
        criteria(11)="IM";
    end
    
    %LSS1 [11]
    if vec<6.87
        criteria(12)="BCC";
    elseif 8<=vec
        criteria(12)="FCC";
    else
        criteria(12)="BCC+FCC";
    end
    
    %LSS2 [5]
    if 7.5<vec && 1.6<ea && ea<1.8
        criteria(13)="FCC";
    elseif vec<7.5 && 1.8<ea && ea<2.3
        criteria(13)="BCC";
    else
        criteria(13)="Uncertain";
    end
    
    %LSS3 [7]
    value=Tinput/tm;
    if 0.9<=value
        if 7.84<vec
            criteria(14)="FCC(0.9<T/Tm)";
        elseif vec<6.87
            criteria(14)="BCC(0.9<T/Tm)";
        else
            criteria(14)="BCC+FCC(0.9<T/Tm)";
        end
    elseif 0.5<=value && value<0.9
        if 7.8<vec
            criteria(14)="FCC(0.5<T/Tm<0.9)";
        elseif vec<6
            criteria(14)="BCC(0.5<T/Tm<0.9)";
        else
            criteria(14)="BCC+FCC(0.5<T/Tm<0.9)";
        end
            
    else
        criteria(14)="Uncertain";
    end
    
    %LSS4 [12]
    if phi_bcc>=20 || phi_fcc>=20
        if 7.5<vec && vec<9.5
            criteria(15)="FCC";
        elseif 4.3<vec && vec<5.7
            criteria(15)="BCC";
        elseif 2.6<vec && vec<3
            criteria(15)="HCP";
        else
            criteria(15)="Uncertain";
        end
    else
        criteria(15)="Uncertain";
    end
    
    %LSS5 [13]
    if 1080<Tinput && Tinput<1660 && 7.8<vec && vec<8.96 && 1.26<r/100 && r/100<1.28 && 0.066<D_xa && D_xa<0.105 && 1.5<D_vec && D_vec<2.26
        criteria(16)="FCC";
    elseif 1330<Tinput && Tinput<1690 && 6.45<vec && vec<7.55 && 1.28<r/100 && r/100<1.3 && 0.071<D_xa && D_xa<0.096 && 1.46<D_vec && D_vec<2.1
        criteria(16)="BCC";
    else
        criteria(16)="Uncertain";
    end
    
    %ImF1 [14]
    if 0.133<D_xp 
        criteria(17)="TCP Phase";
    elseif D_xp<0.117
        criteria(17)="TCP Free";
    else
        criteria(17)="Uncertain";
    end
    
    %ImF2 [15][16]
    if 6.88<=vec && vec<=7.84 && 42.5<psfe
        criteria(18)="Sigma Phase";
    elseif vec<6.88 || 7.84<vec || psfe<22.5
        criteria(18)="Sigma Free";
    else
        criteria(18)="Uncertain";
    end
    
    %ImF3 [17]
    if 7<d_xa && 5<d_r
        criteria(19)="Laves Phase";
    else
        criteria(19)="Laves Free";
    end
    
    %FMP1 [18]
    if vec<4.4
        criteria(20)="Ductile";
    elseif 4.6<vec
        criteria(20)="Brittle";
    else
        criteria(20)="Uncertain";
    end
    
    criteria(1)=num2str(Tinput);
end

%Reference
%[1]    Y. Zhang, Y.J. Zhou, J.P. Lin, G.L. Chen, P.K. Liaw, Adv. Eng. Mater. 10 (2008) 534–538.
%[2]	X. Yang, Y. Zhang, Mater. Chem. Phys. 132 (2012) 233–238.
%[3]	S. Guo, Q. Hu, C. Ng, C.T. Liu, Intermetallics 41 (2013) 96–103.
%[4]	Z. Wang, Y. Huang, Y. Yang, J. Wang, C.T. Liu, Scr. Mater. 94 (2014) 28–31.
%[5]	M.G. Poletti, L. Battezzati, Acta Mater. 75 (2014) 297–306.
%[6]	A.K. Singh, N. Kumar, A. Dwivedi, A. Subramaniam, Intermetallics 53 (2014) 112–119.
%[7]	Z. Wang, S. Guo, C.T. Liu, Jom 66 (2014) 1966–1972.
%[8]	Y.F. Ye, Q. Wang, J. Lu, C.T. Liu, Y. Yang, Scr. Mater. 104 (2015) 53–55.
%[9]	M.C. Troparevsky, J.R. Morris, P.R.C. Kent, A.R. Lupini, G.M. Stocks, Phys. Rev. X 5 (2015) 1–6.
%[10]	O.N. Senkov, D.B. Miracle, J. Alloys Compd. 658 (2016) 603–607.
%[11]   S. Guo, C. Ng, J. Lu, C.T. Liu, J. Appl. Phys. 109 (2011).
%[12]	Y.F. Ye, Q. Wang, J. Lu, C.T. Liu, Y. Yang, Mater. Today 19 (2016) 349–362.
%[13]	Y. Zeng, M. Man, K. Bai, Y.-W. Zhang, Mater. Des. 202 (2021) 109532.
%[14]   Y. Dong, Y. Lu, L. Jiang, T. Wang, T. Li, Intermetallics 52 (2014) 105–109.
%[15]   M.H. Tsai, K.C. Chang, J.H. Li, R.C. Tsai, A.H. Cheng, Mater. Res. Lett. 4 (2015) 90–95
%[16]   M.H. Tsai, K.Y. Tsai, C.W. Tsai, C. Lee, C.C. Juan, J.W. Yeh, Mater. Res. Lett. 1 (2013) 207–212.
%[17]	N. Yurchenko, N. Stepanov, G. Salishchev, Mater. Sci. Technol. (United Kingdom) 33 (2017) 17–22.
%[18]   S. Sheikh, S. Shafeie, Q. Hu, J. Ahlström, C. Persson, J. Veselý, J. Zýka, U. Klement, S. Guo, J. Appl. Phys. 120 (2016).
