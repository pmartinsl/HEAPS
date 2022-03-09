%This fuction corrects the syntaxis of the alloy entered as input, allowing
%to calculate the semi-empiric parameters and properties alloy. Also, 
%indentify if any element entered by the user isn't in the data base or if 
%any elements is repeating%
%It is expressed as:
%[hea,f_atom,elements,index,n_elements,noelem,repeat]=alloy_validator(alloy)
%
%Input: import -> 
%       Tfix -> Entered temperature
%       amount_type -> Selected composition scale
%Output:heas -> iterated alloys
%       m_parameters -> Matrix with the propierties and parameters of the iterated alloys
%       m_critieria -> Matrix with the criteria of the iterated alloys
%       m_element -> Number of candidate constituent elements of the computed alloys
%       m_Afractions -> Matrix with atomic fractions of the elements of the iterated alloys
%       m_Wfractions -> Matrix with weight fractions of the elements of the iterated alloys
%       m_properties -> Properties of the elements in the alloy
%
function [m_heas,m_parameters,m_criteria,m_element,m_Afraction,m_Wfraction,m_properties]=import_function(import,Tfix,amount_type)
    
    element_data=["Li" "Be" "B" "C" "N" "Na" "Mg" "Al" "Si" "P" "K" "Ca" "Sc" "Ti" ...
        "V" "Cr" "Mn" "Fe" "Co" "Ni" "Cu" "Zn" "Ga" "Ge" "As" "Rb" "Sr" "Y" "Zr" ...
        "Nb" "Mo" "Tc" "Ru" "Rh" "Pd" "Ag" "Cd" "In" "Sn" "Sb" "Cs" "Ba" "La" "Hf" ... 
        "Ta" "W" "Re" "Os" "Ir" "Pt" "Au" "Hg" "Tl" "Pb" "Bi"];
    
    at=length(amount_type);
    str="";
    for i=1:1:at
        str=str+amount_type(i);
    end
    amount_type=str;
    %import = readcell('Import.xlsx');
    data=length(element_data)+1;
    n_import=height(import(:,1));
    
    P=0;
    m_element=[""];
    m_Afraction=[];
    m_Wfraction=[];
    m_properties=[];
    
    heas=import(:,1);
    Tinput=import(:,2);
    
    for i=1:1:n_import
        if isnan(heas{i})
            continue
        else    
            [element,amount]=regexp(heas{i},'[A-Z][a-z]?','match','split');
        end
        hea=""; %Saves the alloy corrected without spaces and with atomic fraction equal to 1 when this parameter isn't entered by the user.%

%Transforms the number saved as string into a numeric value. The strings that doesn't contain a associated number are saved as an index without value (NaN). Also, it eliminates the first column of the 'regexp' function output.%
        amount=str2double(amount(2:end));

%Every index without value (NaN) is transformed into a 1. It considers the equimolar alloys written as: 'CuNiCo', 'Cu2NiCo' y 'CuNi2Co2'.%
        if amount_type=='Molar ratio'
            amount(isnan(amount))=1;
        else
            amount(isnan(amount))=10;
        end
        n_element=length(element); %Number of the alloy components.%
        
        if ischar(Tinput{i})
            Tinput{i}=Tfix;
        elseif isnan(Tinput{i})
            Tinput{i}=Tfix;
        end
        
        no_calculate=false;
%the following loop validates the alloy entered, against elements that are not in the database.
        for j=1:1:n_element
            for k=1:1:data
                if k==data
                    no_calculate=true; %Saves the elements as string that aren't contained in the database array.
                    hea=hea+element(j)+amount(j);
                elseif element(j)==element_data(k)
                    hea=hea+element(j)+amount(j); 
                    break %If it's found any element that is contained in the database array, it continues to the next index in the 'elements' vector.% 
                end
            end
        end
    %The following condition identifies repeating elements in the alloy.
        if n_element~=length(unique(element(:)))
            no_calculate=true;
        end
        
        m_heas(i,1)=hea;
        if no_calculate
            m_parameters(i,:)=NaN;
            m_criteria(i,:)="";
            m_Afraction(i,:)=NaN;
            m_Wfraction(i,:)=NaN;
        else
            [properties,enthalpy,imenthalpy]=data_base(element,n_element);
            [a_fraction,w_fraction]=fraction_calculate(n_element,amount,amount_type,properties);
            [parameters]=parameters_calculate(a_fraction,n_element,properties,enthalpy,imenthalpy,Tinput{i});
            [criteria]=alloy_criteria(parameters);
            parameters=parameters(1:31);
            m_parameters(i,:)=parameters;
            m_criteria(i,:)=criteria;
            for j=1:1:n_element
                if element(j)~=m_element
                    P=P+1;
                    m_element(P)=element(j);
                end
            end
            
            for j=1:1:n_element
                for k=1:1:P
                    if element(j)==m_element(k)
                        m_Afraction(i,k)=round(a_fraction(j),2);
                        m_Wfraction(i,k)=round(w_fraction(j),2);
                        m_properties(:,k)=round(properties(:,j),2);
                        break
                    end
                end
            end
        end  
    end
end
