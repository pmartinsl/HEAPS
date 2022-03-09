%This function estimates the normalized atomic and weight fractions of the
%elements in the alloy%
%This expressed as:
%[a_fraction,w_fraction]=fraction_calculate(n_element,amount,amount_type,atomic_weight)
%
%Input: n_element -> Number of candidate constituent elements of the computed alloys
%       amount -> Stores the amount of each element
%       amount_type -> Selected composition scale
%       atomic_weight -> Atomic weight of elements in the alloy
%Output:a_fraction -> Normalized atomic fraction of elements in the alloy
%       w_fraction -> Normalized weight fraction of elements in the alloy
%
function [a_fraction,w_fraction]=fraction_calculate(n_element,amount,amount_type,atomic_weight)

    if sum(amount)==0
        a_fraction=zeros(1,n_element);
        w_fraction=zeros(1,n_element);
        
    elseif amount_type=="weight percent"
        if sum(amount)~=100
            amount(1)=100-sum(amount(2:n_element));
        end
        
        w_fraction=amount/sum(amount);
        for i=1:1:n_element
            mol_ratio(i)=w_fraction(i)/atomic_weight(i);
        end
        a_fraction=mol_ratio/sum(mol_ratio);
        
    elseif amount_type=="atomic percent" 
        if sum(amount)~=100
            amount(1)=100-sum(amount(2:n_element));
        end
        a_fraction=amount/sum(amount);
        for i=1:1:n_element
            weight(i)=atomic_weight(i)*a_fraction(i);
        end
        w_fraction=weight/sum(weight);
        
    else 
        a_fraction=amount/sum(amount);
        for i=1:1:n_element
            weight(i)=atomic_weight(i)*a_fraction(i);
        end
        w_fraction=weight/sum(weight);
    end

end