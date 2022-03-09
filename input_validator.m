%This function validates the numerical input related to the actions
%'Calculate'%. Additionally, causes a pop-up window to appear in the event of errors in 
%the inputs entered by the user.
%
%Input: n_element -> Number of candidate constituent elements of the computed alloys
%       amount -> Stores the amount 'Fixed' of each element in the alloy
%       amount_type -> Selected composition scale
%       fix -> Vector containing the logical value of the ith 'Fixed' checkbox 
%       from -> Vector containing the lower value of the range in which the
%       ith constituent element may vary
%       stepsize -> Corresponds to the composition step in each numerical iteration
%       to -> Vector containing the upper value of the range in which the
%       ith constituent element may vary
%       filters -> Logic value that stores the state of the 'Apply
%       restrictions' checkbox
%       Prst -> Vector that stores the minimum value of the ith restriction in filters panel
%       PRst -> Vector that stores the maximum value of the ith restriction in filters panel
%       at_fracc -> Vector that stores the minium and maximum allowed
%       atomic fraction of the constituent elements
%       n_comp -> Vector that stores the minimum and maximum allowed number
%       of constituents
%Output:limit -> Vector containing the inconsistencies of the user-entered
%       selections
function [limit]=input_validator(n_element,amount,amount_type,fix,from,stepsize,to,filters,Prst,PRst,at_fracc,n_comp)
    limit=[""];
    i=0; %variable that counts the rows of errors related to the limits of the inputs.
    
    if amount_type=="Composition scale"
        i=i+1;
        limit(i,1)="-You must choose a composition scale";
    elseif amount_type=="Molar ratio"
    else
        if amount(1)<0
            i=i+1;
            limit(i,1)="-The sum of the atomic or weight percentages must be less than or equal to 100";
        end
    end
    
    for j=1:1:n_element
        if ~fix(j)
            if to(j)<=from(j)
                i=i+1;
                limit(i,1)="-'To' values cannot be equal or smaller that 'From' values.";
                break
            end
        end
    end
    
    for j=1:1:n_element
        if ~fix(j)
            var=to(j)-from(j);
            if var<stepsize
                i=i+1;
                limit(i,1)="-'Stepsize' value cannot be bigger that 'To'-'From' value.";
                break
            end
        end
    end
    
    if filters

        if n_comp(1)>n_element%
            i=i+1;
            limit(i,1)="-'Number of components' restrictions cannot be smaller that number of components in the alloy.";%Esto esta malo, porque es el número mínimo de elementos que quieres tenga tu aleación
        end
        
        restrictions=["Melting point [Tm]","Density [ρ]","Radius","χ Pauling","χ Allen","VEC","e/a","δr","δχ Pauling","Δχ Pauling","δχ Allen","Δχ Allen","ΔVEC","ΔHᵐ","ΔSᵐ","ΔGᵐ","ΔHᵐ","Gamma [γ]","Omega [Ω]","Lambda [Λ]","Phi [φ]","Eta [η]","Δk","PSFE"];
        n_rst=length(Prst);
        error="";
        for j=1:1:n_rst
            if PRst(j)<=Prst(j)
              error=error+restrictions(j)+", ";
            end
        end
        
        if at_fracc(2)<at_fracc(1)
            i=i+1;
            error=error+"Atomic fraction, ";
        end
        if n_comp(2)<n_comp(1)
            i=i+1;
            error=error+"Number of components, ";
        end
        %The following condition builds the message related to the inputs associated with the restrictions of the action 'Explore'
        if strlength(error)~=0
            i=i+1;
            limit(i,1)="-The superior limit of cannot be equal or smaller that its inferior limit: "+error+".";
        end
    end
    
    i=0;%variable that counts the rows of errors 
    error=[""]; %saves the errors indicated by the function 'input_validator'
    
    if strlength(limit(1))~=0
        n_rst=length(limit(:,1));
        for j=1:1:n_rst
            i=i+1;
            error(i,1)=limit(j,1);
        end
        %popup window creation
        msgbox(error, 'Exception!','error');
    end
    
end
