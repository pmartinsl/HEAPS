%This function defines the numerical values according to the calculation
%parameters entered by the user
%This expressed as:
%amount_explorer=amount_variator(n_element,amount,amount_type,fix,from,stepsize,to)
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
%
%Output:amount_explorer -> Matrix that stores the composition (quantity of
%       the ith constituent element) of each jth computed alloy
%
function amount_explorer=amount_variator(n_element,amount,amount_type,fix,from,stepsize,to)

    variator=ones(1,n_element);
    if amount_type=="Molar ratio"
        start=1;
    else
        start=2;
    end

    for i=start:1:n_element
        if fix(i)
        else
            variator(i)=floor((to(i)-from(i))/stepsize+1);
        end
    end
    iterations=prod(variator);

    if iterations==1
        amount_explorer=amount(1:n_element);
        return
    end

    explorer=zeros(iterations,n_element);

    for i=n_element:-1:1
            if variator(i)==1
                explorer(:,i)=amount(i);
            else
                repetitions=prod(variator(i+1:end));
                loop=prod(variator(1:i-1));
                inicio=1;
                for k=1:1:loop
                    in=inicio;
                    for j=1:1:variator(i)
                        en=in+repetitions-1;
                        explorer(in:en,i)=from(i)+stepsize*(j-1);
                        in=en+1;
                    end
                    inicio=k*repetitions*variator(i)+1;
                end
            end   
    end
    amount_explorer=unique(explorer,'rows');
    if amount_type=="Molar ratio"
    else
        k=0;
        amount_exploreri=[];
        for i=1:1:iterations
            ss=sum(amount_explorer(i,2:end));
            if ss<=100
                k=k+1;
                amount_explorer(i,1)=100-ss;
                amount_exploreri(k,:)=amount_explorer(i,:);
            end
        end
        amount_explorer=amount_exploreri;
    end

end