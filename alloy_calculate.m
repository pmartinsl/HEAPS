%This function allows to extract the values that are shown in the results
%display window, when the calculate button is pressed.%
%
%Input: element -> Vector that stores the ith constituent elements of the
%       computed alloys
%       n_element -> Number of candidate constituent elements of the computed alloys
%       amount_explorer -> Matrix that stores the composition (quantity of
%       the ith constituent element) of each jth computed alloy
%       amount_type -> Selected composition scale
%       filters -> Logic value that stores the state of the 'Apply
%       restrictions' checkbox
%       Prst -> Vector that stores the minimum value of the ith restriction in filters panel
%       PRst -> Vector that stores the maximum value of the ith restriction in filters panel
%       Crst -> Vector that stores the prediction goal of the ith criterion
%       at_fracc -> Vector that stores the minium and maximum allowed
%       atomic fraction of the constituent elements
%       n_comp -> Vector that stores the minimum and maximum allowed number
%       of constituents
%       Tinput -> Entered temperature
%       Ppstate -> Vector that stores the status of the ith parameter
%       restriction
%
%Output:heas -> Iterated alloys
%       m_parameters -> Matrix containing the ith properties/parameters of the jth computed alloys
%       m_criteria -> Matrix containing the ith criteria evaluation of the jth computed alloys
%       m_Afraction -> Matrix containing ith elements atomic fraction of the jth computed alloys
%       m_Wfraction -> Matrix containing ith elements weight fraction of the jth computed alloys
%       properties -> Matrix containing ith elemental properties of each
%       jth constituent element
%
function [heas,m_parameters,m_criteria,m_Afraction,m_Wfraction,properties]=alloy_calculate(element,n_element,amount_explorer,amount_type,filters,Prst,PRst,Crst,at_fracc,n_comp,Tinput,Ppstate)

    iterations=length(amount_explorer(:,1));
    if 1<iterations
        h_text=iterations+" alloys will be evaluated. Initializing ...";
        h=waitbar(0,h_text,'Name','Exploring Alloys','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
        setappdata(h,'canceling',0);
        p=uicontrol('Style', 'pushbutton', 'String', 'Continue','Position', [210 10 60 23],'Callback', 'setappdata(gcbf,''continue'',1)');
        setappdata(h,'continue',0);
        step = 4;
        for k = 1:step  
            if getappdata(h,'canceling')
                iterations=0;
                break
            end
            if getappdata(h,'continue')
                break
            end
            pause(1);
            j=step-k;
            h_text=iterations+" alloys will be evaluated. Initializing "+j+"...";
            waitbar(0,h,h_text);
        end
        delete(p);
    end
    
    %Call data_base function
    [properties,enthalpy,imenthalpy]=data_base(element,n_element);
    
    evaluation=0; %variable that indicates the quantity of alloys that will be shown on the display
    no_evaluation=0; %variable that indicates the quantity of alloys that will not be shown on the display
    
    for i=1:1:iterations
        
        %Creates a pop-up window showing the progress of the iterations
        if 1<iterations
            if i==1
                s=uicontrol('Style', 'pushbutton', 'String', 'Stop','Position', [210 10 60 23],'Callback', 'setappdata(gcbf,''stop'',1)');
            end
            %fraction related to the waitbar that updates iteration progress
            h1=i/iterations;
            %popup window text that updates iteration progress
            h_text=i+" of "+iterations+" iterations";
            %popup window update
            waitbar(h1,h,h_text);
        end

        %Variable that contains the quantities of the alloy to iterate
        i_amount=amount_explorer(i,:);
        
        if i_amount(1)<0
            eye=0;
        else
            eye=1;
            %Call fraction_calculate function
            [a_fraction,w_fraction]=fraction_calculate(n_element,i_amount,amount_type,properties(9,:));
            %Call parameters_calculate function
            [parameters]=parameters_calculate(a_fraction,n_element,properties,enthalpy,imenthalpy,Tinput);
            %Call alloy_criteria function
            [criteria]=alloy_criteria(parameters);
            parameters=parameters(1:31);
        end
        hea=""; %variable that updates the simplified name of the iterated alloy
        rdx=0;
        if eye
            for k=1:1:n_element
                if filters
                    if at_fracc(2)<a_fraction(k) || a_fraction(k)<at_fracc(1)
                        eye=0;
                    end
                end
                hea=hea+element(k)+i_amount(k); 
            end
            %Returns a 1's and 0's array. If two compared elements are different, it returns a 1. Else, if the element are equals, it returns a 0.%
            idx = a_fraction~=0;
            %Column vector that contains the addition of each row of idx.
            rdx = sum(idx,2);
        end
        
        %variable that updates the number of evaluations
        L=evaluation;
        %Conditions to display the values related to the iterated alloy
        if ~filters && rdx>=3 && eye
            evaluation=evaluation+1;
            
        elseif n_comp(1)<=rdx && rdx<=n_comp(2) && eye
            prst=[parameters(2),parameters(4),parameters(7),parameters(8),...
                parameters(9),parameters(10),parameters(11),parameters(12),...
                parameters(13),parameters(14),parameters(15),parameters(16),...
                parameters(17),parameters(18),parameters(19),parameters(20),...
                parameters(21),parameters(24),parameters(25),parameters(26),...
                parameters(27),parameters(28),parameters(29),parameters(30),parameters(31)];
            crst=criteria(2:end);
            
            for j=1:1:length(Crst)
                if Crst{j}=="Select"

                elseif Crst{j}~=crst(j)
                    eye=0;
                    break
                end
            end
            
            for j=1:1:length(Prst)
                if Ppstate(j)=='on'
                    
                    if Prst(j)<=prst(j) && prst(j)<=PRst(j) && eye  %In paper

                    else
                        eye=0;
                        break
                    end
                    
                       
                    
                end
            end
            if eye
                evaluation=evaluation+1;
            end
        end
        %If the alloy is evaluated, it saves the values in the respective output
        if L<evaluation
            heas(evaluation,:)=hea;
            m_parameters(evaluation,:)=parameters;
            m_criteria(evaluation,:)=criteria;
            m_Afraction(evaluation,:)=round(a_fraction,2);
            m_Wfraction(evaluation,:)=round(w_fraction,2);
        else
            no_evaluation=no_evaluation+1;
        end
        if 1<iterations
            if getappdata(h,'canceling')
                iterations=0;
                break
            end
            if getappdata(h,'stop')
                iterations=evaluation+no_evaluation;
                break
            end
        end
    end
    %modifies the popup window when iterations are performed
    if 1<iterations
        delete(h);
        if iterations==no_evaluation
            h_text=iterations+" alloys were evaluated, but none of them satisfy the constraints";
            heas=char.empty;
            m_parameters=char.empty;
            m_criteria=char.empty;
            m_Afraction=char.empty;
            m_Wfraction=char.empty;
        else
            h_text=+iterations+" alloys were evaluated";
            if filters
                h_text=h_text+", "+evaluation+" alloys are shown";
                if no_evaluation>0
                    h_text=h_text+" and "+no_evaluation+" were discarted";
                end
            end
        end
        waitbar(1,h_text,'Name','Resume');
    elseif iterations==0
        delete(h);
        h_text="Cancel Evaluating";
        heas=char.empty;
        m_parameters=char.empty;
        m_criteria=char.empty;
        m_Afraction=char.empty;
        m_Wfraction=char.empty;
        properties=char.empty;
        waitbar(1,h_text,'Name','Resume');
    end
    
end


    