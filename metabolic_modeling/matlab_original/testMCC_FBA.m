function [Simuluation_results,Predictions,MCC]=testMCC_FBA(model,Data,met_Csource,met_Nsource)

PMs=Data(:,1);
mets=Data(:,6);
Experiments=Data(:,5);
Csource_ex=['EX_',met_Csource,'_e'];
Nsource_ex=['EX_',met_Nsource,'_e'];
growth_results_csource={};
growth_results_nsource={};
for i=1:length(mets)
    PM=PMs{i};
    met=mets{i};
    met_ex=['EX_',met,'_e'];
    met_rxn_model=findRxnIDs(model,met_ex);
    Exp=Experiments{i};
    if contains(PM,'PM1') || contains(PM,'PM2')
        if met_rxn_model == 0
            growth_results_csource{end+1,1}=PM;
            growth_results_csource{end,2}=met;
            growth_results_csource{end,3}= Exp;
            growth_results_csource{end,4}= 'No exchange rxn';
            growth_results_csource{end,5}= '';
            growth_results_csource{end,6}= 0;
            growth_results_csource{end,7}= 0;
        else
            modelPMc=changeRxnBounds(model,{Csource_ex,Nsource_ex,met_ex},[0,-5,-5],'l');
            FBA=optimizeCbModel(modelPMc);
            growth_results_csource{end+1,1}=PM;
            growth_results_csource{end,2}=met;
            growth_results_csource{end,3}= Exp;
            if FBA.f > 0
                growth_results_csource{end,4}= '+';
                if contains(Exp,'-')
                   growth_results_csource{end,5}= 'FP';
                else
                   growth_results_csource{end,5}= 'TP';
                end
            else
                growth_results_csource{end,4}= '-';
                if contains(Exp,'-')
                   growth_results_csource{end,5}= 'TN';
                else
                   growth_results_csource{end,5}= 'FN';
                end
            end
            growth_results_csource{end,6}=FBA.f;
            growth_results_csource{end,7}=FBA.x;
        end
    else
        if met_rxn_model == 0
            growth_results_nsource{end+1,1}=PM;
            growth_results_nsource{end,2}=met;
            growth_results_nsource{end,3}= Exp;
            growth_results_nsource{end,4}= 'No exchange rxn';
            growth_results_nsource{end,5}= '';
            growth_results_nsource{end,6}= 0;
            growth_results_nsource{end,7}= 0;
        else
            modelPMn=changeRxnBounds(model,{Csource_ex,Nsource_ex,met_ex},[-5,0,-5],'l');
            FBA=optimizeCbModel(modelPMn);
            growth_results_nsource{end+1,1}=PM;
            growth_results_nsource{end,2}=met;
            growth_results_nsource{end,3}= Exp;
            if FBA.f > 0
                growth_results_nsource{end,4}= '+';
                if contains(Exp,'-')
                   growth_results_nsource{end,5}= 'FP';
                else
                   growth_results_nsource{end,5}= 'TP';
                end
            else
                growth_results_nsource{end,4}= '-';
                if contains(Exp,'-')
                   growth_results_nsource{end,5}= 'TN';
                else
                   growth_results_nsource{end,5}= 'FN';
                end
            end
            growth_results_nsource{end,6}=FBA.f;
            growth_results_nsource{end,7}=FBA.x;
        end
    end
end

results=vertcat(growth_results_csource,growth_results_nsource);
PMs=results(:,1);
metID=results(:,2);
Experimental_results=results(:,3);
Simulation_results=results(:,4);
Prediction=results(:,5);
FBA_f=results(:,6);
FBA_x=results(:,7);
Simuluation_results=table(PMs,metID,Experimental_results,Simulation_results,Prediction,FBA_f,FBA_x);

% Getting MCC
%General
TPg=length(strmatch('TP',results(:,5)));
TNg=length(strmatch('TN',results(:,5)));
FNg=length(strmatch('FN',results(:,5)));
FPg=length(strmatch('FP',results(:,5)));
Accuracy=(TPg+TNg)./(TPg+TNg+FNg+FPg);
Predictions=table(TPg,TNg,FNg,FPg,Accuracy,'VariableNames',["TP","TN","FN","FP","Accuracy"],'RowNames',"General");


if ((TPg.*TNg)-(FPg.*FNg))==0
    MCC_general=0;
else
    MCC_general=((TPg.*TNg)-(FPg.*FNg))/sqrt((TPg+FPg).*(TPg+FNg).*(TNg+FPg).*(TNg+FNg));
end



%Carbon sources
TPc=length(strmatch('TP',growth_results_csource(:,5)));
TNc=length(strmatch('TN',growth_results_csource(:,5)));
FNc=length(strmatch('FN',growth_results_csource(:,5)));
FPc=length(strmatch('FP',growth_results_csource(:,5)));
Accuracy=(TPc+TNc)./(TPc+TNc+FNc+FPc);
Predictionsc=table(TPc,TNc,FNc,FPc,Accuracy,'VariableNames',["TP","TN","FN","FP","Accuracy"],'RowNames',"CarbonSources");
Predictions = [Predictions;Predictionsc];

if ((TPc.*TNc)-(FPc.*FNc))==0
    MCC_carbon_sources=0;
else
    MCC_carbon_sources=((TPc.*TNc)-(FPc.*FNc))/sqrt((TPc+FPc).*(TPc+FNc).*(TNc+FPc).*(TNc+FNc));
end


%Nitrogen sources
TPn=length(strmatch('TP',growth_results_nsource(:,5)));
TNn=length(strmatch('TN',growth_results_nsource(:,5)));
FNn=length(strmatch('FN',growth_results_nsource(:,5)));
FPn=length(strmatch('FP',growth_results_nsource(:,5)));
Accuracy=(TPn+TNn)./(TPn+TNn+FNn+FPn);
Predictionsn=table(TPn,TNn,FNn,FPn,Accuracy,'VariableNames',["TP","TN","FN","FP","Accuracy"],'RowNames',"NitrogenSources");
Predictions = [Predictions;Predictionsn];

if ((TPn.*TNn)-(FPn.*FNn))==0
    MCC_nitrogen_sources=0;
else
    MCC_nitrogen_sources=((TPn.*TNn)-(FPn.*FNn))/sqrt((TPn+FPn).*(TPn+FNn).*(TNn+FPn).*(TNn+FNn));
end

MCC=table(MCC_general,MCC_carbon_sources,MCC_nitrogen_sources);



end