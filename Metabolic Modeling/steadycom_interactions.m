%% Code in MATLAB

% Code prepared for understanding the interactions of E. bugandensis
% with its co-existing communities in International Space Station (ISS)
% Author: Pratyay Sengupta
% Clear all variables in the workspace before running this!!

% Input variables
% Line: Path to the models
% Line: Path to the medium file with bounds
% Line: Path to the output file

%% Load the models and the media file

% Models in sbml format and media file is an excel sheet with reaction ids 
% in first column and corresponding lower bounds in second column

% Given path of the models and media

% Example: 
% models_file_path = "S:/IITM-JPL_Enterobacter/metabolic modelling/Model files/F3-2P";
models_file_path = "/path/to/directory/with/models"

% Example: 
% media_path = "S:/IITM-JPL_Enterobacter/metabolic modelling/Model files/SteadyCom_media/F3-2P_SC.xlsx";
media_path = "/path/to/media.xlsx"

old_folder=pwd;
cd(models_file_path)

% Listing all the SBML models in the path
files = dir(fullfile(strcat(models_file_path,'/*.sbml')));

% Loading media file
media_table = readtable(media_path);
media = table2cell(media_table);
for i = 1:length(files)
  models{i} = readCbModel(files(i).name);
end
cd(old_folder)

%% Constraining models according to the provided media

% Constrainting models
con_models = models;
for k = 1:length(models)
    % Finding unconstraint growth
    maxbio_unconstrained(k) = optimizeCbModel(models{k});
    [EX1,up1] = findExcRxns(models{k});
    EXRxns1 = models{k}.rxns(up1);
    
    % Constraining all the exchange reactions with the given lower bounds
    for m = 1:length(EXRxns1)
        if ismember(EXRxns1(m),media(:,1))
            con_models{k} = changeRxnBounds(con_models{k},EXRxns1(m),-10,'l');
        else
            con_models{k} = changeRxnBounds(con_models{k},EXRxns1(m),0,'l');
        end
    end
    % Finding constraint growth
    maxBio(k) = optimizeCbModel(con_models{k});
end 

%% Steadycom analysis

% Pairwise steadycom analysis
for k = 1:length(con_models)
    if k~=length(con_models)
        for j = k+1:length(con_models)
            test1 = con_models{k}.rxns(find(contains(con_models{k}.rxns,'bio')));
            test2 = con_models{j}.rxns(find(contains(con_models{j}.rxns,'bio')));
            R_matrix{k,j} = Steadycom_gen(con_models{k},con_models{j});
        end
    end
end

%% Generation of interaction table

% Preparation of the output table
ind_biomass = nchoosek(length(models),2);
T{ind_biomass+1,13} = {};
T(1,:)= {'Organism A','Organism B','Vbio of A','Vbio of B','Vbio of A in AB','Vbio of B in AB','ratio of A in AB','ratio of B in AB','effect of AB on A','effect of AB on B','Significant effect on A','Significant effect on B','Interaction type'};
l=2;
for i = 1:length(models)
    if i ~= length(models)
        for j = i+1:length(models)
            T{l,1} = models{i}.modelID;
            T{l,2} = models{j}.modelID;
            T{l,3} = maxBio(i).f;
            T{l,4} = maxBio(j).f;
            T{l,5} = R_matrix{i,j}.vBM(1);
            T{l,6} = R_matrix{i,j}.vBM(2);
            T{l,7} = R_matrix{i,j}.BM(1);
            T{l,8} = R_matrix{i,j}.BM(2);
            T{l,9} = round((T{l,5}-T{l,3})/T{l,3},3);
            T{l,10} = round((T{l,6}-T{l,4})/T{l,4},3);
            if T{l,9}<0.1 & T{l,9}>-0.1
                T{l,11}=0;
            else
                T{l,11}=T{l,9};
            end
            if T{l,10}<0.1 & T{l,10}>-0.1
                T{l,12}=0;
            else
                T{l,12}=T{l,10};
            end
            if (T{l,11}<0 & T{l,12}<0)
                T{l,13}='Competition';
            elseif (T{l,11}<0 & T{l,12}>0) | (T{l,11}>0 & T{l,12}<0)
                T{l,13}='Parasites';
            elseif (T{l,11}==0 & T{l,12}>0) | (T{l,11}>0 & T{l,12}==0)
                T{l,13}='Commensals';
            elseif (T{l,11}<0 & T{l,12}==0) | (T{l,11}==0 & T{l,12}<0)
                T{l,13}='Amensals';
            elseif (T{l,11}==0 & T{l,12}==0) 
                T{l,13}='Neutral';
            elseif (T{l,11}>0 & T{l,12}>0)
                T{l,13}='Mutualism';
            end
            l=l+1;
        end
    end
end

% Path to the output
% Example
% xlswrite("S:/IITM-JPL_Enterobacter/metabolic modelling/Model files/SteadyCom_media/F3-2P_Results.xlsx",T);
xlswrite("path/to/output_file.xlsx",T);

%% SteadyCom

% Steadycom between two members 
function R = Steadycom_gen(A,B)
Joint = createMultipleSpeciesModel({A;B},{'A';'B'});
[Joint.infoCom,Joint.indCom] = getMultiSpeciesModelId(Joint,{'A';'B'});
% To find biomass rxns in both models, change the term 'biomass' according
% to models (it varies as Biomass/bio/Bio)

biomassrxns=strcat({'A'; 'B'},[A.rxns(find(contains(A.rxns,'bio')));B.rxns(find(contains(B.rxns,'bio')))]); 
biomassid=findRxnIDs(Joint,biomassrxns);
Joint.infoCom.spBm=biomassrxns;
Joint.indCom.spBm=biomassid;

% Find the exchange rxns of a community in a Joint model and find their 
% lb in actual constrained models
Exrxns_comm=strrep(Joint.infoCom.EXcom,'[u]','_e0'); % Change the nomenclature accordingly
EXRxns1=A.rxns(findExcRxns(A));
EXRxns2=B.rxns(findExcRxns(B));
for i=1:length(Exrxns_comm)
    if ismember(Exrxns_comm(i),EXRxns1)
        excom_lb(i)=A.lb(find(contains(A.rxns,Exrxns_comm(i))));
    elseif ismember(Exrxns_comm(i),EXRxns2)
        excom_lb(i)=B.lb(find(contains(B.rxns,Exrxns_comm(i))));
    end
end

% Change the exchnage rxns of community in joint model as same
% as actual models
Joint.lb(Joint.indCom.EXcom)=excom_lb; % Community

% To find the commuunity biomass using steadycom
[S,R]=SteadyCom(Joint); 
%R.BM % ratio of biomass of two species in a community
%R.vBM % biomass of two species in a community
end
