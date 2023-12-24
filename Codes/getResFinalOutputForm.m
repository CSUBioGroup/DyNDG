function [ResStat,ResVecAllDis,ResScalarAllDis, n_Tables] = ...
    getResFinalOutputForm(  TableDiseases_cut, ResTablePerfCurve_AllDiseaseMeanCV, ResTablePerfCurve_AllDiseaseMeanCV_std, nCV_list, methodset  )
% Input
% TableSet is a cell matrix, each cell is a table
    %%
    % % average all scalars and curves over all diseases 
    [meanDISmeanCVTable, stdDISmeanCVTable, n_disease_in_Table] = getResMeanFromTableSetNew(  ResTablePerfCurve_AllDiseaseMeanCV, ResTablePerfCurve_AllDiseaseMeanCV_std, nCV_list  ) ; 
    
% % % % % % % % % % % % % % % % % % %     
    ResStat.mean = meanDISmeanCVTable; 
    ResStat.std  = stdDISmeanCVTable; 
% %     rownames = ResStat.mean.Properties.RowNames; 
    is_scalar = logical( ResStat.mean{:,'is_scalar'} ); 
    t_ResStat = ResStat; 
    ResStat.mean(is_scalar,:) = []; 
    ResStat.std(is_scalar,:)  = []; 
    %
    t_ResStat.mean(~is_scalar,:) = []; 
    t_ResStat.std(~is_scalar,:)  = []; 
    
    rownames = strcat(t_ResStat.mean{:,'type'}, t_ResStat.mean{:,'value'} ); 
    m_tb   = t_ResStat.mean{is_scalar,'table'}; 
    std_tb = t_ResStat.std{is_scalar,'table'};
    for ii = 1:length( m_tb )
        m_tb{ii}.Properties.VariableNames = ['ref',methodset ]; 
        std_tb{ii}.Properties.VariableNames = ['ref',methodset ]; 
    end  
    m_tb   = cat(1,m_tb{:});   m_tb.Properties.RowNames= rownames  ; 
    std_tb = cat(1,std_tb{:}); std_tb.Properties.RowNames= rownames  ;  
    
    ResStat.mean = [ResStat.mean; cell2table( {{m_tb}, 'scalar','all', 1} ,'VariableNames',ResStat.mean.Properties.VariableNames ,'RowNames',{'scalar'} ) ];
    ResStat.std = [ResStat.std; cell2table(  {{std_tb}, 'scalar','all', 1},'VariableNames',ResStat.std.Properties.VariableNames ,'RowNames',{'scalar'} )];
    %% 
    % extract curves and scalar for list of diseases 
    ResScalarTablePerfCurve_AllDisease     = ResTablePerfCurve_AllDiseaseMeanCV; 
    ResScalarTablePerfCurve_AllDisease_std = ResTablePerfCurve_AllDiseaseMeanCV_std; 
    %
    ResVecTablePerfCurve_AllDisease     = ResTablePerfCurve_AllDiseaseMeanCV; 
    ResVecTablePerfCurve_AllDisease_std = ResTablePerfCurve_AllDiseaseMeanCV_std; 
    %
    is_scalar = logical( ResScalarTablePerfCurve_AllDisease{1}{:,'is_scalar'} ); 
    for ii_dis = 1:n_disease_in_Table 
        ResScalarTablePerfCurve_AllDisease{ii_dis}(~is_scalar,:)     = []; 
        ResScalarTablePerfCurve_AllDisease_std{ii_dis}(~is_scalar,:) = []; 
        %
        ResVecTablePerfCurve_AllDisease{ii_dis}(is_scalar,:)     = []; 
        ResVecTablePerfCurve_AllDisease_std{ii_dis}(is_scalar,:) = []; 
    end
    diseases = TableDiseases_cut;
    % for save curves of each diseases  
    ResVecTablePerfCurve_AllDisease     = cell2table(ResVecTablePerfCurve_AllDisease,'RowNames',diseases{:,1}','VariableNames',{'table'}); 
    ResVecAllDis.mean= [ResVecTablePerfCurve_AllDisease, diseases]     ;  
    ResVecTablePerfCurve_AllDisease_std = cell2table(ResVecTablePerfCurve_AllDisease_std,'RowNames',diseases{:,1}' ,'VariableNames',{'table'}); 
    ResVecAllDis.std = [ResVecTablePerfCurve_AllDisease_std, diseases]     ;  

    %
    TableSet = ResScalarTablePerfCurve_AllDisease ; 
    TableSet_std = ResScalarTablePerfCurve_AllDisease_std ; 
    meanTable      = TableSet{1} ; 
    stdTable       = TableSet_std{1} ; 
    n_measures = size(meanTable,1); 
    n_Tables = length( TableSet ) ;
    for ii=1:n_measures
%         subTable = TableSet{1}{ii,1}{1} ; 
%         [n_row,n_col] = size( subTable );
        subTableSet =repmat({table},n_Tables,1); 
        for jj=1:n_Tables
            subTableSet{jj} = TableSet{jj}{ii,1}{1} ;           
        end
        subTable = cat(1, subTableSet{:} ); 
        subTable.Properties.RowNames = reshape(diseases{:,1},1,[]);
        meanTable{ii,'table'} = {subTable}; 
        %
        subTableSet =repmat({table},n_Tables,1); 
        for jj=1:n_Tables
            subTableSet{jj} = TableSet_std{jj}{ii,1}{1} ;           
        end
        subTable = cat(1, subTableSet{:} ); 
        subTable.Properties.RowNames = reshape(diseases{:,1},1,[]);
        stdTable{ii,'table'} = {subTable};  
    end
    ResScalarAllDis.mean = meanTable; 
    ResScalarAllDis.std  = stdTable; 
end

 


