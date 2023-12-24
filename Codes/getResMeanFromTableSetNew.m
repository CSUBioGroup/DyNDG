function [meanTable, stdTable, n_Tables] = getResMeanFromTableSetNew(  TableSet, TableSet_std, n_CV  )
% Input
% TableSet is a cell matrix, each cell is a table
% SpecifiedVarSet  specify the variables being outputed. 
% isdebug    disp warning and other info
% Ouput
% meanTable
% stdTable
    % 对多个table 的数值结果计算平均值 方差
    if ~exist('TableSet_std','var')
        TableSet_std = []  ;
    end    
    if ~exist('n_CV','var')
        n_CV = []  ;
    end    
    
    if ~exist('isdebug','var')
        isdebug = 1;
    end  
    
    if ~exist('SpecifiedVarSet','var')
        SpecifiedVarSet = []; 
    end       
    
    % remove empty cells
    idxSet_non_empty = ~cellfun(@isempty,TableSet) ;   
    TableSet = TableSet(idxSet_non_empty) ; 
    n_Tables = length(TableSet);  %table数量
    if n_Tables==1 && isdebug
        warning('There is only one table. The statistics may be unneccessary.'); 
    elseif n_Tables==0 
        meanTable =table;
        stdTable  =table;
        warning('TableSet is empty. There is no effetive table variable.'); 
        return; 
    end
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
    % 确定将输出的 变量（集合）
    % (1) for specified colmuns/variablenames
    % (2) statistics for scalar quantities: fast
    meanTable      = TableSet{1} ; 
    stdTable       = TableSet{1} ; 
    VarAllTableSet = meanTable.Properties.VariableNames; 
    RowAllTableSet = meanTable.Properties.RowNames  ;
    n_measures = size(meanTable,1); 
    for ii=1:n_measures
        subTable = TableSet{1}{ii,1}{1} ; 
        [n_row,n_col] = size( subTable );
        mat = zeros(n_row,n_col,n_Tables );  
        for jj=1:n_Tables
            subTable = TableSet{jj}{ii,1}{1} ; 
            mat(:,:,jj) = subTable{:,:};        
        end
        %
        mat_std = []; 
        if ~isempty( TableSet_std )
            mat_std = zeros(n_row,n_col,n_Tables ); 
            for jj=1:n_Tables
                subTable = TableSet_std{jj}{ii,1}{1} ;           
                mat_std(:,:,jj) = subTable{:,:};  
            end
        end
     
        % %         mat = mean(mat,3); 
%         if ismember(RowAllTableSet{ii},{'AUROC','AUPRC'})
%             meanTable{ii,1}{1}{:,:} =  mean(mat,3) ; 
%             stdTable{ii,1}{1}{:,:}  =  std(mat,0,3) ;              
%         else
        if ~isempty( TableSet_std ) 
            [m_mat,s_mat ,std_type] = getIntegratedMeanAndStd_Matrix_IN( 0, mat, mat_std, n_CV ,  3 ) ;
        else
            m_mat = mean(mat,3); 
            s_mat = std(mat,0,3) ;   
        end
        meanTable{ii,1}{1}{:,2:end} =  m_mat(:,2:end) ; 
        stdTable{ii,1}{1}{:,2:end}  =  s_mat(:,2:end) ; 
%         end
    end
%     meanTable
%     stdTable
    
% %     meanTable.Properties.RowNames = TableSet{1}.Properties.RowNames;
% %     stdTable.Properties.RowNames  = TableSet{1}.Properties.RowNames;   
% % % disp('tetttttttttttttttttttttt');
end

% % % % % % % % % % % % % % % % % % % % %  
function [meam_Integrated,std_Integrated ,std_type] = getIntegratedMeanAndStd_Matrix_IN( std_type, mean_M, std_M, n_sample_in_dataset_rowvec ,  dim_dataset )
% % S = std(A,w) specifies a weighting scheme for any of the previous syntaxes. When w = 0 (default), 
% % S is normalized by N-1. When w = 1, S is normalized by the number of observations, N. w also can be 
% % a weight vector containing nonnegative elements. In this case, the length of w must equal the length of the dimension over which std is operating.
%% 将分组计算的 均值 和 方差 合并成综合的均值和方差，
% 输入： 需记录下每组数据的均值、方差（及其类型：有偏或无偏）、每组数据的样本数目
%  mean_M      mean values  矩阵 最后一个维度是要累积综合的维度        neccesary
%  std_M       std  values   
%  n_sample_rowvec   # samples      neccesary 
%  std_type     type of std    similar to w in std function； 0  unbias 
% % % Output
% meam_Integrated    col vec
% std_Integrated     col vec 
% Ju Xiang
% 2019-11-27 
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% %     mean_M     = mean_M(:);
% %     n_sample_rowvec = n_sample_rowvec(:); 
%     [n_row, n_dataset ] = size( mean_M ); 
    n_dataset= size( mean_M,  dim_dataset ) ;  
    % %          
    dims = [ones(1,dim_dataset-1 ), n_dataset  ] ; 
    if length(n_sample_in_dataset_rowvec)==1  
        n_sample_in_dataset_rowvec = n_sample_in_dataset_rowvec*ones(1,n_dataset )  ; 
        n_sample_in_dataset_rowvec = reshape( n_sample_in_dataset_rowvec, dims) ; 
    else
        n_sample_in_dataset_rowvec = reshape( n_sample_in_dataset_rowvec, dims) ; 
    end
    % % 
    N = sum(n_sample_in_dataset_rowvec(:)) ;   
    meam_Integrated = sum(mean_M.*n_sample_in_dataset_rowvec, dim_dataset)/( N ); 
    % 
    if exist('std_M','var') && ~isempty(std_M)
        if std_type==0  % unbias 
% %             std_Integrated = (   (std_M.^2)*(n_sample_in_dataset_rowvec-1)' + (( meam_Integrated-mean_M ).^2)*n_sample_in_dataset_rowvec' ) /( N -1 +eps); 
            std_Integrated = (   sum((std_M.^2).*(n_sample_in_dataset_rowvec-1),dim_dataset) + sum((( meam_Integrated-mean_M ).^2).*n_sample_in_dataset_rowvec,dim_dataset) ) /( N -1 +eps ); 

        elseif std_type==1
% %             std_Integrated = (   (std_M.^2).*(n_sample_in_dataset_rowvec) +  (( meam_Integrated-mean_M ).^2).*n_sample_in_dataset_rowvec ) /( N ); 
            std_Integrated = (   sum((std_M.^2).*(n_sample_in_dataset_rowvec),dim_dataset) + sum((( meam_Integrated-mean_M ).^2).*n_sample_in_dataset_rowvec,dim_dataset) ) /( N ); 

        else
            error('no definition')
        end
        std_Integrated = sqrt( std_Integrated );
    else
        std_Integrated = []; 
        std_type       = []; 
        
    end

end



