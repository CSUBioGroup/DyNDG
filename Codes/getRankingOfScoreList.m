function [Ranks, ord ]= getRankingOfScoreList(   ScoreList, sorttype, OP_IdenticalValues )
% by Xiang 
% 2019-2
% % Ranks = getRankingOfScoreList(   ScoreList,  'descend', true )
if ~exist('ScoreList','var')
    sorttype= 'descend'; 
    ScoreList = rand( 5,3) ; 
    ScoreList([ 3 4 5],:)=2;
    ScoreList 
    OP_IdenticalValues ='RandPermutation'; 
%     OP_IdenticalValues ='MeanRank'; 
%     OP_IdenticalValues ='none'; 
%     ScoreList=sparse(ScoreList);
%     class( ScoreList )IsMeanRank
end
    if ~exist('OP_IdenticalValues','var') || isempty(OP_IdenticalValues)
    % 对于相同值的元素赋予平均的rank 值
       OP_IdenticalValues = 'MeanRank';  
    end
    %
% %     if strcmpi(OP_IdenticalValues,'RandPerturbation')  % 对于相同元素值，如果不采用平均 rank，则添加随机扰动
% %         ScoreList = ScoreList +  ( ( abs(ScoreList) +realmin/10^-15 )*10.^-10 ).*rand( size(ScoreList ) );   
% %     end
    [~, ord ] = sort(ScoreList , 1,  sorttype );  
    
    IsSparse_ScoreList = issparse( ScoreList );
    if IsSparse_ScoreList  
        % 对于稀疏矩阵 会出现 N-dimensional indexing allowed for full matrices only.
        % 如果 对大小等于1 的、其它维度索引 
        ScoreList = full( ScoreList ); 
    end
    
    IDlist = [1: size( ScoreList, 1 ) ]'; 
    Ranks = zeros( size( ScoreList ) ); 
    rank_t = zeros( size(IDlist) ); 
    for d2=1:size( ScoreList ,2 )
        for d3 = 1:size( ScoreList ,3)
            for d4 = 1:size( ScoreList ,4)
                score_t    = ScoreList(:,d2,d3,d4); 
                rank_t(ord(:,d2,d3,d4)) = IDlist; 
                % Ranks( ord(:,d2,d3,d4),   d2,d3,d4) = IDlist ; 
                % 
                if ~strcmpi(OP_IdenticalValues, 'None')   % 对于相同元素值，采取操作 
% %                 if strcmpi(OP_IdenticalValues, 'MeanRank') || strcmpi(OP_IdenticalValues, 'RandPermutation')  % 对于相同元素值，采用平均 rank 
                    [uniqueScores, ~,ic] = unique( score_t ) ;
                    if length( ic  ) ~= length( uniqueScores  )
                        for ii_uniqueScorese = 1: length( uniqueScores  ) 
                            idx = ( ic== ii_uniqueScorese ) ; 
                            n_thisscore = nnz( idx )   ;
                            if n_thisscore>1
                                if strcmpi(OP_IdenticalValues, 'MeanRank') 
                                    % Ranks( idx,   d2,d3,d4 ) = mean(Ranks(idx,   d2,d3,d4), 1);  %对于相同分值的元素赋予相同的所有ranks的均值
                                    rank_t( idx ) = mean(rank_t(idx), 1);  %对于相同分值的元素赋予相同的所有ranks的均值
                                    % % sum( labels_ord( idx ) )/nnz( idx )
                                elseif strcmpi(OP_IdenticalValues, 'RandPermutation') 
                                    ind = find(idx); 
                                    ind_randperm = ind(  randperm( n_thisscore  )  ); 
                                    rank_t( ind ) = rank_t( ind_randperm ); 
                                else
                                    error('There is no definition of OP_IdenticalValues');
                                end
                            end  
                            % %sum(labels_ord) 
                        end
                    end
                end
                Ranks( : ,   d2,d3,d4) = rank_t;  
                % % 
            end
        end
    end 
    
    if IsSparse_ScoreList   %如果输入为稀疏矩阵，输出rank 也转化成 稀疏形式 
        Ranks = sparse( Ranks );  
    end

end
