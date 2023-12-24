function [Ranks, ord ]= getRankingOfScoreList(   ScoreList, sorttype, OP_IdenticalValues )
if ~exist('ScoreList','var')
    sorttype= 'descend'; 
    ScoreList = rand( 5,3) ; 
    ScoreList([ 3 4 5],:)=2;
    ScoreList 
    OP_IdenticalValues ='RandPermutation'; 
end
    if ~exist('OP_IdenticalValues','var') || isempty(OP_IdenticalValues)
       OP_IdenticalValues = 'MeanRank';  
    end
    [~, ord ] = sort(ScoreList , 1,  sorttype );  
    
    IsSparse_ScoreList = issparse( ScoreList );
    if IsSparse_ScoreList  
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
                % 
                if ~strcmpi(OP_IdenticalValues, 'None')  
                    [uniqueScores, ~,ic] = unique( score_t ) ;
                    if length( ic  ) ~= length( uniqueScores  )
                        for ii_uniqueScorese = 1: length( uniqueScores  ) 
                            idx = ( ic== ii_uniqueScorese ) ; 
                            n_thisscore = nnz( idx )   ;
                            if n_thisscore>1
                                if strcmpi(OP_IdenticalValues, 'MeanRank') 
                                    rank_t( idx ) = mean(rank_t(idx), 1); 
                                elseif strcmpi(OP_IdenticalValues, 'RandPermutation') 
                                    ind = find(idx); 
                                    ind_randperm = ind(  randperm( n_thisscore  )  ); 
                                    rank_t( ind ) = rank_t( ind_randperm ); 
                                else
                                    error('There is no definition of OP_IdenticalValues');
                                end
                            end  
                        end
                    end
                end
                Ranks( : ,   d2,d3,d4) = rank_t;  
            end
        end
    end 
    
    if IsSparse_ScoreList   
        Ranks = sparse( Ranks );  
    end

end
