function [ TableScores ] = A_DGP_DyNDG(Mdynetset,Mstatic,PStage0, P0, restart, pro_jump_DyNet, pro_jump_Static, tau_Static, InterlayerType )
    % % % % % % % % % % % % % % % % % % 
    if ~exist('restart','var') || isempty (restart)
        restart = 0.7 ; 
    end       
    NormalizationType = 'ProbabilityNormalizationColumn';
if ~isempty(Mdynetset)      
    [A_nLxnL,N_node, L_net] = getMultiplexMatrixFromAdjSet_IN(Mdynetset, pro_jump_DyNet,   'None',InterlayerType  ) ;   
    %
	if  ~isempty( Mstatic )
		delta             = pro_jump_Static;
		DiffusionConstant = 1-pro_jump_Static ;
	     Aeye_col = repmat( (delta/(L_net)).*speye( N_node, N_node), 1,    L_net); 
	     Aeye_row = repmat( (delta/(L_net)).*speye( N_node, N_node), L_net,1    ); 
		if ~isempty( Mstatic )
			A_nLxnL = [ DiffusionConstant.*Mstatic,      Aeye_col   ;...      
						 Aeye_row,                       A_nLxnL   ]; 
			L_net   = L_net + 1 ; 
		end 
    end    
    
else
	A_nLxnL = Mstatic; 
	N_node = length(Mstatic); 
	L_net = 1 ;	
end 
     
    A_nLxnL = sparse( A_nLxnL ) ;
    A_nLxnL = sparse(  getNormalizedMatrix_IN(A_nLxnL, NormalizationType, 1  )  );  
	P0 = reshape(P0,[],1); 
    if ~isempty(Mdynetset)
        stage=size(PStage0);
        P0_DyNet=[];
        for i=1:stage
            P0_DyNet=[P0_DyNet;PStage0{i,1}{1}./sum(PStage0{i,1}{1})];
        end
        if ~isempty( Mstatic )
            P0_L     = [ tau_Static*P0./sum(P0+eps); (1-tau_Static)*P0_DyNet./sum(P0_DyNet+eps) ] ;
        else
            P0_L     = [ (1-tau_Static)*P0_DyNet ] ;
        end
    else
        P0_DyNet = repmat(P0, L_net-1, 1);
        P0_L     = [ tau_Static*P0; (1-tau_Static)*P0_DyNet ] ;
    end
    
      
    N_max_iter     =100; 
    Eps_min_change =10^-6; 
    Pt = P0_L;  %%% size(A_nLxnL), size(P0_L)
    for T = 1: N_max_iter
        Pt1 = (1-restart)*A_nLxnL*Pt + restart*P0_L;
        if all( sum( abs( Pt1-Pt )) < Eps_min_change )
            break;
        end
        Pt = Pt1;
    end
    Pt = full(Pt);     
    % % % % % % % % % % % % % % % %     % 
    Pset =reshape(Pt, N_node, L_net) ;  
	Ranks  = getRankingOfScoreList_IN(   Pset, 'descend' ) ; 	
	
    TableScores = table;  
if ~isempty(Mdynetset)     
	if ~isempty(Mstatic); 
		TableScores.(['DyNDGgm_Dy',InterlayerType])     = 1./geomean( Ranks(:,[2:end]) , 2  );  
		TableScores.(['DyNDGgm',InterlayerType])     = 1./geomean( Ranks , 2  );   
	else
		TableScores.(['DyNDGgm_Dy',InterlayerType])     = 1./geomean( Ranks , 2  );  
	end  
end

if  ~isempty( Mstatic )
    TableScores.(['DyNDGstatic',InterlayerType])  = Pset(:,1);  
end        
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function [Ranks, ord ]= getRankingOfScoreList_IN(   ScoreList, sorttype, OP_IdenticalValues )
    if ~exist('OP_IdenticalValues','var') || isempty(OP_IdenticalValues)
       OP_IdenticalValues = 'MeanRank';  % % for elements with the same values 
    end
    % 
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
                % % 
            end
        end
    end 
    
    if IsSparse_ScoreList   % 
        Ranks = sparse( Ranks );  
    end

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function [A_nLxnL,N_node, L_net] = getMultiplexMatrixFromAdjSet_IN(AdjSet, pro_jump, NormalizationType, InterlayerType, DiffusionConstant)
    if ~exist('NormalizationType','var') || isempty(NormalizationType)
        NormalizationType = 'None'; 
    end 
    SetIsolatedNodeSelfLoop = true;    
    if ~exist('DiffusionConstant','var') || isempty(DiffusionConstant)
        DiffusionConstant = 1-pro_jump; 
    end

    delta =pro_jump;
    %
    if isa(AdjSet,'struct') || isa(AdjSet,'table')
        if isa(AdjSet,'struct')
            fieldnameset = fieldnames( AdjSet); 
        elseif isa(AdjSet,'table')
            fieldnameset = AdjSet.Properties.VariableNames; 
        end
        L_net  =  length( fieldnameset ); 
        N_node =  length( AdjSet.(fieldnameset{1}));
        NxL = N_node*L_net ; 
        %
        if L_net==1;  delta = 0; end %% no jumping for only single layer. 
        %      
        switch  InterlayerType
            case 'Multiplex'
                A_nLxnL = repmat( (delta/(L_net-1)).*speye( N_node, N_node), L_net,L_net);            
                for ii_net = 1: L_net
                    idx = N_node*(ii_net-1)+[1: N_node ] ; 
                    if strcmpi(NormalizationType,'None') 
                        A_nLxnL(idx,idx)=  (DiffusionConstant).*AdjSet.(fieldnameset{ii_net}); 
                    else 
                        A_nLxnL(idx,idx)=  (DiffusionConstant).*getNormalizedMatrix_IN( AdjSet.(fieldnameset{ii_net}) , NormalizationType, SetIsolatedNodeSelfLoop );  
                    end                   
                end  
             
            case 'LINE'     
                A_nLxnL = repmat(  sparse( N_node, N_node), L_net,L_net);       
                ii = [N_node+1:NxL       ]' ;
                jj = [       1:NxL-N_node   ]' ;
                A_nLxnL(sub2ind(size(A_nLxnL),ii,jj)) = delta ;
                A_nLxnL(sub2ind(size(A_nLxnL),jj,ii)) = delta ; 
                %
                for ii_net = 1: L_net
                    idx = N_node*(ii_net-1)+[1: N_node ] ; 
                    if strcmpi(NormalizationType,'None') 
                        A_nLxnL(idx,idx)=  (DiffusionConstant).*AdjSet.(fieldnameset{ii_net}); 
                    else 
                        A_nLxnL(idx,idx)=  (DiffusionConstant).*getNormalizedMatrix_IN( AdjSet.(fieldnameset{ii_net}) , NormalizationType, SetIsolatedNodeSelfLoop );  
                    end                   
                end  
                 
            case {'None','Star'} 
                A_nLxnL = repmat(  sparse( N_node, N_node), L_net,L_net);            
                for ii_net = 1: L_net
                    idx = N_node*(ii_net-1)+[1: N_node ] ; 
                    if strcmpi(NormalizationType,'None') 
                        A_nLxnL(idx,idx)=  (DiffusionConstant).*AdjSet.(fieldnameset{ii_net}); 
                    else 
                        A_nLxnL(idx,idx)=  (DiffusionConstant).*getNormalizedMatrix_IN( AdjSet.(fieldnameset{ii_net}) , NormalizationType, SetIsolatedNodeSelfLoop );  
                    end                   
                end  
                
            otherwise
                error(['LinkType is wrong: ', InterlayerType ]);
                
        end 
          
    else
        error(['AdjSet is wrong. It should be a cell matrix or struct.' ]);
        
    end 
        
end
% % % % % % % % % % % % % 
function WAdj = getNormalizedMatrix_IN(Adj, NormalizationType, SetIsolatedNodeSelfLoop  )   
    if ischar(NormalizationType)
        switch  lower( NormalizationType )
            case lower( { 'column','col',  ...
                    'ProbabilityNormalizationColumn','ProbabilityNormalizationCol',...
                    'ProbabilityColumnNormalization','ProbabilityColNormalization',...
                    'NormalizationColumn','NormalizationCol' , ...
                    'ColumnNormalization','ColNormalization'   })
                NormalizationName = 'ProbabilityNormalization' ;  %  'Random Walk'  
                dim =1;
            case lower({ 'row' ,'ProbabilityNormalizationRow' ,'NormalizationRow' ,'ProbabilityRowNormalization' ,'RowNormalization'   })
                NormalizationName = 'ProbabilityNormalization' ;  %  'Random Walk'  
                dim =2; 
            case lower({'none', 'None', 'NONE'})
                WAdj = Adj; 
                return; 
            otherwise
                error(['There is no type of normalization: ',char( string(NormalizationType) )] );
        end
         
    elseif isempty( NormalizationType )
        WAdj = Adj; 
        return;  
        
    else; error('There is no defintion of NormalizationType')
    end 
    switch lower( NormalizationName )
        case lower( 'ProbabilityNormalization' )
            degrees = sum(Adj,dim);
            if any( degrees~=1)
                WAdj = Adj./ ( degrees+eps  );            
            else
                WAdj = Adj; 
            end
            % 
            if SetIsolatedNodeSelfLoop  && size(Adj,1)==size(Adj,2) 
                ii = find( ~degrees ); 
                idx = sub2ind( size(Adj), ii,ii ); 
                WAdj(idx) = 1;  % set to be 1 for isolated nodes, 
            end 
        case lower( {'None','none'} )
            WAdj = Adj;   
        otherwise
            error(['NormalizationName is wrong: ',char(string(NormalizationName) )   ]);
    end
    
 
 
 
end
 
