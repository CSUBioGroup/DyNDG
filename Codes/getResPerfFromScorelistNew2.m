% % function  [ResScalartable, ResCurvetable ,ResTopKtable, n_gene,n_dis] = getResPerfFromScorelist(M_gene_dis_postestlabel,Score_gene_dis, topk_set, topk_AP)
function  [Restable, n_gene,n_method] = getResPerfFromScorelistNew2(testsetlabel,TableScores, ScalarSet, CurveSet )
    %% calculate AUC for one disease each time ....  % % % % % % % % % % % % % % % % % % % % % % %
if ismember('all', ScalarSet )|| ismember('ALL', ScalarSet )
    ScalarSet = {'AUROC','AUPRC', 'MRank', 'MRR',      'AURecall', 'AUPrec',   'topk_recall','topk_prec','topk_recall_p','topk_prec_p'} ;
end
if ismember('all', CurveSet )|| ismember('ALL', CurveSet )
    CurveSet = {'CurveTopKvsRecall','CurveTopKvsPrec',      'CurveTopKvsRecall_p','CurveTopKvsPrec_p',    'ROC','PRC'}; 
end
  
    %
    methodnames = TableScores.Properties.VariableNames ; 
    [n_gene,n_method ] = size(TableScores) ;      
    XVals = [0:0.001:1]; n_XVals = length(XVals); 
    ROCy  = zeros(n_XVals ,n_method ) ; 
    PRCy  = zeros(n_XVals ,n_method ) ;   
    AUROCset = 0.5.*ones( 1, n_method    );
    AUPRCset = 0.5.*ones( 1, n_method    );
    
    % top k curves for count 
    top_k_vec = (1:200)' ; len_topk  = length( top_k_vec  ); 
    CurveRecall_topk = zeros(len_topk ,n_method ) ; 
    CurvePrec_topk   = zeros(len_topk ,n_method ) ; 
    %top k curves for percent 
    top_k_vec_percent = [0.01:0.01:0.5]' ; len_topk_percent  = length( top_k_vec_percent  ); 
    CurveRecall_topk_percent = zeros(len_topk_percent ,n_method ) ; 
    CurvePrec_topk_percent   = zeros(len_topk_percent ,n_method ) ; 
    
    % top k scalar for count 
    Use_topk_percent = 0;  topk_scalarset = [1  5   10      50   100   200]' ;         
    n_topk_scalar = length(topk_scalarset);   
    ScalarSet_Recall_topk = zeros(n_topk_scalar ,n_method ) ; 
    ScalarSet_Prec_topk   = zeros(n_topk_scalar ,n_method ) ;  
    %top k scalar for percent 
    Use_topk_percent = 1;  topk_scalarset_percent = [0.01  0.05    0.10  0.20     0.50     ]' ;  
    n_topk_scalar_percent = length(topk_scalarset_percent);   
    ScalarSet_Recall_topk_percent = zeros(n_topk_scalar_percent ,n_method ) ; 
    ScalarSet_Prec_topk_percent   = zeros(n_topk_scalar_percent ,n_method ) ;       
    %
    % % Dividing the rank of a test gene by the total number of test and control genes in a validation
    % % run, we obtained the rank ratio of the test gene. 
    % % Averaging rank ratios of all test genes, we obtained a criterion called the mean rank ratio (MRR).
    MRRset  = -1.*ones( 1, n_method    );    % mean rank ratio (MRR)
    MRankset = -1.*ones( 1, n_method    );   % mean rank of test gene(s)  
    labels = testsetlabel;
    for ii=1:n_method  
        scores = double(  TableScores{:, ii }    );  
        N_pos  = nnz(  labels   ); 
        % AUROC 
        [X,Y,~,AUC ] = perfcurve( labels , scores ,1 );  
        if 0<X(1) % add  (0,0) if noexisted
            X = [0 ; X ];
            Y = [0;  Y ];
        end 
        X(end) =1;       
        [X,Y] = getUniformXY(X,Y,XVals, [],   [], 'linearperturbation' , 0 ) ; 
        AUROCset(ii) = AUC ;  
        ROCy(:,ii) = Y;    
        % %         ROCy = ROCy + Y;   
         
        %AUPRC 
        [X,Y, ~,AUPRC  ] = perfcurve( labels , scores , 1,  'XCrit','tpr', 'YCrit','prec' );
        if isnan(Y(1)); Y(1) = Y(2) ; end 
        if 0<X(1) 
            X = [0 ; X ];
            Y = [Y(1); Y ];   % N_pos/n_label
        end 
        X(end) =1 ; 
        [X,Y] = getUniformXY(X,Y,XVals, [],   [], 'linearperturbation' , 0 ) ; 
        AUPRCset(ii) = trapz(X, Y);   
        PRCy(:,ii) = Y;       
  
        % top k    reacall   precision 
        top_k_vec_t = top_k_vec; top_k_vec_t(top_k_vec_t>n_gene)=n_gene;
        [X00,Y00 ] = perfcurve( labels , scores , 1,  'XCrit','tp+fp', 'YCrit','tp' ) ;  
        [~,Y] = getUniformXY(X00,Y00,top_k_vec_t, [],   [], 'linearperturbation' , 0 ) ; 
        CurveRecall_topk(: ,ii) = Y./N_pos;       
        CurvePrec_topk(: ,ii)   = Y./top_k_vec;  
        
        %
        [~,Y] = getUniformXY(X00./length(X00),Y00,top_k_vec_percent, [],   [], 'linearperturbation' , 0 ) ; 
        CurveRecall_topk_percent(: ,ii) = Y./N_pos;       
        CurvePrec_topk_percent(: ,ii)   = Y./(length(X00).*top_k_vec_percent+eps);   
        
        % top k  scalar output   
        topk_scalarset_t = topk_scalarset; topk_scalarset_t(topk_scalarset_t>n_gene)=n_gene;
        [~,Y] = getUniformXY(X00,Y00,topk_scalarset_t, [],   [], 'linearperturbation' , 0 ) ; 
        ScalarSet_Recall_topk(: ,ii) = Y./N_pos;
        ScalarSet_Prec_topk(: ,ii)   = Y./( topk_scalarset);   
        % top k  scalar _percent
        [~,Y] = getUniformXY(X00./length(X00),Y00,topk_scalarset_percent, [],   [], 'linearperturbation' , 0 ) ;             
        ScalarSet_Recall_topk_percent(: ,ii) = Y./N_pos;
        ScalarSet_Prec_topk_percent(: ,ii)   = Y./(length(X00).*topk_scalarset_percent);   
 
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         % % % % % % 
        sorttype= 'descend';
        [~, ord ] = sort(scores , 1,  sorttype ); 
        IDlist = [1: n_gene ]';  
        rank_t = zeros( size(IDlist) );   
        rank_t(ord) = IDlist;     
        %sumrank = sum( rank_t( logical(labels) ) )  -   sum( 0:(N_pos-1) ); 
        %MRankset(ii)  = (sumrank )./N_pos            % mean rank (MRank) of test gene(s)  
        %MRRset(ii)    = (sumrank )./N_pos./n_gene    % mean rank ratio (MRR)  
        MRankset(ii)  = mean( rank_t( logical(labels) ) ) - (N_pos-1)/2 ;   % mean rank (MRank) of test gene(s), corrected for multiple test genes 
        MRRset(ii)    = MRankset(ii)./( n_gene - (N_pos-1) ) ;    % mean rank ratio (MRR), corrected for multiple test genes 		
        % % % % % % % % % % % % % % % % % % % % % 
        
    end 
     
    %%  [ ResAUROC, ResAUPRC ] 
    is_scalar = []; kk = 0; 
    % save to table 
    Restable = table;  
    
    %
    type = 'AUROC'; rowname = type; 
    ResAUROC = array2table([-1, AUROCset], 'VariableNames',['none',strcat('AUROC_',methodnames) ]);  
    Restable{rowname,{'table','type','value'}} = {ResAUROC,type,'all' }; 
    kk = kk + 1; is_scalar(kk,1) = true; 
    
    %
    type = 'AUPRC'; rowname = type; 
    ResAUPRC = array2table([-1,AUPRCset], 'VariableNames',['none',strcat('AUPRC_',methodnames)] );
    Restable{rowname,{'table','type','value'}} = {ResAUPRC,type,'all' }; 
    kk = kk + 1; is_scalar(kk,1) = true; 

    % mean rank (MRank) of test gene(s)  
    type = 'MRank'; rowname = type; 
    ResMRank = array2table([-1, MRankset], 'VariableNames',['none',strcat('MRank_',methodnames) ]);  
    Restable{rowname,{'table','type','value'}} = {ResMRank,type,'all' }; 
    kk = kk + 1; is_scalar(kk,1) = true; 
    
    % mean rank ratio (MRR)  
    type = 'MRR'; rowname = type; 
    ResMRR = array2table([-1, MRRset], 'VariableNames',['none',strcat('MRR_',methodnames) ]);  
    Restable{rowname,{'table','type','value'}} = {ResMRR,type,'all' }; 
    kk = kk + 1; is_scalar(kk,1) = true; 
     
    % top-k AURecall  CurveRecall_topk          'AURecall' 'AUPrec'    
    for topk=[50, 100 , 200]   
        ind = top_k_vec<=topk; maxK = length( top_k_vec ); 
        AURecall = sum( CurveRecall_topk( ind ,:),1   )  ; 
        type = ['AURecall',num2str(topk)]; rowname = type; 
        ResAU = array2table([topk, AURecall], 'VariableNames',['none',strcat([type,'_'],methodnames) ]);  
        Restable{rowname,{'table','type','value'}} = {ResAU,type,['top1-',num2str(min(maxK,topk))] }; 
        kk = kk + 1; is_scalar(kk,1) = true; 
    end
     
    % top-k AUPrec        CurvePrec_topk 
    for topk=[50, 100 , 200]   
        ind = top_k_vec<=topk; maxK = length( top_k_vec ); 
        AUPrec = sum( CurvePrec_topk( ind ,:),1   )  ; 
        type = ['AUPrec',num2str(topk)]; rowname = type; 
        ResAU = array2table([topk, AUPrec ], 'VariableNames',['none',strcat([type,'_'],methodnames) ]);  
        Restable{rowname,{'table','type','value'}} = {ResAU,type,['top1-',num2str(min(maxK,topk))] }; 
        kk = kk + 1; is_scalar(kk,1) = true; 
    end
     
    type = 'topk_recall'; 
    if ~isempty( ScalarSet ) &&  ismember(type, ScalarSet )
        for ii=1:n_topk_scalar
            strtopk = num2str(topk_scalarset(ii));
             rowname = [type,'_',num2str(ii)] ; 
            tb = array2table([topk_scalarset(ii), ScalarSet_Recall_topk(ii ,:)], 'VariableNames',['topk',strcat('Recall_',methodnames)] );
            Restable{rowname,{'table','type','value'}} = {tb,type,strtopk }; 
            kk = kk + 1; is_scalar(kk,1) = true; 
        end
    end
    %
    type = 'topk_prec';
    if ~isempty( ScalarSet ) && ismember(type, ScalarSet )
        for ii=1:n_topk_scalar
            strtopk = num2str(topk_scalarset(ii));  
            rowname = [type,'_',num2str(ii)] ; 
            tb = array2table([topk_scalarset(ii),ScalarSet_Prec_topk(ii ,:)], 'VariableNames',['topk',strcat('Prec_',methodnames)] )  ;
            Restable{rowname,{'table','type','value'}} = {tb,type,strtopk };  
            kk = kk + 1; is_scalar(kk,1) = true; 
        end  
    end
    
    type = 'topk_recall_p';
    if  ~isempty( ScalarSet ) && ismember(type, ScalarSet )
        for ii=1:n_topk_scalar_percent
            strtopk = num2str(topk_scalarset_percent(ii));
              rowname = [type,'_',num2str(ii)] ; 
            tb = array2table([topk_scalarset_percent(ii), ScalarSet_Recall_topk_percent(ii ,:)], 'VariableNames',['topk_p',strcat('Recall_',methodnames)] )  ;
            Restable{rowname,{'table','type','value'}} = {tb,type,strtopk };    
            kk = kk + 1; is_scalar(kk,1) = true; 
        end
    end
    %
    type = 'topk_prec_p';
    if  ~isempty( ScalarSet ) && ismember(type, ScalarSet )
        for ii=1:n_topk_scalar_percent
            strtopk = num2str(topk_scalarset_percent(ii));  
             rowname = [type,'_',num2str(ii)] ; 
            tb = array2table([topk_scalarset_percent(ii),ScalarSet_Prec_topk_percent(ii ,:)], 'VariableNames',['topk_p',strcat('Prec_',methodnames)] )  ;
            Restable{rowname,{'table','type','value'}} = {tb,type,strtopk };    
            kk = kk + 1; is_scalar(kk,1) = true; 
        end  
    end
     
    %%       
    type = 'CurveTopKvsRecall'; rowname = type;  
    if ~isempty( CurveSet ) && ismember(type, CurveSet )
        ResTopKvsRecall = array2table([top_k_vec, CurveRecall_topk], 'VariableNames',['topk', strcat('Recall_',methodnames) ]  );  
        Restable{rowname,{'table','type','value'}} = {ResTopKvsRecall,type,'all' };    
        kk = kk + 1; is_scalar(kk,1) = false; 
    end
    %   
    type = 'CurveTopKvsPrec'; rowname = type; 
    if ~isempty( CurveSet ) && ismember(type, CurveSet )
        ResTopKvsPrec   = array2table([top_k_vec, CurvePrec_topk  ], 'VariableNames',['topk', strcat('Prec_',methodnames) ]);   
        Restable{rowname,{'table','type','value'}} = {ResTopKvsPrec,type,'all' };   
        kk = kk + 1; is_scalar(kk,1) = false; 
    end 
    type = 'CurveTopKvsRecall_p'; rowname = type; 
    if  ~isempty( CurveSet ) && ismember(type, CurveSet )
        ResTopKvsRecall = array2table([top_k_vec_percent, CurveRecall_topk_percent], 'VariableNames',['topk_p', strcat('Recall_',methodnames) ]  );  
        Restable{rowname,{'table','type','value'}} = {ResTopKvsRecall,type,'all' }; 
        kk = kk + 1; is_scalar(kk,1) = false; 
    end
    %
    type = 'CurveTopKvsPrec_p'; rowname = type; 
    if  ~isempty( CurveSet ) && ismember(type, CurveSet )
        ResTopKvsPrec   = array2table([top_k_vec_percent, CurvePrec_topk_percent  ], 'VariableNames',['topk_p', strcat('Prec_',methodnames) ]);   
        Restable{rowname,{'table','type','value'}} = {ResTopKvsPrec,type,'all' }; 
        kk = kk + 1; is_scalar(kk,1) = false; 
    end
    
    %%  curves
    type = 'ROC'; rowname = type; 
    if  ~isempty( CurveSet ) && ismember(type, CurveSet )
        ResROC = array2table([reshape(XVals,[],1), ROCy],   'VariableNames',['ROCx', strcat('ROCy_',methodnames) ]  );  
        Restable{rowname,{'table','type','value'}} = {ResROC,type,'all' }; 
        kk = kk + 1; is_scalar(kk,1) = false; 
    end
    %
    type = 'PRC'; rowname = type; 
    if  ~isempty( CurveSet ) && ismember(type, CurveSet )
        ResPRC = array2table([reshape(XVals,[],1), PRCy  ], 'VariableNames',['PRCx', strcat('PRCy_',methodnames) ]);   
        Restable{rowname,{'table','type','value'}} = {ResPRC,type,'all' };  
    % %     Restable{rowname,'table'} = {ResPRC}; 
    % %     Restable{rowname,'type'}  = {type}; 
    % %     Restable{rowname,'value'} = {'all'};  
    % %     Restable{rowname,'scalar'}= false;  
        kk = kk + 1; is_scalar(kk,1) = false;  
    % %     Restable{:,'is_scalar'} = is_scalar; 
    end
    Restable.is_scalar = is_scalar; 

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [XVals,YVals ] = getUniformXY(X,Y,XVals, funiform, y0, UniqueMethod , IsExtrap ) 
if ~exist('XVals','var') || isempty( XVals  )  
    XVals = []; 
end

if ~exist('funiform','var') || isempty( funiform  )
    funiform = @mean; 
end

if ~exist('y0','var') || isempty( y0  )
    y0 = [] ; 
end

if ~exist('UniqueMethod','var') || isempty( UniqueMethod  )
    UniqueMethod = 'linearperturbation' ;  
end

if ~exist('IsExtrap','var') || isempty( IsExtrap  )
    IsExtrap = false ;  
end

X= X(:);
Y= Y(:);  

if nargin<3 || strcmp( XVals, 'unique') || isempty( XVals )
    [X_u,~,ic] = unique(X) ; 
    Y_u = zeros( size(X_u) ); 
    for ii = 1: length( X_u )
        % Y_u( ii ) = mean( Y(ic==ii) ) ;
        Y_u( ii ) = funiform( Y(ic==ii) ) ;
    end
    XVals = X_u;
    YVals = Y_u;   
    if nargout>2
        warning('Output1-2 are used for unique X; Output3-4 are unneccesary!');
    end
    return; 

elseif isnumeric(XVals) 
    XVals =XVals(:) ;
    switch UniqueMethod  
        case 'linearperturbation' 
            if  ~( any( X==0 ) ) && ~isempty( y0  )  % set y0 at x=0 if no point of x=0
                % any( X==0 )
                X = [0; X ];
                Y = [y0;Y ];
            end
            % check if not extrap
            if ~IsExtrap
                X_max = max( X ) ;  
                X_min = min( X ) ;    
                if ~( all( X_min<=XVals )  && all( XVals <=X_max ) )
                    error('IsExtrap is false. The data cannot be extraped.');
                end            
            end
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
            X_rescale = max( abs(X) ); 
            if X_rescale~=0 
               if X_rescale ~=1     % normalization for more effective op 
                    X     = X./X_rescale;     
                    XVals = XVals./X_rescale;     
               end             
            else
                error('X_rescale: x_max==0'); 
            end 
            pEPS     = 10^-12 ;
            length_X = length(X); 
            if length_X*pEPS>1/length(XVals)/10  
                warning(['List of X is too long: length(X)=',num2str(length_X)]); 
            elseif length_X<2 
                error(['List of X is too short: length(X)=',num2str(length_X)]); 
            end          
            X = X + [0:length_X-1]'.*pEPS;  
            try 
                if IsExtrap
                    YVals = interp1(X ,Y ,XVals,'linear','extrap') ; 
                else
                    YVals = interp1(X ,Y ,XVals,'linear' ) ;             
                end 
            catch 
                disp( min(X) )
                disp( max(X) )
                disp( min(XVals) )
                disp( max(XVals) )
                [ min(XVals)<min(X)        max(X)<max(XVals)]
            end
            %
            if X_rescale~=1 
                XVals = XVals*X_rescale; 
            end
            return; 

        % % % % % % % % % % % % % % % % % % % % %       
        case 'unique'
            [X_u,~,ic] = unique(X) ; 
            Y_u = zeros( size(X_u) ); 
            for ii = 1: length( X_u )
                Y_u( ii ) = funiform( Y(ic==ii) ) ;
            end
            % % set y0 at x=0 if no point of x=0
            if  ~( any( X_u==0 ) ) && ~isempty( y0  )  
                X_u = [0;X_u ]; 
                Y_u = [y0;Y_u ];
            end
            %
            if IsExtrap
                YVals = interp1(X_u ,Y_u ,XVals,'linear','extrap') ; 
            else
                YVals = interp1(X_u ,Y_u ,XVals,'linear' ) ;             
            end 

            return; 

        otherwise
            error( 'UniqueMethod is wrong. It should be linearperturbation or unique.' );
    end    
    
else 
	error('Input parameter: XVals is wrong for getUniformXY(X,Y,XVals)!') ;     
end
 
end
