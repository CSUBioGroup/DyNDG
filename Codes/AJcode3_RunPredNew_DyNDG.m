runlabel='RunDyDNet-variant-'
fdatestrForfile = datestr(now,'yyyy.mmm.dd-HH.MM.SS') ; 
fdatestr =datestr(now,'yyyy.mmm.dd')
Disease=[{'AML','CML','CLL'}];
Disease_DyNetSet=struct('AML',{'GSE122917_ksigma2'},'CML',{'GSE47927_ksigma2'},'CLL',{'GSE2403_ksigma2'});
 
netnameSet = { 'STRINGa','HumanNetXC'};


for i_Dis=1:length(Disease)
    disease=Disease{i_Dis};
    Root=['Dis_',disease];
    FormatMatData = [Root,filesep,disease,'_FormatMatData' ];
    ResDir        = [Root,filesep,runlabel,fdatestr];  if ~exist(ResDir,'dir'); mkdir(ResDir); end
    DyNetTypeSet = cellstr(Disease_DyNetSet(1).(disease)); 
    for i_DyNetType = 1:length( DyNetTypeSet )
        DyNetType = DyNetTypeSet{i_DyNetType}
        DyNetTypeSplit = strsplit( DyNetType,'_');
        DyExpDataType =  DyNetTypeSplit{1};
        DyNetMethod   = DyNetTypeSplit{2} ;
        fin = []; knn    = []; pre_fndata =[];
        prelog = false ;
        switch DyExpDataType

            case 'GSE122917'
                DIR_gExpDataSet      = 'AML_GSE122917_gExp3Stage';
                MarkStr_gExpDataFile = 'GeneExpDataStage' ;
                logstr = 'log'
                prelog = true
            case 'GSE47927'
                DIR_gExpDataSet      = 'CML_GSE47927_gExp4Stage';
                MarkStr_gExpDataFile = 'GeneExpDataStage' ;
                logstr = 'log'
                prelog = true
            case 'GSE2403'
                DIR_gExpDataSet      = 'CLL_GSE2403_gExp3Stage';
                MarkStr_gExpDataFile = 'GeneExpDataStage' ;
                logstr = 'log'
                prelog = true


            otherwise
        end
        DIR_gExpOrNetDataSet = strjoin({DIR_gExpDataSet,logstr, DyNetMethod},'_');
    	pre_fndata = DIR_gExpOrNetDataSet;
    	%%%%%%
    	fin = [FormatMatData,filesep,'OK-DataSet-',pre_fndata,'-NetMatrixSet_StaticNet&DyNet_.mat' ]
    	load( fin )
    	disset = DataSet.DisGeneSet.Properties.RowNames 

    	% %
    	fdatestr   = datestr(now,'yyyy.mmm.dd-HH.MM.SS') ;
    	resdir = [ResDir ] , if ~exist(resdir,'dir'); mkdir(resdir); end
    	%
    	nCVTimes  = 1,   n_fold    = 1 ,   CVtype = ['PredNewIndTest'] ,  MinSizeDisGeSet = n_fold ;

    	for i_net = 1:length( netnameSet )
    		netname = netnameSet{i_net}
    		fndata = [pre_fndata,'_', netname]
    		if ~isfield(DataSet, fndata);
    			disp(repmat(['*'],10,1)  )
    			error(['There is no fieldname: ', fndata] ) ;
    			continue;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    		Mstatic        = DataSet.(fndata).Mstatic;
    		Mdynetset      = DataSet.(fndata).Mdynetset;
    		TableNode_gene = DataSet.(fndata).TableNode_gene;
    		geneset_ref    = DataSet.(fndata).geneset_ref;
    		%%%
    		filtering = 0
    		if filtering
    			fnset = fieldnames(Mdynetset);
    			ksum  = 0 ;
    			for i_net = 1:length( fnset )
    				Mtp  = Mdynetset.(fnset{i_net});
    				ksum = ksum + sum( Mtp>0,2);
    			end
    			ind = ksum>0;
    			Mstatic       = Mstatic(ind,ind);   netsize = size( Mstatic )
    			for i_net = 1:length( fnset )
    				Mtp  = Mdynetset.(fnset{i_net});
    				Mdynetset.(fnset{i_net}) = Mtp(ind,ind);
    			end
    			TableNode_gene = TableNode_gene(ind,:);
    			geneset_ref    = geneset_ref(ind,:);
    		end
    		netsize        = size( Mstatic )


    		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    		fChrom_ID = 'update2022-3-24_TableGeneLocation.mat'
    		datatemp = load(fChrom_ID); TableGeneLocation = datatemp.TableGeneLocation; datatemp =[];
    		[~,ia,ic]= intersect(TableGeneLocation.symbol, geneset_ref ); TableGeneLocation = TableGeneLocation(ia,:);
    		SizeInTableGenes_genelocation = size(TableGeneLocation)
    		[ mat_genedistance, matlogic_coChrom , Chrom_ID ,GeneStart, GeneEnd,geneset00   ]  = getGeneChromDistance2022(TableGeneLocation, [],[], true ) ;
    		%
    		mat_genedistance =   ( getMap2RefMatrix( mat_genedistance,geneset00 ,geneset00, geneset_ref,geneset_ref  ) ) ;
    		matlogic_coChrom =   ( getMap2RefMatrix( matlogic_coChrom,geneset00 ,geneset00, geneset_ref,geneset_ref  ) ) ;
    		is_sparse_mat_genedistance = issparse(mat_genedistance) ;
    		n_neighbors_artificial_linkage_interval = 99  ;

    		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    		% % % % % % % % % % % % % %
    		n_gene_all = length( geneset_ref )
    		n_disease  = 1

    		%
    		SubNetworkNames =  fieldnames(Mdynetset )';
    		disp( SubNetworkNames  )
    		%
    		%
    		ind_disset = 1: length( disset )
    		for ii=ind_disset
    			DisSetName = disset{ii}
    			ac_gene_dis00_train  = DataSet.(fndata).DisGeneVecSet{DisSetName,1}{1};
                PStage0=DataSet.(fndata).DEGVecSet;
    			if filtering; ac_gene_dis00_train = ac_gene_dis00_train(ind); end
    			n_disgene = nnz( ac_gene_dis00_train )

    			if n_disgene<MinSizeDisGeSet; continue; end
    			gScorelistTB_discell = cell(1,1 );

                i_cv=1; i_fold = 1 ;

                disp(['i_cv-',num2str(i_cv),'; i_fold-',num2str(i_fold)])
                P0                = ac_gene_dis00_train;
                P0_G  = P0;
                AdjGfD=P0; ID_dis= 1;

                tic
                TableScores = table ;

                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                MatrixSet_gene2gene_t         = Mdynetset;
                MatrixSet_gene2gene_t.Mstatic = Mstatic;
                plusSN = '_SN'  ;

    			if  00
    				InterlayerType = 'LINE';
    				for restart=0.1:0.2:0.9
        				for pro_jump_DyNet =0:0.2:1
            				for pro_jump_Static=0:0.2:1
            					TableScores_t = A_DGP_DyNDG(Mdynetset,Mstatic, PStage0,P0, restart, pro_jump_DyNet, pro_jump_Static, tau_Static, InterlayerType );
            					parastr =['r_',num2str(restart),'_PJd_',num2str(pro_jump_DyNet),'_PJs_',num2str(pro_jump_Static)];
            					parastr = strrep(parastr,'.','p');
            					TableScores_t.Properties.VariableNames =  strcat( [parastr,'_'], TableScores_t.Properties.VariableNames );
            					TableScores = [TableScores, TableScores_t(:,end)] ;
            				end
        				end
    				end

    			end
    			if  111
    				restart = 0.1; pro_jump_DyNet=1;  pro_jump_Static=0.5;  tau_Static=0.5;
                    stage=size(PStage0);
                    P0=[];
                    for i=1:stage
                        P0=[P0,PStage0{i,1}{1}./sum(PStage0{i,1}{1})];
                    end
    		        P0= mean(P0,2);
                    InterlayerType = 'LINE';
                    TableScores_t = []; 
                    TableScores_t = A_DGP_DyNDG(Mdynetset,Mstatic,PStage0, P0, restart, pro_jump_DyNet, pro_jump_Static, tau_Static, InterlayerType );
    				TableScores_t.Properties.VariableNames =  strcat( TableScores_t.Properties.VariableNames, '_SetTmean');
                    TableScores = [TableScores, TableScores_t] ;
                    Mdynetset_t =[]; Mdynetset_t.eye = speye( size(Mstatic) ) ;
    				InterlayerType = 'LINE';
                    TableScores_t = []; 
    				TableScores_t = A_DGP_DyNDG([],Mstatic,PStage0, P0, restart, pro_jump_DyNet, pro_jump_Static, tau_Static, InterlayerType );
    				TableScores_t.Properties.VariableNames =  strcat( TableScores_t.Properties.VariableNames, '_StaticSetTmean');
    				TableScores = [TableScores, TableScores_t] ;
                    Mstatic_t = speye( size(Mstatic) ) ; 
				    InterlayerType = 'LINE';
				    TableScores_t = A_DGP_DyNDG(Mdynetset,[],PStage0, P0, restart, pro_jump_DyNet, pro_jump_Static, tau_Static, InterlayerType ); 
				    TableScores_t.Properties.VariableNames =  strcat( TableScores_t.Properties.VariableNames, '_DyNetSetTmean');
				    TableScores = [TableScores, TableScores_t] ;   
                end
                if  111      
					plus_method_set = { 'RWR' , 'KS' } ; 
                    stage=size(PStage0);
                    P0=[];
                    for i=1:stage
                        P0=[P0,PStage0{i,1}{1}./sum(PStage0{i,1}{1})];
                    end
    		        P0= mean(P0,2);
                    AdjGfD=P0;
					[TableScores_t ] =  A_DGP_H_methodset(Mstatic,AdjGfD,[], ID_dis, plus_method_set  ) ;       
					TableScores = [ TableScores, TableScores_t ] ;		
                    P0(:) = 1;
					[TableScores_RWRM1 ] = A_RWRMcp1(MatrixSet_gene2gene_t, P0, [],[], 'col') ; 
					TableScores.(['RWRMP1',plusSN]) = TableScores_RWRM1.RWRM1;   
					[TableScores_RWRM2 ] = A_RWRMcp2(MatrixSet_gene2gene_t, P0, [],[], 'col') ; 
					TableScores.(['RWRMP2',plusSN]) = TableScores_RWRM2.RWRM2;  
					 [Pt ] = A_RWRplus_Multigraph(MatrixSet_gene2gene_t, P0, 0.7, [], [], [], 'ColNormalization_MultiGraph_LI'); 
					 TableScores.(['RWRMG',plusSN]) = Pt  ;  Pt =[];  
					 Score_Method=[];AggregationMethod=[]; Score_Method.name = 'RWR'; AggregationMethod.name = 'DRS';
					 [ Pagg ] = A_DGP_AggregationMethodset_MultipleGeneNets_plus(MatrixSet_gene2gene_t, P0, Score_Method,AggregationMethod );
					 TableScores.(['RWRDRS',plusSN]) = Pagg  ;  Pagg =[];
					 Score_Method=[];AggregationMethod=[];Score_Method.name = 'Endeavour' ; AggregationMethod.name = 'NDOS' ;
					 [ Pagg ] = A_DGP_AggregationMethodset_MultipleGeneNets_plus(MatrixSet_gene2gene_t, P0, Score_Method,AggregationMethod );
					 TableScores.(['Endeavour',plusSN]) = Pagg  ;  Pagg =[];
					disp('CP  OK')  ; toc 
                end 
    			toc
    			% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    			ii_dis = 1  ;
    			n_disease_in_Table = 1
    			gScorelistTB_discell{ii_dis,1}  = TableScores ;
    			% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    			methodset = TableScores.Properties.VariableNames
    			methodset_cell  =repmat({[]},n_disease_in_Table, 1);
    			methodset_cell{ii_dis} = methodset;
                %%%%%
    			for jj = ind_disset
    				%%%
    				DisSetName_TEST = disset{jj}
    				ac_gene_dis00_TEST  = DataSet.(fndata).DisGeneVecSet{DisSetName_TEST,1}{1};
    				if filtering; ac_gene_dis00_TEST = ac_gene_dis00_TEST(ind); end
    				%
    				n_fold_real   = n_fold ;
    				streff = '';
    				idx_pos_test  = find( ac_gene_dis00_TEST & ~ac_gene_dis00_train );
    				if isempty( idx_pos_test ); streff ='TestSetEmptyIneff'; idx_pos_test= find(  ac_gene_dis00_train ); ; end
    				n_pos_test =length(idx_pos_test);
    				%
    				ac_gene_dis00 = ac_gene_dis00_train | ac_gene_dis00_TEST ;
    				idx_pos       = find( ac_gene_dis00 );
    				idx_neg       = find( ~ac_gene_dis00 );
    				%%%%
    				TableDiseases_cut = table;
    				TableDiseases_cut{disease,'Dis'}={disease};

    				% % % % % % % % % % % % % % % % % % % % % % % % % % % %
    				ScalarSet = {'AUROC','AUPRC','topk_recall','topk_prec'};
    				CurveSet = {'CurveTopKvsRecall','CurveTopKvsPrec' }  ;

    				% % % % % % % % % % % % % % % % % % % % % %
    				AUCmeanset   = table;
    				AUPRCmeanset = table;
    				% % % % % % % % % % % % % % % %
    				n_disease_in_Table = 1
    				ResTablePerfCurve_AllDiseaseMeanCV     = repmat({table},n_disease_in_Table, 1);
    				ResTablePerfCurve_AllDiseaseMeanCV_std = repmat({table},n_disease_in_Table, 1);
    				%
    				ResTablePerfCurve_AllDiseaseMeanCV_X100_RC     = repmat({table},n_disease_in_Table, 1);
    				ResTablePerfCurve_AllDiseaseMeanCV_X100_RC_std = repmat({table},n_disease_in_Table, 1);
    				%
    				ResTablePerfCurve_AllDiseaseMeanCV_X100_ALI     = repmat({table},n_disease_in_Table, 1);
    				ResTablePerfCurve_AllDiseaseMeanCV_X100_ALI_std = repmat({table},n_disease_in_Table, 1);
    				%
    				nCV_list   = zeros( n_disease_in_Table, 1 );
    				tAUROC = table;   tAUPRC = table;
    				tAUROC_X100_RC =[table]; tAUPRC_X100_RC =[table];
    				tAUROC_X100_ALI=[table]; tAUPRC_X100_ALI=[table];

    				%

    				tic
    				mat_genedistance_copy = mat_genedistance;
    				matlogic_coChrom_copy = matlogic_coChrom;
    				ResTablePerfCurve_1Disease_AllCV          =repmat({table},n_fold_real, nCVTimes);
    				ResTablePerfCurve_1Disease_AllCV_X100_RC  =repmat({table},n_fold_real, nCVTimes);
    				ResTablePerfCurve_1Disease_AllCV_X100_ALI =repmat({table},n_fold_real, nCVTimes);

    				%
    				idx_neg_test_WG   = idx_neg ;  n_neg_test_all = length(  idx_neg_test_WG  ) ; 
    				idx_neg_test_X100_RC  = idx_neg_test_WG(   randperm(n_neg_test_all,  min(n_pos_test*n_neighbors_artificial_linkage_interval,n_neg_test_all)  ) ); n_neg_test_x100 = length(  idx_neg_test_X100_RC  ) ;
    				idx_neg_test_X100_ALI = getGeneNeigbors_ArtitificalLinkageInterval(mat_genedistance_copy,matlogic_coChrom_copy,idx_pos_test, idx_pos, n_neighbors_artificial_linkage_interval   ) ;
    				%
    				idx_test_pos_neg_WG        = [idx_neg_test_WG; idx_pos_test ] ;
    				idx_test_pos_neg_X100_RC   = [idx_neg_test_X100_RC; idx_pos_test ] ;
    				idx_test_pos_neg_X100_ALI  = union( idx_neg_test_X100_ALI, idx_pos_test ) ;
    				%
    				[Restable, n_gene,n_method] = getResPerfFromScorelistNew2(ac_gene_dis00(idx_test_pos_neg_WG),TableScores(idx_test_pos_neg_WG,:)   , ScalarSet, CurveSet  );
    				ResTablePerfCurve_1Disease_AllCV{i_fold, i_cv} = Restable  ;
    				[RestableX100_RC, n_gene,n_method] = getResPerfFromScorelistNew2(ac_gene_dis00(idx_test_pos_neg_X100_RC),TableScores(idx_test_pos_neg_X100_RC,:)   , ScalarSet, CurveSet  );
    				ResTablePerfCurve_1Disease_AllCV_X100_RC{i_fold, i_cv} = RestableX100_RC  ;

   				[RestableX100_ALI, n_gene,n_method] = getResPerfFromScorelistNew2(ac_gene_dis00(idx_test_pos_neg_X100_ALI),TableScores(idx_test_pos_neg_X100_ALI,:)   , ScalarSet, CurveSet  );
                ResTablePerfCurve_1Disease_AllCV_X100_ALI{i_fold, i_cv} = RestableX100_ALI  ;


                aaa = [Restable{'AUPRC',1}{1} ; RestableX100_RC{'AUPRC',1}{1}; RestableX100_ALI{'AUPRC',1}{1}  ] ;
                aaa{:,:}
                disp(DyNetType); disp(netname);
                % save mean and std
                [meanCVTable, stdCVTable, nCV] = getResMeanFromTableSetNew(  ResTablePerfCurve_1Disease_AllCV  ) ;
                nCV_list(ii_dis) = nCV;
                ResTablePerfCurve_AllDiseaseMeanCV{ii_dis,1}     = meanCVTable; 
                ResTablePerfCurve_AllDiseaseMeanCV_std{ii_dis,1} = stdCVTable;  
                tAUROC = [tAUROC; meanCVTable{'AUROC',1}{1} ]
                tAUPRC = [tAUPRC; meanCVTable{'AUPRC',1}{1} ]

                [meanCVTable, stdCVTable, nCV] = getResMeanFromTableSetNew(  ResTablePerfCurve_1Disease_AllCV_X100_RC  ) ;
                nCV_list(ii_dis) = nCV;
                ResTablePerfCurve_AllDiseaseMeanCV_X100_RC{ii_dis,1}     = meanCVTable;  
                ResTablePerfCurve_AllDiseaseMeanCV_X100_RC_std{ii_dis,1} = stdCVTable;   
                tAUROC_X100_RC = [tAUROC_X100_RC;  meanCVTable{'AUROC',1}{1} ]
                tAUPRC_X100_RC = [tAUPRC_X100_RC; meanCVTable{'AUPRC',1}{1}  ]

                [meanCVTable, stdCVTable, nCV] = getResMeanFromTableSetNew(  ResTablePerfCurve_1Disease_AllCV_X100_ALI  ) ;
                nCV_list(ii_dis) = nCV;
                ResTablePerfCurve_AllDiseaseMeanCV_X100_ALI{ii_dis,1}     = meanCVTable; 
                ResTablePerfCurve_AllDiseaseMeanCV_X100_ALI_std{ii_dis,1} = stdCVTable;   
                tAUROC_X100_ALI = [tAUROC_X100_ALI;  meanCVTable{'AUROC',1}{1} ]
                tAUPRC_X100_ALI = [tAUPRC_X100_ALI; meanCVTable{'AUPRC',1}{1}  ]

                %
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                methodset = methodset_cell{1};
                %
                [ResMeanDis,ResVecAllDis,ResScalarAllDis, n_Tables] = getResFinalOutputForm(  TableDiseases_cut, ResTablePerfCurve_AllDiseaseMeanCV_X100_ALI, ResTablePerfCurve_AllDiseaseMeanCV_X100_ALI_std, nCV_list, methodset  ) ;
                Res_X100_ALI=[];
                Res_X100_ALI.ResMeanDis      = ResMeanDis   ;
                Res_X100_ALI.ResVecAllDis    = ResVecAllDis   ;
                Res_X100_ALI.ResScalarAllDis = ResScalarAllDis ;
                Res_X100_ALI.n_Tables        = n_Tables  ;
                Res_X100_ALI.DisSetName = DisSetName  ;
                Res_X100_ALI.methodset       = methodset  ;
                Res_X100_ALI.nCV_list        = nCV_list  ;
                Res_X100_ALI.date            = datestr(now,'yyyy.mmm.dd-HH.MM.SS')  ;
                Res_X100_ALI.note            = 'Artificial linkage interval: select 99 neighboring genes for each of pos-test genes, plus themselves, resulting a test set with X100 genes' ;
                %
                [ResMeanDis,ResVecAllDis,ResScalarAllDis, n_Tables] = getResFinalOutputForm(  TableDiseases_cut, ResTablePerfCurve_AllDiseaseMeanCV_X100_RC, ResTablePerfCurve_AllDiseaseMeanCV_X100_RC_std, nCV_list, methodset  ) ;
                Res_X100_RC=[];
                Res_X100_RC.ResMeanDis      = ResMeanDis   ;
                Res_X100_RC.ResVecAllDis    = ResVecAllDis   ;
                Res_X100_RC.ResScalarAllDis = ResScalarAllDis ;
                Res_X100_RC.n_Tables        = n_Tables  ;
                Res_X100_RC.DisSetName = DisSetName  ;   ;
                Res_X100_RC.methodset       = methodset  ;
                Res_X100_RC.nCV_list        = nCV_list  ;
                Res_X100_RC.date            = datestr(now,'yyyy.mmm.dd-HH.MM.SS')  ;
                Res_X100_RC.note            = 'Random control: #of genes in control gene set equal to #of pos-test genes X100' ;
                %
                [ResMeanDis,ResVecAllDis,ResScalarAllDis, n_Tables] = getResFinalOutputForm(  TableDiseases_cut, ResTablePerfCurve_AllDiseaseMeanCV, ResTablePerfCurve_AllDiseaseMeanCV_std, nCV_list, methodset  ) ;
                Res_WG=[];
                Res_WG.ResMeanDis      = ResMeanDis   ;
                Res_WG.ResVecAllDis    = ResVecAllDis   ;
                Res_WG.ResScalarAllDis = ResScalarAllDis ;
                Res_WG.n_Tables        = n_Tables  ;
                Res_WG.DisSetName = DisSetName  ;   ;
                Res_WG.methodset       = methodset  ;
                Res_WG.nCV_list        = nCV_list  ;
                Res_WG.date            = datestr(now,'yyyy.mmm.dd-HH.MM.SS')  ;
                Res_WG.note            = 'Whole genome control: All unknown genes are used as control gene set' ;
                para=[];
                para.fin    = fin;
                para.netname    = netname;
                para.DisSetName = DisSetName;
                para.SubNetworkNames = SubNetworkNames;
                para.knn        = knn;
                para.Date       = fdatestrForfile;
                para.MinSizeDisGeSet     =  MinSizeDisGeSet ;
                para.n_disgene     =  n_disgene ;
                para.n_gene     =  n_gene_all ;
                para.n_disease    =  n_disease ;
                para.CVtype    =  CVtype ;
                para.CVtime    =  nCVTimes ;
                para.n_neighbors_artificial_linkage_interval    =  n_neighbors_artificial_linkage_interval ;
                para.ScalarSet    =  ScalarSet ;
                para.CurveSet    =  CurveSet ;
                date_complete = datestr(now,'yyyy.mmm.dd-HH.MM.SS')  ;
                parastr = sprintf('CVtype=%s_CVtime=%d_MSDGS%d',CVtype  ,  nCVTimes, MinSizeDisGeSet )
                outfile = [resdir,filesep, DyNetType,'-DyNet-Res3Ctrl_DGP_',netname,'_',DisSetName,'-Pred-',DisSetName_TEST,'-',streff,'.mat']
                save([outfile],  'Res_X100_ALI','Res_X100_RC','Res_WG','para', 'date_complete' )   ;
                disp(['Complete: ',outfile ])
                disp([ '             ' ,date_complete])

    			end
    			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                dir_RES =['RES-',disease,'-DyDNet-PredAll-',fdatestrForfile]; if ~exist(dir_RES,'dir'); mkdir(dir_RES);end
                TableScores_t = gScorelistTB_discell{1,1} ;
                MethodNameSet = TableScores_t.Properties.VariableNames; % method info
                disIDstr_cut = TableDiseases_cut.Properties.RowNames;
                %
                n_method = length( MethodNameSet )
                n_dis    = length(gScorelistTB_discell)
                mat_gscore_dis_method = cell(n_dis,1);
                for ii_dis = 1:n_dis
                    TB_gscore_method =  gScorelistTB_discell{ii_dis};
                    IDstr_dis = disIDstr_cut{ii_dis} ;
                    mat_gscore_dis_method{ii_dis} = reshape( full( TB_gscore_method{:,:} ), [],1, n_method) ;
                end
                disp('Complete preprocessing of score list----1')
                mat_gscore_dis_method = cat(2,mat_gscore_dis_method{:} );
                disp('Complete preprocessing of score list----2')
                % % %
                AdjGfD_effdis = ac_gene_dis00_train;
                %

                ind_pos = find(AdjGfD_effdis); AdjGfD_eff_infvalues = double( AdjGfD_effdis ); AdjGfD_eff_infvalues(ind_pos) = -inf;
                mat_gscore_dis_method = mat_gscore_dis_method  +  AdjGfD_eff_infvalues  ;
                sorttype              = 'descend'
                OP_IdenticalValues    ='MeanRank'
    			% %     OP_IdenticalValues ='none';
                [mat_Rank_dis_method, ord ]= getRankingOfScoreList(   mat_gscore_dis_method, sorttype, OP_IdenticalValues );
                SZ_XYZ=size( mat_gscore_dis_method )
                n_gene =size( mat_gscore_dis_method, 1 )
    			% % %
                methodnames_formal = MethodNameSet
                TBmethod_scoremat_gene_dis = table;
                mat_gene_dis = zeros(SZ_XYZ(1), SZ_XYZ(2) );
                for ii_method=1:n_method

                    mat_gene_dis = mat_Rank_dis_method(:,:,ii_method);
                    mat_gene_dis(ind_pos) = 0 ;
                    mat_gene_dis = 1-mat_gene_dis./(SZ_XYZ(1)- sum(AdjGfD_effdis>0,1) +eps ) ;  
                    %
                    dir_RES_method =[dir_RES,filesep,MethodNameSet{ii_method}]; if ~exist(dir_RES_method,'dir'); mkdir(dir_RES_method);end
                    digitss = ceil( log10( n_gene ) ) +2;

                    for jj_dis=1:n_dis
                        tb = table;
                        tb.gene = geneset_ref(ord(:,jj_dis,ii_method));
                        %%% geneset_ref
                        score = round( mat_gene_dis(:,jj_dis),  digitss   );
                        tb.score = score(ord(:,jj_dis,ii_method));      
                        fout = [dir_RES_method,filesep, DisSetName,'-',DyNetType,'-',netname,'-',disIDstr_cut{jj_dis},'.txt'] ;
                        writetable(tb,fout,'WriteVariableNames', false, 'FileType', 'text' ,'delimiter','\t')
                    end
                end
                disp(['Complete saving seperate files!!!!!!!!!!!!!!!!!!!!!!!!!'])



    		end
        end
    end
end