Disease=[{'AML','CML','CLL'}];
DIR_StaticPPI       = 'StaticPPI' ;
Disease_DyNetSet=struct('AML',{'GSE122917'},'CML',{'GSE47927'},'CLL',{'GSE2403'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tb_fStaticPPIset = table;   %%%  file, score_threshold
tb_fStaticPPIset{'Science2015',1:2} = {[DIR_StaticPPI,filesep, 'Science2015_static_network_PPI.txt'  ], []   };
tb_fStaticPPIset{'STRING',1:2}      = {[DIR_StaticPPI,filesep, 'STRING_static_network_PPI.txt'  ], [0.4] };
tb_fStaticPPIset{'HPRD',1:2}        = {[DIR_StaticPPI,filesep, 'HPRD_static_network_PPI.txt'  ] , []};
tb_fStaticPPIset{'BioGrid',1:2}     = {[DIR_StaticPPI,filesep, 'BIOGRID_static_network_PPI.txt'  ] , []};
tb_fStaticPPIset{'HIPPIE',1:2}     = {[DIR_StaticPPI,filesep, 'HIPPIE_static_network_PPI.txt'  ] , []};
tb_fStaticPPIset{'HumanNetXC',1:2}     = {[DIR_StaticPPI,filesep, 'HumanNetXC_static_network_PPI.txt'  ] , []};
tb_fStaticPPIset{'HumanNetFN',1:2}     = {[DIR_StaticPPI,filesep, 'HumanNetFN_static_network_PPI.txt'  ] , []};
tb_fStaticPPIset{'HumanNetPI',1:2}     = {[DIR_StaticPPI,filesep, 'HumanNetPI_static_network_PPI.txt'  ] , []};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for i_Dis=1:length(Disease)
    DisName=Disease{i_Dis};%要处理的疾病
    Root=['Dis_',DisName];
    DIRdisgene          = sprintf('%s_DisGeneSet', DisName ), MarkStr_DisGeneSetFile  = sprintf('%s_DisGeneSet' , DisName )
    DIR_FormatMatData   = sprintf('%s_FormatMatData', DisName) , if ~exist([Root,filesep,DIR_FormatMatData],'dir'); mkdir([Root,filesep,DIR_FormatMatData]); end
    %%%% generate multi-stage networks
    logstr = 'log'
    DyExpDataTypeSet = cellstr(Disease_DyNetSet(1).(DisName)); %选择疾病对应想要处理的GSE
    DyNetMethod = 'ksigma2'
    tbError = table;

    for i_DyExpDataType = 1:length( DyExpDataTypeSet )
    	DyExpDataType = DyExpDataTypeSet{i_DyExpDataType};
    	switch DyExpDataType
            case 'GSE122917'
                DIR_gExpDataSet      = 'AML_GSE122917_gExp3Stage'; 
			    MarkStr_gExpDataFile = 'GeneExpDataStage' ;
    		case 'GSE47927'
    			DIR_gExpDataSet      = 'CML_GSE47927_gExp4Stage';
    			MarkStr_gExpDataFile = 'GeneExpDataStage' ;
            case 'GSE2403'
                DIR_gExpDataSet      = 'CLL_GSE2403_gExp3Stage'; 
			    MarkStr_gExpDataFile = 'GeneExpDataStage' ;

    		otherwise
    			error(['!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'])
    	end
    	%%%%%%
    	fprintf( ['\n\n\n'])
    	if ~exist([Root,filesep,DIR_gExpDataSet],'dir'); error( ['There is no directory: ' ,DIR_gExpDataSet] ); end

    	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    	% % % % % % % % % % % % % % % % % % % % % % % % % %
    	% load multistage expdata    必须确保 每个表达数据文件的 行基因顺序都是一样的，否则还要重新对其表达数据的行
    	disp( repmat('*',1, 50) ) %打印50个*
    	disp( [ 'Load multistage expdata from ', Root,filesep,DIR_gExpDataSet,filesep, MarkStr_gExpDataFile,'*'] )%打印
    	fgExpDataSet       = dir( [Root,filesep,DIR_gExpDataSet,filesep,[MarkStr_gExpDataFile,'*'] ] );%拼接地址，获得指定文件夹下的所有文件
    	struct2table(fgExpDataSet) %将结构化数组转化为table
    	n_stage            = size(fgExpDataSet,1) ;%返回其行数即有几个stage
    	fDynNetSet         = []  ;
    	gsExp_set          = [] ;
        GeneSymbol_set     = [] ;
    	for i_stage = 1:n_stage
    		fgExpDataSet   = dir( [Root,filesep,DIR_gExpDataSet,filesep,MarkStr_gExpDataFile,num2str(i_stage),'*' ]  );%获取某stage到结构化数组中
    		disp(['Loading ',DIR_gExpDataSet,'--i_stage=', num2str(i_stage)])
    		%
    		if size(fgExpDataSet,1)>1;
    			disp(fgExpDataSet);
    			error('There are more than 1 network at this stage.');
    		end
    		%
    		fDynNetSet{i_stage,1} = [Root,filesep,DIR_gExpDataSet,filesep,fgExpDataSet.name];%fgExpDataSet.name此时读取的文件名
    		gsExp_set{i_stage,1}  = readtable(fDynNetSet{i_stage,1},  'ReadVariableNames',true, 'ReadRowNames', true, 'delimiter','\t'  );
            GeneSymbol_set{i_stage,1}=gsExp_set{i_stage,1}.Properties.RowNames;
    	end
    	fprintf( ['Load--OK \n\n\n'] )
    	%
    	% scheme 1    median
    	% scheme 2    all
    	disp( ['Pre-Process multistage expdata ' ] )
    	gsExp_MatrixSet = [];
    	for i_stage = 1:n_stage
    		disp(['Pre-Processing ', DIR_gExpDataSet,'--i_stage=',num2str(i_stage)])
    		tb = gsExp_set{i_stage,1}; %tb是该i_stage阶段的样本原始矩阵
            Mexp = (tb{:,:}) ;  %不做log变换
    		Mexp(isnan(Mexp)) = 0;
    		Mexp(isinf(Mexp)) = 0;
    		switch  DyNetMethod
    			case {'3sigma','3sigmaA'}
    				Mexp = median( Mexp, 2 );   %求同一阶段同一基因不同样本的中位数，Mexp是该中位数列表（只有一列）
    			case { '3sigmaB', '3sigma2','ksigma2'}
    				% Mexp = Mexp;
    			otherwise
    				error(['There is no ',DyNetMethod]);
    		end
    		gsExp_MatrixSet{ i_stage, 1} = [Mexp] ;
    	end
    	gsExp_Matrix = cat(2, gsExp_MatrixSet{:} );
    	gExp_mean = mean( gsExp_Matrix, 2 ) ;%基因对应的表达的均值
    	gExp_std  = std( gsExp_Matrix,0, 2 ) ;%基因对应的表达标准差
        gExp_threshold3=gExp_mean+3.*gExp_std.*(1-1./(1+gExp_std.*gExp_std));%基因对应的表达阈值--k=3
        gExp_threshold2=gExp_mean+2.*gExp_std.*(1-1./(1+gExp_std.*gExp_std));%--k=2
        gExp_threshold1=gExp_mean+1.*gExp_std.*(1-1./(1+gExp_std.*gExp_std));%--k=1
    	fprintf( ['Pre-Processing--OK \n\n\n'] )

    	% generate DyNet
    	NetStastics = table  ;
    	StaticPPIsetNames  = tb_fStaticPPIset.Properties.RowNames;%PPI网络的名字们
    	for i_PPI = 1: length( StaticPPIsetNames )
    		StaticPPIname  = StaticPPIsetNames{i_PPI};
    		disp([ repmat('*',1, 50) ])
    		disp(['Generate DyNet by ',StaticPPIname ,' & ', DIR_gExpDataSet])
    		disp(['Loading ',StaticPPIname])
    		% load static network into table % % % % % % % % % % % % % % % % %
    		fStaticPPI   = tb_fStaticPPIset{i_PPI,1}{1};  %当前的PPI网络的地址
    		threshold    = tb_fStaticPPIset{i_PPI,2}{1};  %当前的PPI网络的阈值，例如STRING的0.4
    		tbStaticPPI  = readtable(fStaticPPI,'ReadVariableNames',false);  %读取当前PPI网络：基因1 基因2 边权重
            n_gene = length( unique( tbStaticPPI{:,[1,2]} )   ); %当前PPI网络的节点数
    		n_edge = size(tbStaticPPI,1); %当前PPI网络的边数
    		NetStastics{[ StaticPPIname,'-staticnet',],{'N_node','N_edge','km'}} =[ n_gene , n_edge, n_edge*2/n_gene ] ;

    		%
    		DIRout =[Root,filesep, strjoin({DIR_gExpDataSet,logstr,DyNetMethod},'_') ]; if ~exist(DIRout,'dir'); mkdir(DIRout); end
    		for i_stage = 1:n_stage
    			fprintf(['Processing ',StaticPPIname, ' & ' DIR_gExpDataSet,'--i_stage=',num2str(i_stage), '   ',DyNetMethod, ' \n'])
    			disp(DyNetMethod)
    			switch  DyNetMethod
    				case  {'3sigma','3sigmaA'}
    					gsExp_stage   = gsExp_Matrix(:,i_stage);gExp_threshold=gExp_threshold3;
    					ind_gsExp_com_threshold=(gsExp_stage>gExp_threshold);
    					genes_active  = GeneSymbol_set{i_stage,1}(ind_gsExp_com_threshold) ;
                        %存储表达值大于阈值的基因名
    					tbDyPPI       = tbStaticPPI;
    					[uGenes, ia, ib] = unique( [tbDyPPI{:,1}; tbDyPPI{:,2}] ) ;
                        %[tbDyPPI{:,1}; tbDyPPI{:,2}]：tbStaticPPI一二列合并在一起
                        % ia：uGenes在后者的位置；
                        % ib：后者中的在uGenes中的位置
                        %uGenes：在staticPPI的唯一基因列表
    					[~, ic, id] = intersect( uGenes, genes_active );
                        %ic：staticPPI和活跃表达基因的交集在uGenes的位置
                        %id：staticPPI和活跃表达基因的交集在活跃基因的位置
    					uIND        = false( size(uGenes) );
    					uIND(ic)    = 1;
    					ind2        = uIND(ib)  ;  ind2  = reshape( ind2,[],2) ;
    					%% ind_exist   = ind2(:,1) & ind2(:,2) ;
    					ind_exist   = all( ind2 ,2);

    					% % all( ind_exist==ind_exist2)
    					tbDyPPI = tbDyPPI(ind_exist,:);
    				case  { '3sigmaB', '3sigma2'}
    					tbDyPPI       = tbStaticPPI;
    					[uGenes, ia, ib] = unique( [tbDyPPI{:,1}; tbDyPPI{:,2}] ) ;
                        %[tbDyPPI{:,1}; tbDyPPI{:,2}]：tbStaticPPI一二列合并在一起
                        % ia：uGenes在后者的位置；
                        % ib：后者中的在uGenes中的位置
                        %uGenes：在staticPPI的唯一基因列表
    					gsExp_stage = gsExp_MatrixSet{ i_stage, 1} ;
                        %该阶段的表达值
    					n_sample    = size(gsExp_stage,2) ;
                        %
    					scores      = zeros(n_edge,1 ); gExp_threshold=gExp_threshold3;
    					for i_sp = 1: n_sample
    						ind_gsExp_com_threshold=(gsExp_stage(:,i_sp)>gExp_threshold);
    						genes_active  = GeneSymbol_set{i_stage,1}(ind_gsExp_com_threshold) ;
    						[~, ic, id] = intersect( uGenes, genes_active );
    						uIND        = false( size(uGenes) );
    						uIND(ic)    = 1;
    						ind2        = uIND(ib)  ;  ind2  = reshape( ind2,[],2) ;
    						%% ind_exist   = ind2(:,1) & ind2(:,2) ;
    						ind_exist   = all( ind2 ,2);
    						scores(ind_exist) = scores(ind_exist)+tbDyPPI{ind_exist, 3};
    					end
    					tbDyPPI{:,3} = scores./n_sample;
    					tbDyPPI      = tbDyPPI(scores>0,:);
    				case {'ksigma2'}
                        %表达值大于阈值1的即存在，具体的概率体现在边的权重上
                        %计算每个节点的邻居数目，计算出边的权重*原来的权重
                        tbDyPPI =tbStaticPPI;
                        [uGenes,ia,ib]=unique([tbDyPPI{:,1};tbDyPPI{:,2}]);
                        gsExp_stage=gsExp_MatrixSet{i_stage,1};
                        n_sample=size(gsExp_stage,2);
                        scores=zeros(n_edge,1);
                        for i_sp=1:n_sample
                            ind_gsExp_com_threshold=(gExp_threshold1<gsExp_stage(:,i_sp)<gExp_threshold2);
                            genes_active1=GeneSymbol_set{i_stage,1}(ind_gsExp_com_threshold);
                            [~,ic,id]=intersect(uGenes,genes_active1);
                            uIND=zeros(size(uGenes));
                            uIND(ic)=0.68;
                            ind_gsExp_com_threshold=(gExp_threshold2<gsExp_stage(:,i_sp)<gExp_threshold3);
                            genes_active2=GeneSymbol_set{i_stage,1}(ind_gsExp_com_threshold);
                            [~,ic,id]=intersect(uGenes,genes_active2);
                            uIND(ic)=0.95;
                            ind_gsExp_com_threshold=(gExp_threshold3<gsExp_stage(:,i_sp));
                            genes_active3=GeneSymbol_set{i_stage,1}(ind_gsExp_com_threshold);
                            [~,ic,id]=intersect(uGenes,genes_active3);
                            uIND(ic)=0.99;
                            ind2=uIND(ib);ind2=reshape(ind2,[],2);%存着边的2个节点的概率值
                            ind_exist=all(ind2,2);
                            scores(ind_exist)=scores(ind_exist)+tbDyPPI{ind_exist,3}.*ind2(ind_exist,1).*ind2(ind_exist,2);
                        end
                        tbDyPPI{:,3}=scores./n_sample;
                        tbDyPPI=tbDyPPI(scores>0,:);
    				otherwise; error('!!!!!!!!!!!!!!!!!!!');
    			end

    			fout = [DIRout, filesep, StaticPPIname,'_net_Stage',num2str(i_stage) ];
    			writetable(tbDyPPI , fout, 'WriteVariableNames',false, 'WriteRowNames', false, 'delimiter','\t' )
    			%%%fprintf( ['Processing--OK--i_stage=',num2str(i_stage),'\n\n\n'] )
    			nd_gene = length( unique( tbDyPPI{:,[1,2]} )   );
    			nd_edge = size(tbDyPPI,1);
    			NetStastics{[ StaticPPIname,'-',DIR_gExpDataSet,'-dynet-',num2str(i_stage)],{'N_node','N_edge','km'}} =[ nd_gene , nd_edge, nd_edge*2/nd_gene ] ;


    		end
    		fprintf('---OK---\n\n\n');
    		%
    	end
    	fINFOout = [DIRout, filesep, StaticPPIname,'_dynet.info.mat'  ];
    	save(fINFOout,'fDynNetSet','NetStastics', 'GeneSymbol_set')



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end

        tbError


end