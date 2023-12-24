Disease=[{'AML','CML','CLL'}];
DIR_StaticPPI       = 'StaticPPI' ;
Disease_DyNetSet=struct('AML',{'GSE122917_ksigma2'},'CML',{'GSE47927_ksigma2'},'CLL',{'GSE2403_ksigma2'});
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

for i_Dis=1:length(Disease)
    DisName=Disease{i_Dis};%要处理的疾病
    Root=['Dis_',DisName];
    DIRdisgene          = sprintf('%s_DisGeneSet', DisName ), MarkStr_DisGeneSetFile  = sprintf('%s_DisGeneSet' , DisName )
    DIR_FormatMatData   = sprintf('%s_FormatMatData', DisName) , if ~exist([Root,filesep,DIR_FormatMatData],'dir'); mkdir([Root,filesep,DIR_FormatMatData]); end
    DEG_Stage=sprintf('%s_DEG_', DisName );
    DyExpDataTypeSet = cellstr(Disease_DyNetSet(1).(DisName)); %选择疾病对应想要处理的GSE,目前只能选一个
    for i_DyExpDataType = 1:length( DyExpDataTypeSet )
        DyNetType = DyExpDataTypeSet{i_DyExpDataType};
        switch DyNetType

            case 'GSE122917_ksigma2'
                IsNet = true, MarkStr_gExpOrNetDataFile = 'net_Stage',logstr = 'log',DyNetMethod = 'ksigma2'
                DIR_gExpOrNetDataSet = 'AML_GSE122917_gExp3Stage_log_ksigma2',
            case 'GSE47927_ksigma2'
                IsNet = true, MarkStr_gExpOrNetDataFile = 'net_Stage',logstr = 'log',DyNetMethod = 'ksigma2'
                DIR_gExpOrNetDataSet = 'CML_GSE47927_gExp4Stage_log_ksigma2',
            case 'GSE2403_ksigma2'
                IsNet = true, MarkStr_gExpOrNetDataFile = 'net_Stage',logstr = 'log',DyNetMethod = 'ksigma2'
                DIR_gExpOrNetDataSet = 'CLL_GSE2403_gExp3Stage_log_ksigma2',


    		otherwise
    			error(['!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'])
        	end
        	%%%%%%
        	disp('')
        	if ~exist([Root,filesep,DIR_gExpOrNetDataSet],'dir'); error( ['There is no directory: ' ,DIR_gExpOrNetDataSet] ); end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            DataSet=[];
            DisGeneSet   = table ;
            DataSet.DisGeneSet = table;
            DEG_Stage_Set=table;
            DataSet.DEG_Stage_Set=table;

            % % % % % % % % % % % % % % % %
            SeedGeFileSet   = dir( [Root,filesep,DIRdisgene,filesep,'*',MarkStr_DisGeneSetFile,'*'])
            n_SeedGeFile    = size(SeedGeFileSet,1)
            if n_SeedGeFile<1; error('There is no effective Seed-gene set.'); end
            DEG_set=dir(['*',Root,filesep,DIRdisgene,filesep,DEG_Stage,'*']);
            DisGene_set=setdiff({SeedGeFileSet.name},{DEG_set.name})
            n_DEGFile=size(DEG_set,1);
            if n_DEGFile<1; error('There is no effective DEG set.'); end
            stage=ceil(sqrt(n_DEGFile));
            for i_stage=1:stage
                n_DEG_set=dir([Root,filesep,DIRdisgene,filesep,DEG_Stage,num2str(i_stage),'vs','*.txt']);
                n_DEG_File    = size(n_DEG_set,1);
                dg=table;
                for i_DEG_set =1:n_DEG_File
                    fname    = n_DEG_set(i_DEG_set).name;
                    ind      = strfind(fname,'.');
                    rowname  = fname(1:ind(end)-1);
                    fDisGene = [Root,filesep,DIRdisgene,filesep,fname ];
                    dg_new   = readtable( fDisGene ,'ReadVariableNames',false );
                    dg_new=dg_new(:,[1 4 5]);
                    dg_new(:,2)=array2table(abs(table2array(dg_new(:,2))));%T统计量取绝对值
                    dg_new.Properties.VariableNames = {'Gene' ['T_',num2str(i_DEG_set)] ['Pvalue_',num2str(i_DEG_set)]};
                    if(isempty(dg))
                        dg       =dg_new;
                    else
                        dg = outerjoin(dg,dg_new,'MergeKeys',true);
                    end
                end
                k=dg.Properties.VariableNames;
                dg.Mean = mean(dg{:,k(2:2:n_DEG_File*2)},2);
                dg.Mean=mapminmax(dg.Mean',0,1)';
                DEG_Stage_Set{ {[DisName,'_Stage_',num2str(i_stage)]} ,{[DisName,'_DisGeneSet']} } =  {dg};
            end
            n_DisGeFile=length(DisGene_set)
            if n_DisGeFile<1; error('There is no effective disease-gene set.'); end
            for i_DisGeFile = 1:n_DisGeFile
                fname    = cell2mat(DisGene_set(i_DisGeFile));
                ind      = strfind(fname,'.');
                rowname  = fname(1:ind(end)-1);
                fDisGene = [Root,filesep,DIRdisgene,filesep,fname ];
                dg       = readtable( fDisGene ,'ReadVariableNames',false ); dg = dg{:,1}';
                DisGeneSet{ {rowname} ,{[DisName,'_DisGeneSet']} } =  {dg} ;
            end
            DataSet.DisGeneSet = DisGeneSet;
            DataSet.DEG_Stage_Set=DEG_Stage_Set;
            fDynNetSet=[];
            % % % % % % % % % % % % % % % % % % % % % % % % % %

            if contains(DyNetMethod, 'sigma' )  || contains(DyNetMethod, 'sigma2' )
                % generate co-exp networks
                NamesPPIset = tb_fStaticPPIset.Properties.RowNames;
                for i_PPI = 1: length( NamesPPIset ) 	% % % % % % % % % % % % % % % % % % % % % % % % % %
                    Name_StaticPPI = NamesPPIset{i_PPI}
                    fStaticPPI   = tb_fStaticPPIset{i_PPI,1}{1};
                    threshold    = tb_fStaticPPIset{i_PPI,2}{1};
                    MarkStr      = [Name_StaticPPI,'_',MarkStr_gExpOrNetDataFile]
                    fgExpOrNetDataSet  = dir( [Root,filesep,DIR_gExpOrNetDataSet,filesep, [ MarkStr,'*'] ]  )
                    if length( fgExpOrNetDataSet ) ==0
                        warning( ['There is no network data for ', Name_StaticPPI, ' in ', [Root,filesep,DIR_gExpOrNetDataSet,filesep, [ MarkStr,'*'] ] ] );
                        tbError{size(tbError,1)+1,:}={DyNetType,Name_StaticPPI, [Root,filesep,DIR_gExpOrNetDataSet,filesep, [ MarkStr,'*'] ],fStaticPPI, exist(fStaticPPI,'file')}
                        continue
                    end
                    n_stage            = size(fgExpOrNetDataSet,1)
                    fDynNetSet         = []
                    for i_stage = 1:n_stage
                        fgExpOrNetDataSet   = dir( [Root,filesep,DIR_gExpOrNetDataSet,filesep, [MarkStr,num2str(i_stage),'*'] ]  );
                        if size(fgExpOrNetDataSet,1)>1; disp(fgExpOrNetDataSet); error('There are more than 1 network at this stage.');
                        else
                            fDynNetSet{i_stage,1} = [Root,filesep,DIR_gExpOrNetDataSet,filesep,fgExpOrNetDataSet.name];
                        end
                    end
                    %
                    [DataSet.([DIR_gExpOrNetDataSet,'_', Name_StaticPPI ]) ]   = getNetMatrixSetFromFiles_MultipleStageCoExpOrDyNet(fDynNetSet , fStaticPPI,  DisGeneSet,DEG_Stage_Set, threshold);
                end
            end


            %    HCC_GSE6764_exp_8Stage
            fdatestr   = datestr(now,'yyyy.mmm.dd-HH.MM.SS') ;
            fdatestr   = '';
            fout = [Root,filesep,DIR_FormatMatData,filesep,'OK-DataSet-', DIR_gExpOrNetDataSet,'-NetMatrixSet_StaticNet&DyNet_', fdatestr,'.mat']
            save(fout, 'DataSet')


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            tbError

    end
end