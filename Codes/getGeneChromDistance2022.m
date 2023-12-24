function [ mat_distance, matlogic_coChrom, Chrom_ID ,GeneStart, GeneEnd , geneset00 ]  = getGeneChromDistance(Chrom_ID, GeneStart,GeneEnd, Disturbance ) 
% 为明确location的基因， 位置记为 inf 
% % Chrom_ID = 'update2020_Menche_Science2015PPI_TableGenes_genelocation&gene1000neigbors.mat'
% Chrom_ID 可输入TableGenes_genelocation 或者数据文件名     
    if isempty(Chrom_ID)
        Chrom_ID = 'update2022-3-24_TableGeneLocation.mat'  ; 
    end
    if isa(Chrom_ID,'table')
        TableGeneLocation = Chrom_ID ; 
        Chrom_ID  = TableGeneLocation.Chrom_ID ;   
        GeneStart = TableGeneLocation.GeneStart ;     
        GeneEnd   = TableGeneLocation.GeneStart ; 
        GeneMean  = TableGeneLocation.GeneMean ; 
		geneset00 = TableGeneLocation.symbol; 
		
    elseif ischar(Chrom_ID) 
        data = load(Chrom_ID); TableGeneLocation = data.TableGeneLocation; data =[]; 
        Chrom_ID  = TableGeneLocation.Chrom_ID ;   
        GeneStart = TableGeneLocation.GeneStart ;     
        GeneEnd   = TableGeneLocation.GeneStart ; 
		GeneMean  = TableGeneLocation.GeneMean ; 
        geneset00 = TableGeneLocation.symbol; 
    end
    % % % % % % % % % % 
% %     n_genes = length( Chrom_ID ) ; 
    % % mat_distance_base = randi([5,10],[5 6])
    % 计算 所有基因对 距离，不区分染色体  
   mat_distance = abs(  (GeneStart-GeneEnd')); 
   mat_distance = min(mat_distance, abs(  (GeneEnd  -GeneStart'))  ); 
%%     mat_distance = abs(  (GeneMean-GeneMean')); 
    % % mat_distance( isnan(mat_distance)) = inf; 
    % 标记共同染色体 的基因对  
    % 转成 数值 id 加速计算  
    [uChrom_ID,~,ic] = unique(Chrom_ID); 
    ids = [1:length( uChrom_ID )]'; 
    Chrom_IDnumeric = ids(ic);
    matlogic_coChrom = sparse( Chrom_IDnumeric==Chrom_IDnumeric'  );  
% %     tic
% %     n_genes = length( Chrom_ID ) ; 
% %     matlogic_coChrom = sparse(  strcmp( repmat(Chrom_ID,1,n_genes) , repmat(Chrom_ID', n_genes , 1) )  );   
    mat_distance( ~matlogic_coChrom )  = 0 ; 
    logic_isinf_isnan_and_coChrom = sparse( isinf(mat_distance) | isnan( mat_distance)   ) & matlogic_coChrom ; 
    if Disturbance
        % 扰动便于排序时囊括 不在同染色体的基因 ,其距离更远，随机分布，便于随机抽样 
        max_location = max( GeneStart(~isinf(GeneStart))) ; 
% %         mat_distance( ~matlogic_coChrom )  =  2*max_location+randi( [ 100, max_location], nnz(~matlogic_coChrom),1  ) ; 
        n1 = nnz(  logic_isinf_isnan_and_coChrom );
        if n1>0
            %在相同染色体的，但缺失位置信息的基因对，距离限定扰动 
            mat_distance( logic_isinf_isnan_and_coChrom ) = max_location + randi( [100,max_location], n1,1 );
        end
    else
% % %         mat_distance( ~matlogic_coChrom ) =  0  ;   % 不在统一染色体的 基因对，距离 0  便于矩阵稀疏化 
        mat_distance( logic_isinf_isnan_and_coChrom ) = NaN;  % 在同一染色体，但缺失位置信息的基因（对），距离为无效值 nan 
    end
    % % mat_distance( ~matlogic_coChrom ) =  -1 ; 
    % 相同基因之间 距离 为 0  
    mat_distance( logical( speye( size( mat_distance ) ) )) = 0 ;  
    %%%mat_distance = sparse(double( mat_distance ) ) ;
    
end
