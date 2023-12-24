function [ mat_distance, matlogic_coChrom, Chrom_ID ,GeneStart, GeneEnd , geneset00 ]  = getGeneChromDistance(Chrom_ID, GeneStart,GeneEnd, Disturbance ) 
% Ϊ��ȷlocation�Ļ��� λ�ü�Ϊ inf 
% % Chrom_ID = 'update2020_Menche_Science2015PPI_TableGenes_genelocation&gene1000neigbors.mat'
% Chrom_ID ������TableGenes_genelocation ���������ļ���     
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
    % ���� ���л���� ���룬������Ⱦɫ��  
   mat_distance = abs(  (GeneStart-GeneEnd')); 
   mat_distance = min(mat_distance, abs(  (GeneEnd  -GeneStart'))  ); 
%%     mat_distance = abs(  (GeneMean-GeneMean')); 
    % % mat_distance( isnan(mat_distance)) = inf; 
    % ��ǹ�ͬȾɫ�� �Ļ����  
    % ת�� ��ֵ id ���ټ���  
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
        % �Ŷ���������ʱ���� ����ͬȾɫ��Ļ��� ,������Զ������ֲ�������������� 
        max_location = max( GeneStart(~isinf(GeneStart))) ; 
% %         mat_distance( ~matlogic_coChrom )  =  2*max_location+randi( [ 100, max_location], nnz(~matlogic_coChrom),1  ) ; 
        n1 = nnz(  logic_isinf_isnan_and_coChrom );
        if n1>0
            %����ͬȾɫ��ģ���ȱʧλ����Ϣ�Ļ���ԣ������޶��Ŷ� 
            mat_distance( logic_isinf_isnan_and_coChrom ) = max_location + randi( [100,max_location], n1,1 );
        end
    else
% % %         mat_distance( ~matlogic_coChrom ) =  0  ;   % ����ͳһȾɫ��� ����ԣ����� 0  ���ھ���ϡ�軯 
        mat_distance( logic_isinf_isnan_and_coChrom ) = NaN;  % ��ͬһȾɫ�壬��ȱʧλ����Ϣ�Ļ��򣨶ԣ�������Ϊ��Чֵ nan 
    end
    % % mat_distance( ~matlogic_coChrom ) =  -1 ; 
    % ��ͬ����֮�� ���� Ϊ 0  
    mat_distance( logical( speye( size( mat_distance ) ) )) = 0 ;  
    %%%mat_distance = sparse(double( mat_distance ) ) ;
    
end
