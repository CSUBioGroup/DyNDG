function [ mat_distance, matlogic_coChrom, Chrom_ID ,GeneStart, GeneEnd , geneset00 ]  = getGeneChromDistance(Chrom_ID, GeneStart,GeneEnd, Disturbance ) 
   
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
   mat_distance = abs(  (GeneStart-GeneEnd')); 
   mat_distance = min(mat_distance, abs(  (GeneEnd  -GeneStart'))  ); 
%%     
    [uChrom_ID,~,ic] = unique(Chrom_ID); 
    ids = [1:length( uChrom_ID )]'; 
    Chrom_IDnumeric = ids(ic);
    matlogic_coChrom = sparse( Chrom_IDnumeric==Chrom_IDnumeric'  );  
    mat_distance( ~matlogic_coChrom )  = 0 ; 
    logic_isinf_isnan_and_coChrom = sparse( isinf(mat_distance) | isnan( mat_distance)   ) & matlogic_coChrom ; 
    if Disturbance
        max_location = max( GeneStart(~isinf(GeneStart))) ;  
        n1 = nnz(  logic_isinf_isnan_and_coChrom );
        if n1>0
            mat_distance( logic_isinf_isnan_and_coChrom ) = max_location + randi( [100,max_location], n1,1 );
        end
    else
        mat_distance( logic_isinf_isnan_and_coChrom ) = NaN; 
    end
    
    mat_distance( logical( speye( size( mat_distance ) ) )) = 0 ;  
    
end
