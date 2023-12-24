function Matrix_Disease_Gene=getMap2RefMatrix( M_dis_gene, DisIDset,GeneSymbolSet, refDisIDset,refGeneSymbolSet  )
if ~isempty(DisIDset) && isempty(GeneSymbolSet)
    [~,ia_dis,ib_dis]   = intersect(refDisIDset,DisIDset);  
    n_dis  = length( refDisIDset );
    n_gene = size(M_dis_gene,2);
    Matrix_Disease_Gene = zeros(n_dis, n_gene ); 
    Matrix_Disease_Gene(ia_dis ,: )= M_dis_gene(ib_dis, :);  
    
elseif isempty(DisIDset) && ~isempty(GeneSymbolSet)
    [~,ia_gene,ib_gene] = intersect(refGeneSymbolSet,GeneSymbolSet); 
    n_dis  = size(M_dis_gene,1);
    n_gene = length( refGeneSymbolSet );
    Matrix_Disease_Gene = zeros(n_dis, n_gene ); 
    Matrix_Disease_Gene(: ,ia_gene )= M_dis_gene(:, ib_gene);  
    
elseif ~isempty(DisIDset) && ~isempty(GeneSymbolSet)
    [~,ia_dis,ib_dis]   = intersect(refDisIDset,DisIDset); 
    [~,ia_gene,ib_gene] = intersect(refGeneSymbolSet,GeneSymbolSet); 
    n_dis  = length( refDisIDset );
    n_gene = length( refGeneSymbolSet );
    Matrix_Disease_Gene = zeros(n_dis, n_gene ); 
    Matrix_Disease_Gene(ia_dis ,ia_gene )= M_dis_gene(ib_dis, ib_gene);  
    
else
    error( 'errorerrorerrorerrorerrorerror'  )
    
    
end
 

end
