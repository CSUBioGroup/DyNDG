function Matrix_Disease_Gene=getMap2RefMatrix( M_dis_gene, DisIDset,GeneSymbolSet, refDisIDset,refGeneSymbolSet  )
%  Get the mapping to the reference matrix 
% map a matrix to a reference set
% so that it has the same size of matrix, the same annotated terms and orders for rows
% and columns
% % M_dis_gene, DisIDset,GeneSymbolSet, refDisIDset,refGeneSymbolSet 
% % map2refmatrix
% % refDisIDset = data_all.DataDiseaseNet_OMIM6449DiseaseGeneSet.DiseaseGeneTable{:, 'Phenotype_MIM_Number_unique'};
% % refGeneSymbolSet = data_all.DataGeneNet_1Net.TableGenes{:, 'symbol'};
% % DisIDset     = TableOMIM{:,2}; 
% % GeneSymbolSet= TableSymbol{:,2};  
% % 2020-9-19
% % Xiang 
if ~isempty(DisIDset) && isempty(GeneSymbolSet)
% %     map2refmatrix( M_dis_gene, DisIDset,[], refDisIDset,[]  )
    [~,ia_dis,ib_dis]   = intersect(refDisIDset,DisIDset);  
    n_dis  = length( refDisIDset );
    n_gene = size(M_dis_gene,2);
    Matrix_Disease_Gene = zeros(n_dis, n_gene ); 
    Matrix_Disease_Gene(ia_dis ,: )= M_dis_gene(ib_dis, :);  
    
elseif isempty(DisIDset) && ~isempty(GeneSymbolSet)
% %     map2refmatrix( M_dis_gene, [],GeneSymbolSet, [],refGeneSymbolSet  )
    [~,ia_gene,ib_gene] = intersect(refGeneSymbolSet,GeneSymbolSet); 
    n_dis  = size(M_dis_gene,1);
    n_gene = length( refGeneSymbolSet );
    Matrix_Disease_Gene = zeros(n_dis, n_gene ); 
    Matrix_Disease_Gene(: ,ia_gene )= M_dis_gene(:, ib_gene);  
    
elseif ~isempty(DisIDset) && ~isempty(GeneSymbolSet)
% % map2refmatrix( M_dis_gene, DisIDset,GeneSymbolSet, refDisIDset,refGeneSymbolSet  )
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
