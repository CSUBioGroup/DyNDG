function [idx_neg_test_ArtiLinkInterval] = getGeneNeigbors_ArtitificalLinkageInterval(mat_genedistance_copy,matlogic_coChrom_copy, idx_pos_test, idx_mask_genes, n_neighbors_artificial_linkage_interval          )

% select negative neighbors only 

    mat_genedistance_copy =  mat_genedistance_copy(idx_pos_test,:); 
    matlogic_coChrom_copy = matlogic_coChrom_copy(idx_pos_test,:); 
    %
    n_nonCoChrom = nnz(~matlogic_coChrom_copy); 
    Dmax         = max(  mat_genedistance_copy(:)  ) ; 
    mat_genedistance_copy( ~matlogic_coChrom_copy ) = Dmax +randi( [ 100, Dmax],  n_nonCoChrom ,1  ) ; 
    %
    mat_genedistance_copy(:, [ idx_mask_genes'] ) = inf ; 
    [~,testgeneneigbors] = sort(mat_genedistance_copy,2,  'ascend') ; 
    testgeneneigbors = testgeneneigbors(:,1:n_neighbors_artificial_linkage_interval);  
    idx_neg_test_ArtiLinkInterval =  unique( testgeneneigbors ); 
  
end

