function [idx_neg_test_ArtiLinkInterval] = getGeneNeigbors_ArtitificalLinkageInterval(mat_genedistance_copy,matlogic_coChrom_copy, idx_pos_test, idx_mask_genes, n_neighbors_artificial_linkage_interval          )
% %  idx_mask_genes  include  ind_pos_test and ind_pos_train so as to
% select negative neighbors only 

%idx_mask_genes = idx_pos  
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


% % % % 
% % % % mat_genedis =  mat_genedistance_copy(ind_pos_test,:); 
% % % % mat_genedis(:, idx ) = inf ; 
% % % % [~,testgeneneigbors] = sort(mat_genedis,2,  'ascend') ; 
% % % % testgeneneigbors = testgeneneigbors(:,1:n_neighbors_artificial_linkage_interval);  
% % % % ind_neg_test_ArtiLinkInterval =  unique( testgeneneigbors ); 
