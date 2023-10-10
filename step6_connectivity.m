function [connectivity, unconnected_pores]=step6_connectivity(neighbour, throat_body_new)
% This function computes the connectivity of the microstructure and outputs the label of unconnected or isolated pore bodies
% ------
% Author: Nimisha Roy
% e-mail: nroy9@gatech.edu
% ------
%------------------------------------------------------------------------------START CODE------------------------------------------------------------------------------
nodes=reshape(neighbour(:,1:2),[2*size(neighbour,1),1]);%existing pore bodies
x=unique(nodes,'rows');%pore body index
freq=tabulate(nodes);
freq=freq(:,2);
unconnected_voxels=0;unconnected_pores=[];
for i=1:size(freq,1)
    if freq(i,1)== 1% unconnected voxels due to dead end pores
        unconnected_voxels=unconnected_voxels+size(throat_body_new(throat_body_new==x(i)),1);% unconnected due to dead end pore
        unconnected_pores=[unconnected_pores;x(i)];
    end
    %disp([ num2str(i),' processed out of ',num2str(size(freq,1))])
end

% to find unconnected pores due to isolated pores...i.e. ones that are not in neighbour but present in throat_body
sequence=linspace(1,x(end,1),x(end,1)).';%all pore indices
k=setdiff(sequence,x);%pore indices not in neighbour
if size(k,1)>0 && sum( arrayfun (@(a) (size(throat_body_new(throat_body_new==a),1)>0),k))== size(k,1)%pore indices not inneighbour but there in throat_body are isolated pores
    unconnected_voxels=unconnected_voxels+sum( arrayfun (@(a) (size(throat_body_new(throat_body_new==a),1)),k));% unconnectenumberOfParticlesd due to dead end pore
    unconnected_pores=[unconnected_pores;k];
end
connectivity=(size(throat_body_new,1)-unconnected_voxels)*100/size(throat_body_new,1);

% G = graph(neighbour(:,1),neighbour(:,2),neighbour(:,3)) ;
% figure, plot(G)
%----------------------------------------------------------------- END CODE--------------------------------------------------------------------------
