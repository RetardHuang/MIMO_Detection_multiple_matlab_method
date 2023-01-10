function [pos_out] = GTA(TxRx,HTH,HTy,N0,Iter,Constellations)

N = 2*TxRx.Ntx;
L = length(Constellations);
sigma = sqrt(N0);
err=0;
Es = sum(Constellations.^2)/L;
% ////////////////////////
% /////// Tree BP ////////
% ////////////////////////

% Compute weights of total graph 
C=sigma^2*inv(HTH+(sigma^2/Es)*eye(N)); % Covariance matrix 
w=C.^2./(diag(C)*diag(C)');

% As we look for maximum weights, rewrite diagonals of w as -1
% (we are not interested in weights of nodes with theirselves)

w_aux=w-2*diag(diag(w));

% Compute maximum spanning

parent=zeros(N,1); % Keep father of each node, i.e. parent(i) 
% is the parent node of i  

% Assume that root node is node 1
parent(1)=0;
inTree=1; % Keep nodes added to the tree 
Tree=zeros(N); % Edges of tree are keep as '1'

for k=1:N-1
   % We look for maximal weights of nodes inside the tree
   [max_ws,inds]=max(w_aux(inTree,:),[],2);
   % We only save the maximum weigth of previous weights 
   [max_w,ind]=max(max_ws);
   % Write as -1 the edges weigth inside the tree (to not consider
   % them in future iterations) 
   new_node=inds(ind);
   w_aux(inTree,new_node)=-1;
   w_aux(new_node,inTree)=-1;
   % Add new node to the tree 
   parent(new_node)=inTree(ind);
   inTree=[inTree new_node];
   Tree(new_node,parent(new_node))=1;
   Tree(parent(new_node),new_node)=1;
end

% MMSE estimator
%z=inv(H'*H+(sigma^2/Es)*eye(N))*H'*y;
z=sigma^(-2)*C*HTy; % This way avoids to do an inverse again

% f(x1; z, C)
f_x1=@(x1)(exp(-0.5*((x1-z(1)).^2)/C(1,1)));

% f(xi|xj; z, C) (where xj=x_p(i), i.e, the parent node of i)
f_xi_xj=@(xi,zi,xj,zj,Cij,Cii,Cjj)(exp(-0.5*(((xi-zi)-(Cij/Cjj)*(xj-zj)).^2)/(Cii-(Cij^2/Cjj))));

% Look for the tree leaves and its edges 
n_edges=sum(Tree,1)'; % Number of edges of each node 
n_children=[n_edges(1);n_edges(2:N)-1]; % Number of chilgren of each node 
leaf_node=n_children==0; % Leaves are kept as '1' 

% ////////////////////
% /////// BP /////////
% ////////////////////

% //////////////////////////////////////////////
% ///// DOWNWARD BP MESSAGES (m_i-->p(i)) //////
% //////////////////////////////////////////////

messages_down=ones(N,L); % Each column is associated with an alphabet symbol
% Initialized as '1' to merge (25) and (26). This means we
% consider messages from leaves to children as '1'

i_nodes=find(leaf_node==1); % We start from leaves 

n_messages_to_father=zeros(N,1); % Keep the number os messages sended to each parent node

while (isempty(i_nodes)==0)
    i_nodes_sig=[];
    for k=1:length(i_nodes)
        i=i_nodes(k);
        pi=parent(i);
        xi=Constellations;
        children_i=find(parent==i);
        messages_children=prod(messages_down(children_i,:),1); % multiply messages from children to parent node i

        for kk=1:L
            xpi=Constellations(kk);
            messages_down(i,kk)=sum(f_xi_xj(xi,z(i),xpi,z(pi),C(i,pi),C(i,i),C(pi,pi)).*messages_children);
        end

        n_messages_to_father(pi)=n_messages_to_father(pi)+1;

        if (n_edges(pi)-1==n_messages_to_father(pi) && pi~=1) % If parent pi is the root node, it can not added as next node because downward messages end in this node
            i_nodes_sig=[i_nodes_sig; pi]; % pi is added as next node because it has already received all messages from its children and it can send a message to its parent 
        end
    end
    i_nodes=i_nodes_sig;
end


% ////////////////////////////////////////////
% ///// UPWARD BP MESSAGES (m_p(i)-->i) //////
% ////////////////////////////////////////////

messages_up=ones(N,L); % Each column is associated with an alphabet symbol
% Initialized as '1' to merge (27) and (28). This means we
% consider messages from parent of root node to root node as '1'

i_nodes=find(parent==1); % We start from root node to its children 

while(isempty(i_nodes)==0)
    i_nodes_sig=[];
    for k=1:length(i_nodes)
        i=i_nodes(k);
        pi=parent(i);
        xpi=Constellations;

        children_pi=find(parent==pi);
        ind_children=find(children_pi~=i);
        children_pi=children_pi(ind_children);% We remove child of pi, i.e, i, because we are doing it now 
        messages_children=prod(messages_down(children_pi,:),1); % multiply messages from children of node p(i) to p(i) (less i) 

        ppi=parent(pi); % Parent of parent node of i p(p(i))
        message_parent=messages_up(pi,:);

        for kk=1:L
            xi=Constellations(kk);
            messages_up(i,kk)=sum(f_xi_xj(xi,z(i),xpi,z(pi),C(i,pi),C(i,i),C(pi,pi)).*messages_children.*message_parent);
        end

        children_i=find(parent==i);
        i_nodes_sig=[i_nodes_sig; children_i];
    end
    i_nodes=i_nodes_sig;
end



% ////////////////////////////
% ///// BELIEF VARIABLE //////
% ////////////////////////////

% After message passing, we compute the belief at each variable

belief_x=zeros(N,L);
for i=1:N
    xi=Constellations;
    children_i=find(parent==i);
    if i==1
        belief_x(i,:)=f_x1(xi).*prod(messages_down(children_i,:),1);
    else
        belief_x(i,:)=messages_up(i,:).*prod(messages_down(children_i,:),1);
    end
end


% ////////////////////
% ///// DECODER //////
% ////////////////////

% Decoding

[max_belief,ind]=max(belief_x,[],2);
x_decod=Constellations(ind)';

mean = reshape(x_decod,N/2,[]);
symb = mean(:,1) + 1j*mean(:,2);
[~,pos] = min(abs(symb*ones(1,2^L) - ones(N/2,1)*TxRx.Constellations)');
pos_out = pos'*ones(1,Iter);
end