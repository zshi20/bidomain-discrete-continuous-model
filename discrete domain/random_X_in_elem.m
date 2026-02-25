function X = random_X_in_elem(node, elem, which_elem)

X_vert = node(elem(which_elem,:)',:);

J = [X_vert(2,:)-X_vert(3,:); X_vert(1,:)-X_vert(3,:)]';

while true
ksi1 = rand; ksi2 = rand; ksi3 = 1-ksi1-ksi2;
if ksi3 > 0
    break;
end
end

X = X_vert(3,:)' + J*[ksi1; ksi2];

end