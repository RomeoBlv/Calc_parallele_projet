clear all
close all
nb_de = 9; %--> j'ai fixé le niveau de décomposition arbitrairement ici.
for k = 1:nb_de
  name = sprintf('final_u_%d.dat',k-1);
  fid = fopen(name,'r','l');
  mydata=fread(fid,'double');  
   for(i=1:510) %-->nombre de lignes de l'image sans tenir compte de la première et de la dernière ligne donc
        for(j=1:510) %-->nombre de colonnes de l'image sans tenir compte de la première et de la dernière colonne donc

             A(i,j,k)=mydata((i-1)*510+j);
        end
    end 
end

for(k=1:nb_de)
    figure;imshow(sum(A(:,:,1:k),3),[]);
end