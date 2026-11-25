
format long g;
%prije ovog treba u command window napisati ovo, onda se brzo učitava
%G=load(dico.gph);
%save('dico.mat', G);
load('dico.mat');


fid_idx=fopen('index.txt', 'r'); %otvara index.txt, 'r' znači samo za čitanje
sve_rijeci=textscan(fid_idx, '%s'){1}; %niz riječi iz rječnika
fclose(fid_idx);

stop_rijeci={'a', 'an', 'the', 'and', 'or', 'but', 'if', 'then', 'else', 'when', 'at', 'by', 'for', 'with', 'about', 'against', 'between', 'into', 'through', 'during', 'before', 'after', 'above', 'below', 'to', 'from', 'in', 'out', 'on', 'off', 'over', 'under', 'is', 'are', 'was', 'were', 'be', 'been', 'being', 'have', 'has', 'had', 'do', 'does', 'did', 'of', 'that', 'this', 'it', 'as', 'which', 'than', 'whose'};
stop_indeksi=find(ismember((sve_rijeci), (stop_rijeci))); %traži indekse stop rijeci
G=G(~ismember(G(:, 1), stop_indeksi)& ~ismember(G(:, 2), stop_indeksi), :);%zadržava u G retke u kojima niti prva niti druga riječ nisu stop riječi

disp("----------------------------\n");
%r = 54962; % kick
%r=47339; %huge
%r = 11897; % boss
%r = 27893; % disappear
%r = 71282; % parallelogram
rijec=input('Unesite zeljenu rijec: ', 's');
[postoji, r]=ismember(rijec, sve_rijeci);
if ~postoji
  disp("Nema trazene rijeci u rjecniku.")
  return
endif

I = G(G(:, 2)==r, 1); %vrhovi koji pokazuju na r
J = G(G(:, 1)==r, 2); %vrhovi na koje r pokazuje

K = unique([r; I(:); J(:)]); %svi vrhovi u podgrafu

bridovi=ismember(G(:,1), K) & ismember(G(:,2), K); %vektor duljine kao G, 0 ako taj par nije u K, 1 ako je

sum(bridovi);
G_r=G(bridovi, :); %matrica kao G, ali samo s bridovima među vrhovima iz K

dim=size(K, 1); %broj vrhova u podmatrici susjedstva
L=zeros(dim, dim);

br_bridova=size(G_r, 1);

for i=1:dim
  for j=1:dim
    par=find(G_r(:, 1)==K(i) & G_r(:, 2)==K(j));
    %par je indeks retka u G_r u kojem se nalazi par (K(i), K(j))
    if ~isempty(par) %ako taj redak u G_r postoji
      L(i, j)=1;
    end
  end
end

[v, d] = eigs(L' * L + L * L', 1); % ima normu 1
v = abs(v);
[v, I] = sort(v, 'descend');

ind = K(I);

% ovaj blok koda mice sve rijeci prije zadane
%while (ind(1) ~= r)
  %ind(1) = [];
%endwhile
max = length(ind);

% ispisuje sinonime po izracunatoj slicnosti
for j = 1 : max
  disp(sve_rijeci{ind(j)});
endfor

disp("----------------------------\n");



