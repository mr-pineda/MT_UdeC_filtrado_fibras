b="13959";

b_2d_file = fopen("2d_"+b+".bundles.bin","r");
b_3d_file = fopen(b+".bundlesdata","r");

B_2d = {};
B_3d = {};
fiber=1;
while (~feof(b_2d_file))
    A_2d = fread(b_2d_file,1,'int');    
    cord_2d = zeros(A_2d,2);
    
    for i=1:A_2d
        for j=1:2
            cord_2d(i,j) = fread(b_2d_file,1,'float');
        end
    end

    B_2d{fiber} = cord_2d;
    fiber=fiber+1;
    clear cord_2d
end

fiber=1;
while (~feof(b_3d_file))
    A_3d = fread(b_3d_file,1,'int');
    cord_3d = zeros(A_2d,3);

    for i=1:A_3d
        for j=1:3  
            cord_3d(i,j) = fread(b_3d_file,1,'float');
        end   
    end
    B_3d{fiber} = cord_3d;

    fiber=fiber+1;
    clear cord_3d
end

fclose(b_2d_file);
fclose(b_3d_file);

clear b_2d_file b_3d_file A_2d A_3d i j ans fiber

%% Ploteo de 1 fibra

f = 1;

close("all")

figure("Name","vista 2D")
plot(B_2d{f}(:,1),B_2d{f}(:,2));

figure("Name","vista 3D")
plot3(B_3d{f}(:,1),B_3d{f}(:,2),B_3d{f}(:,3));



%% PRUEBAS DE ISOMAP Y MDS
clc
fibra = 1;
n_fib = length(B_3d)-1;

k = 6;
B_mds = cell(1,n_fib);

for f = 1:n_fib
    n_pt = length(B_3d{f});
    m_dist = zeros(n_pt);
%     isomap = zeros(n_pt);   
    
    for i=1:n_pt
        for j=1:n_pt
            m_dist(i,j) = norm(B_3d{f}(i,:)-B_3d{f}(j,:));
        end
    end
    B_mds{f} = cmdscale(m_dist,2);
end
% graph_mat = grmath(m_dist,k);

% for i=1:N
%     isomap(i,:) = dijkstra(graph_mat,i);
% end

% newcord_iso = cmdscale(isomap,2);
newcord_mds = cmdscale(m_dist,2);

%% Ploteo del fascículo

close("all")

figure("Name","bundle 3D")
hold all
for i=1:n_fib
    plot3(B_3d{i}(:,1),B_3d{i}(:,2),B_3d{i}(:,3));
end
hold off

figure("Name","bundle ISOMAP")
hold all
for i=1:n_fib
    plot(B_2d{i}(:,1), B_2d{i}(:,2));
end
hold off

figure("Name","bundle MDS")
hold all
for i=1:n_fib
    plot(B_mds{i}(:,1), B_mds{i}(:,2));
end
hold off

%% FUNCIONES
function graph = grmath(dist,k)
    k=k+1;
    g1 = dist;
    for i=1:length(dist)
        for j=1:length(dist)-k
            [~,idx] = max(g1(i,:));
            g1(i,idx)=-1;
        end
    end
    graph = inf(length(dist),length(dist));
    l1 = g1 ~= -1;
    l2 = g1' ~= -1;
    l3 = or(l1,l2);
    for i=1:length(dist)
        for j=1:length(dist)
            if l3(i,j)
                graph(i,j) = dist(i,j);
            end
        end
    end
end


function distances = dijkstra(map,start)
    N =length(map);
    distances(1:N)=inf;
    visited(1:N) = 0;

    distances(start)=0;
    while sum(visited) < N
        candidates(1:N) = inf;
        for i = 1:N
            if visited(i) == 0
                candidates(i)=distances(i);
            end
        end
        
        [m_val,m_idx]=min(candidates);
        for i = 1:length(map)
            new_dist = m_val + map(m_idx,i);
            if new_dist <distances(i)
                distances(i) = new_dist;
            end
        end
        visited(m_idx)=1;
    end
end

function B_3d = readBundle(name)
b_3d_file = fopen(name+".bundlesdata","r");
B_3d = {};
fiber=1;
    while (~feof(b_3d_file))
        A_3d = fread(b_3d_file,1,'int');
        cord_3d = zeros(A_2d,3);
        
        for i=1:A_3d
            for j=1:3  
                cord_3d(i,j) = fread(b_3d_file,1,'float');
            end   
        end
        B_3d{fiber} = cord_3d;
        
        fiber=fiber+1;
    end
fclose(b_3d_file);
end