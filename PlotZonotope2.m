function pppppp = PlotZonotope2(G,c,plotnum,VecFlag)
% PLOTZONOTOPE(G,plotnum): given a d x n generator matrix G, plot the
% correponding zonotope as a visual aid. d must equal 2 or 3.
G = 2.*G;
if nargin < 4
    VecFlag = 1;
end

if nargin < 3
    plotnum = 1;
end

if nargin < 2
    G = randn(3,5);
end

[d,n] = size(G);
assert(d < 4);

if d == 2
    % Plot the generators as vectors in the plane
    figure(plotnum)
    for i = 1:n
        vector = [0,0;G(:,i)'];
%         plot(vector(:,1),vector(:,2),'-','color','b','linewidth',2)
        hold on
    end

    % Plot all points of the zonotope

    S = CubetoZonotope(G');    
    
%     S(:,1) =  S(:,1) +c(1)-center(1);
%     S(:,2) =  S(:,2); +c(2) - center(2);
%     plot(S(:,1),S(:,2),'.','markersize',20,'color','b')
    hold on
    axis equal

    % Plot the zonotope using the delaunay trianglulation
    DT = delaunayTriangulation(S);
    K = convexHull(DT);
    Vertices = DT.Points(K,:);
    [nn,~] = size(Vertices);
for i_n = 1:nn
    ver(i_n) = Vertices(i_n,1)^2+Vertices(i_n,2)^2;
end
[~,I_n] = max(ver);
if nn==9
    if I_n == 4
        I_n  = I_n+1;
    end
end

    center(1) = 0.5*Vertices(I_n,1);
    center(2) = 0.5*Vertices(I_n,2);
    Vertices(:,1) =  Vertices(:,1) +c(1)-center(1);
    Vertices(:,2) =  Vertices(:,2) +c(2) - center(2);
   pppppp =  plot(Vertices(:,1),Vertices(:,2),'color','b');
%     title(['A Small 2D Zonotope with ' num2str(n) ' generators'])
    
else
    figure(plotnum)
    
    if VecFlag   % if desired, highlight the generator vectors
        for i = 1:n
            vector = [0,0,0;G(:,i)'];
%             plot3(vector(:,1),vector(:,2),vector(:,3),'-','color','r','linewidth',2)
            hold on
        end
    end
    grid on

    % Plot all points of the zonotope
    figure(plotnum)
    S = CubetoZonotope(G');
%     plot3(S(:,1),S(:,2),S(:,3),'.','markersize',20,'color','b')
    hold on
    axis equal

    % Plot the zonotope using the delaunay trianglulation

    DT = delaunayTriangulation(S);
    K = convexHull(DT);
    trisurf(K,DT.Points(:,1),DT.Points(:,2),DT.Points(:,3),'FaceColor','g','FaceAlpha',.2)
    title(['A Small 3D Zonotope with ' num2str(n) ' generators'])
end
    
end
    