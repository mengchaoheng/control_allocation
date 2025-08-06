function ratio = vview_case(B, plim, P, v_1, v_2,v_all, v_21,v_all1,v_22,v_all2)
  
% VVIEW - View the attainable virtual control set.
% 
%  1) vview(B,plim)
%
% Shows the attainable virtual control set considering actuator
% position constraints, given by { v : v = B*u, umin < u < umax }.
%
%  2) ratio = vview(B,plim,P)
% 
% Compares the set of feasible virtual control inputs when
%
%  a) the actuator redundancy is fully utilized (as above) [blue]
%  b) a linear allocation control law u = Pv is used (BP = I) [red]
%
% The second set is given by { v : umin < P*v < umax }.
%
%  Inputs:
%  -------
% B      control effectiveness matrix (k x m)
% plim   control position limits [min max] (m x 2)
% P      virtual control law matrix (m x k)
% 
%  Outputs:
%  --------
% ratio  The ratio between the sizes (areas, volumes, ...)
%        of the two sets
%
% The result is only graphically illustrated for k = 1, 2, or 3.
%
% See also: VVIEW_DEMO
 
% Model dimensions

if nargin < 7
    v_21= [];v_22= [];v_all1=[];v_all2=[];
end

if nargin < 8
    v_22= [];v_all2=[];
end
if nargin < 4
    v_1=[];v_2=[];v_all=[];v_21=[];v_all1=[];v_22=[];v_all2=[];
end

% ===== 配色参数 =====
AMS_fill_color   = [0.95 0.95 1];   
AMS_edge_color   = [0.10 0.35 0.75];   

INV_fill_color   = [1 1 0.9];   
INV_edge_color   = [0.85 0.33 0.10];   

v_1_color  = [0.00, 0.00, 0.00];      
v_2_color  = [0.00, 0.60, 0.60];     
v_21_color = [0.40, 0.00, 0.70];      
v_22_color = [0.70, 0.50, 0.00];       


  [k,m] = size(B);
  
  % ------------------------------------------------
  %  a) Find maximum attainable virtual control set 
  %     considering constraints.
  % ------------------------------------------------
  
  % Generate matrix to index corners of feasible control set.
  idx = zeros(2^m,m);
  M = 1:m;
  for i = 1:2^m;
    cbin = dec2bin(i-1,m); % '001'
    c = str2num(cbin')'; % [0 0 1]
    c = c(end:-1:1); % [1 0 0]
    idx(i,:) = 2*M - c;
  end

  % Generate corner points of the feasible control set.
  plimT = plim';
  U = plimT(idx)';

  % Compute the corresponding points in the virtual control space
  V = B*U;

  if nargin > 2
    
    % ---------------------------------------------
    %  b) Find attainable virtual control set when
    %     a linear control law u=Pv is used.
    % ---------------------------------------------
    
    % We want to determine where the k-dim. hyperplane Pv
    % intersects the m-dim. hyperbox of feasible controls.
    % To get the corner points of this set, solve
    % Pv = x where x has k specified entries.
    % 
    % Example: m=3, k=1 -> points will lie on surfaces
    %          m=3, k=2 -> points will lie on edges
    
    % Generate index matrix for all combinations of min and max indeces
    % in k dimensions.
    sub_idx = idx(1:2^k,1:k);
    
    Ulin = [];
    % Loop over all combinations of dimensions
    i_dim = nchoosek(1:m,k);
    for i = 1:size(i_dim,1)
      % For each combination, compute the intersections with all
      % possible min/max combinations.
      
      % k-dimensional min/max combinations
      sub_plimT = plimT(:,i_dim(i,:));
      sub_u_boundary = sub_plimT(sub_idx)';
      
      % Determine which virtual control sub_u_boundary corresponds to
      sub_P = P(i_dim(i,:),:);
      if rank(sub_P) == k % Avoid "parallel" cases
        % Solve sub_u_boundary = sub_P v for v
	v = sub_P\sub_u_boundary;
	% Determine the full countol vector (contains sub_u_boundary)
	u_boundary = P*v;
	
	% Store feasible points
	i_feas = feasible(u_boundary,plim);
	Ulin = [Ulin u_boundary(:,i_feas)];
      end
    end

    % Compute the corresponing points in the virtual control space
    Vlin = B*Ulin;
  
  end
    
  % Compute and visualize the convex hull of the set(s)
  clf
  switch k
   case 1
    K = [min(V) max(V)];
    if nargin > 2
      Klin = [min(Vlin) max(Vlin)];
      ratio = diff(Klin)/diff(K);
      
      % Illustrate
      plot(K,[0 0],'b-o',Klin,-[0 0],'r-o')
    else
      plot(K,[0 0],'b-o')
    end
    xlabel('v')
    
   case 2
    hold on; axis equal;

    % 初始化 legend 对象数组和标签数组
    legend_handles = [];
    legend_labels  = [];
    
    % 填充区域
    [K,area1] = convhull(V(1,:),V(2,:));
    h_fill1 = fill(V(1,K), V(2,K), AMS_fill_color,'FaceAlpha', 0.7, ...% FaceColor for AMS 
             'EdgeColor', AMS_edge_color, 'LineWidth', 1.2); % EdgeColor for AMS 
    legend_handles(end+1) = h_fill1;
    legend_labels{end+1}  = 'AS';
    
    if nargin > 2 && exist('Vlin', 'var') && ~isempty(Vlin)
        [Klin,area2] = convhull(Vlin(1,:), Vlin(2,:));
        % h_fill2 = fill(Vlin(1,Klin), Vlin(2,Klin), INV_fill_color, ... % FaceColor for inv
        %          'EdgeColor', INV_edge_color, 'LineWidth', 1.2); % EdgeColor for inv
        % legend_handles(end+1) = h_fill2;
        % legend_labels{end+1}  = '$\Pi_{\mathrm{inv}}$';
        ratio = area2 / area1;
    end
    
    
    % ==== 向量绘制 ====
    if exist('v_1','var') && ~isempty(v_1)
        h = quiver(0, 0, v_1(1), v_1(2), 0, '--','Color', v_1_color, ...
            'LineWidth', 1, 'MaxHeadSize', 0.3, 'AutoScale', 'off');
        legend_handles(end+1) = h;
        legend_labels{end+1}  = '$\nu_h$';
    end
    if exist('v_2','var') && ~isempty(v_2)
        h = quiver(v_1(1), v_1(2), v_2(1), v_2(2), 0, '-.', 'Color', v_2_color,...
            'LineWidth', 1, 'MaxHeadSize', 0.3, 'AutoScale', 'off');
        legend_handles(end+1) = h;
        legend_labels{end+1}  = '$\nu_c$';
    end
    if exist('v_all','var') && ~isempty(v_all)
        h = quiver(0, 0, v_all(1), v_all(2), 0, '-','Color', v_2_color, ...
            'LineWidth', 1, 'MaxHeadSize', 0.3, 'AutoScale', 'off');
        legend_handles(end+1) = h;
        legend_labels{end+1}  = '$\nu$';
    end
    
    if exist('v_21','var') && ~isempty(v_21)
        h = quiver(v_1(1), v_1(2), v_21(1), v_21(2), 0, '-.', 'Color', v_21_color,...
            'LineWidth', 1, 'MaxHeadSize', 0.3, 'AutoScale', 'off');
        legend_handles(end+1) = h;
        legend_labels{end+1}  = '$\nu_c^{\,\prime}$';
    end
    if exist('v_all1','var') && ~isempty(v_all1)
        h = quiver(0, 0, v_all1(1), v_all1(2), 0, '-','Color', v_21_color, ...
            'LineWidth', 1, 'MaxHeadSize', 0.3, 'AutoScale', 'off');
        legend_handles(end+1) = h;
        legend_labels{end+1}  = '$\nu^{\,\prime}$';
    
    
    if exist('v_22','var') && ~isempty(v_22)
        h = quiver(v_1(1), v_1(2), v_22(1), v_22(2), 0, '-.','Color', v_22_color, ...
            'LineWidth', 1, 'MaxHeadSize', 0.3, 'AutoScale', 'off');
        legend_handles(end+1) = h;
        legend_labels{end+1}  = '$v_c^{\,\prime\prime}$';
    end
    
    end
    if exist('v_all2','var') && ~isempty(v_all2)
        h = quiver(0, 0, v_all2(1), v_all2(2), 0, '-','Color', v_22_color, ...
            'LineWidth', 1, 'MaxHeadSize', 0.3, 'AutoScale', 'off');
        legend_handles(end+1) = h;
        legend_labels{end+1}  = '$\nu^{\,\prime\prime}$';
    end
    % 设置图例
    legend(legend_handles, legend_labels, 'Location', 'southwest', 'NumColumns', 2, 'Interpreter', 'latex', 'FontSize', 8);

    xlabel('$\nu(1)$', 'Interpreter', 'latex', 'FontSize', 8); ylabel('$\nu(2)$', 'Interpreter', 'latex', 'FontSize', 8);
    
    
    grid on;
    
   otherwise
    [K,vol1]    = convhulln(V');
    if nargin > 2
      [Klin,vol2] = convhulln(Vlin');
      ratio = vol2/vol1;
    end
      
    if k == 3
        % 初始化 legend 句柄与标签
        legend_handles = [];
        legend_labels  = {};

        % Illustrate
        if nargin > 2
            h_inv = polyplot(Klin, Vlin', 1);  % 可能返回多个 patch
            set(h_inv, 'EdgeColor', INV_edge_color, 'FaceColor', INV_fill_color);
            legend_handles(end+1) = h_inv(1);  % 只取第一个 patch 作为图例代表
            legend_labels{end+1}  = '$\Pi_{\mathrm{inv}}$';

            hold on;
            % Fix: Make V wireframe enclose Vlin
            V0 = mean(V, 2);
            V  = 1.0001 * (V - V0) + V0;

            h_ams = polyplot(K, V', 1);
            set(h_ams, 'EdgeColor', AMS_edge_color, 'FaceColor', AMS_fill_color,'FaceAlpha', 0.3);
            legend_handles(end+1) = h_ams(1);
            legend_labels{end+1}  = 'AS';
            hold off;
        else
            h_ams = polyplot(K, V', 1);
            set(h_ams, 'EdgeColor', AMS_edge_color, 'FaceColor', AMS_fill_color);
            legend_handles(end+1) = h_ams(1);
            legend_labels{end+1}  = 'AS';
        end

        % 添加 legend
        legend(legend_handles, legend_labels, 'Location', 'northeast', ...
    'Position',[0.40492232085715 0.78803228224244 0.25 0.118650618374558],...
         'Interpreter', 'latex', 'FontSize', 8);

        xlabel('$\nu(1)$', 'Interpreter', 'latex', 'FontSize', 8)
        ylabel('$\nu(2)$', 'Interpreter', 'latex', 'FontSize', 8)
        zlabel('$\nu(3)$', 'Interpreter', 'latex', 'FontSize', 8)
        view(3);
        axis equal;
        axis vis3d;
        grid on;
    end
       
end
      
function f = feasible(x,plim)
% x   m*n
% lb  m
% ub  m
  
  m = size(x,1);
  
  % Mean point
  x0 = mean(plim,2);
  
  % Make the mean point the origin
  x = x - x0*ones(1,size(x,2));
  lb = plim(:,1) - x0; % < 0
  ub = plim(:,2) - x0; % > 0
  
  % Check for feasibility
  tol = 1e-5;
  f = sum((diag(1./ub)*x <= 1+tol) & (diag(1./lb)*x <= 1+tol)) == m;

function h = polyplot(face,vert,merge)
  
  if merge 
    % Merge adjacent, parallel triangles to get fewer lines that
    % are not edges of the polyhedron.
    face4 = [];
    % Loop over all combinations of triangles
    k = 1;
    while k < size(face,1)
      l = k+1;
      while l <= size(face,1)
	iv = intersect(face(k,:),face(l,:)); % Intersecting vertices
	if length(iv) == 2 % Two common vertices
	  % Are the faces parallel?
	  niv = setxor(face(k,:),face(l,:)); % Non-intersecting vertices
	  % Vectors from first common vertex to remaining three vertices
	  A = [vert(iv(2),:)  - vert(iv(1),:);
	       vert(niv(1),:) - vert(iv(1),:);
	       vert(niv(2),:) - vert(iv(1),:)];
	  if abs(det(A))<100*eps
	    % Vectors lie in same plane -> create patch with four vertices
	    face4 = [face4 ; iv(1) niv(1) iv(2) niv(2)];
	    % ... and remove the two triangles
	    face = face([1:k-1 k+1:l-1 l+1:end],:);
	    k = k-1;
	    break
	  end	  
	end
	l = l+1;
      end % inner loop
      k = k+1;
    end % outer loop
    h = [patch('Faces',face,'Vertices',vert)
	 patch('Faces',face4,'Vertices',vert)];
  else
    % Just plot the polyhedron made up by triangles
    h = patch('Faces',face,'Vertices',vert);
  end
  
