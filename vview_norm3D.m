function ratio = vview_norm3D(B, plim, P, v_cmd, U1, U2, U3, U4)
% VVIEW - View the attainable virtual control set (k = 1,2,3).
%
% 1) vview_norm(B,plim)
% 2) ratio = vview_norm(B,plim,P)
%
% When k = 3, this version also plots v_cmd (and U1..U4) in 3D if provided.
%
% Inputs:
%   B     : control effectiveness (k x m)
%   plim  : actuator pos. limits [min max] (m x 2)
%   P     : allocation (m x k) with B*P ≈ I_k (optional)
%   v_cmd : reference trajectory in v-space (3 x N)   (optional)
%   U1..U4: other trajectories in v-space (3 x N each)(optional)
%
% Output:
%   ratio : size(Π_inv)/size(AS) (length/area/volume), if P provided

  if nargin < 4
      v_cmd=[]; U1=[]; U2=[]; U3=[]; U4=[];
  end

  % ===== Colors =====
  AMS_fill_color = [0.95 0.95 1.00];
  AMS_edge_color = [0.10 0.35 0.75];

  INV_fill_color = [1.00 1.00 0.90];
  INV_edge_color = [0.85 0.33 0.10];

  inv_color = [0.85 0.10 0.40];
  d_color   = [0.00 0.00 0.00];
  pca_color = [0.00 0.60 0.60];
  cmd_color = [0.40 0.00 0.70];
  wls_color = [0.70 0.50 0.00];

  [k,m] = size(B);

  % ---------- (a) Max attainable set: V = B*U, u in box(plim) ----------
  % corners of actuator box
  idx = zeros(2^m,m);
  M = 1:m;
  for i = 1:2^m
      cbin = dec2bin(i-1,m);
      c = str2num(cbin')'; %#ok<ST2NM>
      c = c(end:-1:1);
      idx(i,:) = 2*M - c;
  end
  plimT = plim';
  U = plimT(idx)';          % m x (2^m)
  V = B*U;                  % k x (2^m)

  % ---------- (b) Linear allocation set: v s.t. u = P v in box ----------
  if nargin > 2 && ~isempty(P)
      sub_idx = idx(1:2^k,1:k);
      Ulin = [];
      i_dim = nchoosek(1:m,k);
      for i = 1:size(i_dim,1)
          sub_plimT = plimT(:,i_dim(i,:));
          sub_u_boundary = sub_plimT(sub_idx)'; % k x (2^k)
          sub_P = P(i_dim(i,:),:);             % k x k
          if rank(sub_P) == k
              v = sub_P \ sub_u_boundary;      % k x (2^k)
              u_boundary = P * v;              % m x (2^k)
              i_feas = feasible(u_boundary,plim);
              Ulin = [Ulin u_boundary(:,i_feas)]; %#ok<AGROW>
          end
      end
      Vlin = B*Ulin; % k x n_lin
  end

  % ---------- Plot ----------
  clf
  switch k
      case 1
          K = [min(V) max(V)];
          if exist('Vlin','var') && ~isempty(Vlin)
              Klin = [min(Vlin) max(Vlin)];
              ratio = diff(Klin)/diff(K);
              plot(K,[0 0],'b-o',Klin,-[0 0],'r-o');
          else
              plot(K,[0 0],'b-o');
          end
          xlabel('v');

      case 2
          hold on; axis equal;
          legend_handles = []; legend_labels = {};
          [K,area1] = convhull(V(1,:),V(2,:));
          h_fill1 = fill(V(1,K),V(2,K),AMS_fill_color,'FaceAlpha',0.1,...
              'EdgeColor',AMS_edge_color,'LineWidth',0.1);
          legend_handles(end+1)=h_fill1; legend_labels{end+1}='AS';

          if exist('Vlin','var') && ~isempty(Vlin)
              [Klin,area2] = convhull(Vlin(1,:),Vlin(2,:));
              h_fill2 = fill(Vlin(1,Klin),Vlin(2,Klin),INV_fill_color,...
                  'EdgeColor',INV_edge_color,'LineWidth',0.1,'FaceAlpha',0.1);
              legend_handles(end+1)=h_fill2; legend_labels{end+1}='$\Pi_{\mathrm{inv}}$';
              ratio = area2/area1;
          end

          % Optional command trajectories in 2D (kept for completeness)
          if ~isempty(v_cmd)
              h=plot(v_cmd(1,:),v_cmd(2,:),':','Color',cmd_color,'LineWidth',1.8);
              legend_handles(end+1)=h; legend_labels{end+1}='$\nu$';
          end
          if ~isempty(U1)
              h=plot(U1(1,:),U1(2,:),'--','Color',inv_color,'LineWidth',1.8);
              legend_handles(end+1)=h; legend_labels{end+1}='$\nu_{inv}$';
          end
          if ~isempty(U2)
              h=plot(U2(1,:),U2(2,:),'-.' ,'Color',d_color, 'LineWidth',2);
              legend_handles(end+1)=h; legend_labels{end+1}='$\nu_{d}$';
          end
          if ~isempty(U3)
              h=plot(U3(1,:),U3(2,:),'-'  ,'Color',pca_color,'LineWidth',1.2);
              legend_handles(end+1)=h; legend_labels{end+1}='$\nu_{pca}$';
          end
          if ~isempty(U4)
              h=plot(U4(1,:),U4(2,:),'-'  ,'Color',wls_color,'LineWidth',1.2);
              legend_handles(end+1)=h; legend_labels{end+1}='$\nu_{wls}$';
          end

          legend(legend_handles,legend_labels,'Location','northeast',...
              'NumColumns',1,'Interpreter','latex','FontSize',9);
          xlabel('$\nu(1)$','Interpreter','latex','FontSize',9);
          ylabel('$\nu(2)$','Interpreter','latex','FontSize',9);
          grid on;

      otherwise
          % ---- k == 3 branch with 3D plots + v_cmd/U1..U4 ----
          [K,vol1] = convhulln(V');  % faces for AS
          if exist('Vlin','var') && ~isempty(Vlin)
              [Klin,vol2] = convhulln(Vlin');
              ratio = vol2/vol1;
          end

          hold on; axis equal; grid on; %axis vis3d; 
          legend_handles = []; legend_labels = {};

          if exist('Vlin','var') && ~isempty(Vlin)
              h_inv = polyplot(Klin, Vlin', 1);
              set(h_inv,'EdgeColor',INV_edge_color,'FaceColor',INV_fill_color,'FaceAlpha',0.35);
              legend_handles(end+1) = h_inv(1);  % take one handle
              legend_labels{end+1}  = '$\Pi_{\mathrm{inv}}$';

              % make AS slightly inflated so it encloses Π_inv visually
              V0 = repmat(mean(V,2),1,size(V,2));
              V  = 1.0001*(V - V0) + V0;
              h_as = polyplot(K, V', 1);
              set(h_as,'EdgeColor',AMS_edge_color,'FaceColor','none','LineWidth',1.2);
              legend_handles(end+1) = h_as(1);
              legend_labels{end+1}  = 'AS';
          else
              h_as = polyplot(K, V', 1);
              set(h_as,'EdgeColor',AMS_edge_color,'FaceColor',AMS_fill_color,'FaceAlpha',0.35);
              legend_handles(end+1) = h_as(1);
              legend_labels{end+1}  = 'AS';
          end

          % ---- plot command/estimator trajectories in 3D if provided ----
          % expected size: 3 x N for each
          if ~isempty(v_cmd)
              h = plot3(v_cmd(1,:),v_cmd(2,:),v_cmd(3,:),':','Color',cmd_color,'LineWidth',2);
              legend_handles(end+1)=h; legend_labels{end+1}='$\nu$';
          end
          if ~isempty(U1)
              h = plot3(U1(1,:),U1(2,:),U1(3,:),'--','Color',inv_color,'LineWidth',2);
              legend_handles(end+1)=h; legend_labels{end+1}='$\nu_{inv}$';
          end
          if ~isempty(U2)
              h = plot3(U2(1,:),U2(2,:),U2(3,:),'-.' ,'Color',d_color,'LineWidth',2);
              legend_handles(end+1)=h; legend_labels{end+1}='$\nu_{d}$';
          end
          if ~isempty(U3)
              h = plot3(U3(1,:),U3(2,:),U3(3,:),'-'  ,'Color',pca_color,'LineWidth',1.4);
              legend_handles(end+1)=h; legend_labels{end+1}='$\nu_{pca}$';
          end
          if ~isempty(U4)
              h = plot3(U4(1,:),U4(2,:),U4(3,:),'-'  ,'Color',wls_color,'LineWidth',1.4);
              legend_handles(end+1)=h; legend_labels{end+1}='$\nu_{wls}$';
          end

          xlabel('$\nu(1)$','Interpreter','latex','FontSize',9);
          ylabel('$\nu(2)$','Interpreter','latex','FontSize',9);
          zlabel('$\nu(3)$','Interpreter','latex','FontSize',9);
          view(135,25);

          if ~isempty(legend_handles)
              legend(legend_handles,legend_labels,'Location','northeast',...
                  'NumColumns',1,'Interpreter','latex','FontSize',9);
          end
  end
end

% ---------- helpers ----------
function f = feasible(x,plim)
  % x: m x n
  m = size(x,1);
  x0 = mean(plim,2);
  x  = x - x0*ones(1,size(x,2));
  lb = plim(:,1) - x0; % < 0
  ub = plim(:,2) - x0; % > 0
  tol = 1e-5;
  f = sum((diag(1./ub)*x <= 1+tol) & (diag(1./lb)*x <= 1+tol)) == m;
end

function h = polyplot(face,vert,merge)
  if merge
      face4 = [];
      k = 1;
      while k < size(face,1)
          l = k+1;
          while l <= size(face,1)
              iv = intersect(face(k,:),face(l,:));
              if length(iv)==2
                  niv = setxor(face(k,:),face(l,:));
                  A = [vert(iv(2),:)  - vert(iv(1),:);
                       vert(niv(1),:) - vert(iv(1),:);
                       vert(niv(2),:) - vert(iv(1),:)];
                  if abs(det(A))<100*eps
                      face4 = [face4 ; iv(1) niv(1) iv(2) niv(2)]; %#ok<AGROW>
                      face = face([1:k-1 k+1:l-1 l+1:end],:);
                      k = k-1;
                      break
                  end
              end
              l = l+1;
          end
          k = k+1;
      end
      h = [patch('Faces',face ,'Vertices',vert)
           patch('Faces',face4,'Vertices',vert)];
  else
      h = patch('Faces',face,'Vertices',vert);
  end
end