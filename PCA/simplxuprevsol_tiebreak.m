function [y0, inB, e, itlim, errout] = simplxuprevsol_tiebreak(A, ct, b, inB, h, e, varargin)
% Bounded revised simplex with the minimum deterministic tie-break change.
%
% This file intentionally follows simplxuprevsol.m almost line by line.  The
% only behavioral change is how ties are resolved:
%   1. reduced-cost ties choose the smallest variable number in inD;
%   2. leaving-ratio ties choose the smallest current basic variable in inB.
%
% No projection, feasibility repair, alternative optimality tolerance, or
% target-specific numeric cleanup is done here.  Those changes made the solver
% hard to audit and allowed fixes for one B matrix to perturb simple 3x4 cases.

[opts, varargin] = parse_opts(varargin);
opts = default_tiebreak_opts(opts);
e = logical(e);

switch length(varargin)
    case 0
        itlim = inf;
        [m, n] = size(A);
    case 1
        itlim = varargin{1};
        [m, n] = size(A);
    case 2
        itlim = inf;
        m = varargin{1};
        n = varargin{2};
    case 3
        itlim = varargin{3};
        m = varargin{1};
        n = varargin{2};
    otherwise
        error('simplxuprevsol_tiebreak:InvalidInput', 'Unexpected number of optional arguments.');
end

% Keep the book/MATLAB reference tolerance.  C++ uses its problem.tol because
% the target arithmetic is float, but the algorithmic branches are the same.
tol = 1e-10;

% Index list for non-basic variables.
nind = 1:(n - m);
inD = setdiff(1:n, inB);

% Variables initialized at their upper bound are represented by sign flips.
A(:, ~e) = -A(:, ~e);
ct(~e) = -ct(~e);
b = b + A(:, ~e) * h(~e);

% Initial basic solution.
y0 = A(:, inB) \ b;

done = false;
unbounded = false;

while (~done && ~unbounded) && (itlim > 0)
    itlim = itlim - 1;

    % Relative costs for all non-basic variables.
    lamt = ct(inB) / A(:, inB);
    rdt = ct(inD) - lamt * A(:, inD);

    % Original solver used min(rdt).  The tie version only changes which
    % candidate is selected when several reduced costs are nearly equal.
    [minr, qind] = select_entering(rdt, inD, opts);
    if minr >= 0
        done = true;
        break;
    end

    qel = inD(qind);
    yq = A(:, inB) \ A(:, qel);

    if all(abs(yq) <= tol)
        unbounded = true;
        break;
    end

    % Ratio test copied from simplxuprevsol.m.
    rat = y0 ./ yq;
    hinB = h(inB);
    indm = yq < 0;
    rat(indm) = rat(indm) - hinB(indm) ./ yq(indm);
    rat(abs(yq) <= tol) = inf;

    % Original solver used min(rat).  The tie version only changes which
    % candidate leaves when ratios are nearly equal.
    [minrat, p] = select_leaving(rat, inB, opts);

    % Bland anti-cycling branch is preserved.  The entering variable tie is
    % deterministic here too, but the zero-step condition remains the original.
    if abs(minrat) <= tol
        [~, qind] = select_bland_negative(rdt, inD);
        qel = inD(qind);

        yq = A(:, inB) \ A(:, qel);
        if all(abs(yq) <= tol)
            unbounded = true;
            break;
        end

        rat = y0 ./ yq;
        hinB = h(inB);
        indm = yq < 0;
        rat(indm) = rat(indm) - hinB(indm) ./ yq(indm);
        rat(abs(yq) <= tol) = inf;

        [minrat, p] = select_leaving(rat, inB, opts);
    end

    % Pivot/update logic is unchanged from simplxuprevsol.m.
    if minrat >= h(qel)
        e(qel) = ~e(qel);
        A(:, qel) = -A(:, qel);
        b = b + A(:, qel) * h(qel);
        ct(qel) = -ct(qel);
    elseif yq(p) > 0
        pel = inB(p);
        inB(p) = qel;
        inD(qind) = pel;
    else
        pel = inB(p);
        e(pel) = ~e(pel);
        A(:, pel) = -A(:, pel);
        inB(p) = qel;
        inD(qind) = pel;
        ct(pel) = -ct(pel);
        b = b + A(:, pel) * h(pel);
    end

    y0 = A(:, inB) \ b;
end

errout = unbounded;
end

function [opts, args] = parse_opts(args)
opts = struct();
if isempty(args)
    return;
end

% Preferred form:
%   simplxuprevsol_tiebreak(A,c,b,inB,h,e,m,n,itlim,opts)
if isstruct(args{end})
    opts = args{end};
    args = args(1:end - 1);
    return;
end

% Backward-compatible form:
%   simplxuprevsol_tiebreak(A,c,b,inB,h,e,opts,m,n,itlim)
if isstruct(args{1})
    opts = args{1};
    args = args(2:end);
end
end

function opts = default_tiebreak_opts(opts)
% Only these two knobs remain.  They define "nearly equal" for deterministic
% pivot selection; they do not change simplex feasibility or optimality tests.
if ~isfield(opts, 'tie_rel_tol')
    opts.tie_rel_tol = 1e-5;
end
if ~isfield(opts, 'tie_abs_tol')
    opts.tie_abs_tol = 1e-6;
end
end

function [minr, qind] = select_entering(rdt, inD, opts)
[minr_raw, raw_qind] = min(rdt);
tol = tie_tol(rdt, opts);
candidates = find(abs(rdt - minr_raw) <= tol);

if isempty(candidates)
    qind = raw_qind;
else
    [~, local] = min(inD(candidates));
    qind = candidates(local);
end

minr = rdt(qind);
end

function [minrat, p] = select_leaving(rat, inB, opts)
% The leaving ratio is the feasibility step length.  Unlike reduced-cost
% tie-breaks, we must not choose a candidate with a larger ratio merely
% because it is "near" the minimum; doing that can step past another bound
% and break 0 <= x <= h.  Therefore the deterministic rule is applied only
% to exact equal ratios, matching the original minimum ratio test otherwise.
[minrat_raw, raw_p] = min(rat);

candidates = find(rat == minrat_raw);

if isempty(candidates)
    p = raw_p;
else
    [~, local] = min(inB(candidates));
    p = candidates(local);
end

minrat = rat(p);
end

function [minr, qind] = select_bland_negative(rdt, inD)
candidates = find(rdt < 0);
[~, local] = min(inD(candidates));
qind = candidates(local);
minr = rdt(qind);
end

function tol = tie_tol(x, opts)
scale = max(1, max(abs(double(x(:)))));
tol = opts.tie_abs_tol + opts.tie_rel_tol * scale;
end
