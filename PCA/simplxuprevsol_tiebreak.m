function [y0, inB, e, itlim, errout] = simplxuprevsol_tiebreak(A, ct, b, inB, h, e, opts, varargin)
% Bounded revised simplex with deterministic tie-break rules.
%
% This keeps the structure of simplxuprevsol.m, but replaces ambiguous
% near-equality choices with explicit rules so MATLAB and C++ are less likely
% to diverge because of roundoff.
%
% Default tie-break settings match alloc_cpp/src/ControlAllocation:
%   opts.tie_rel_tol      = 1e-5
%   opts.tie_abs_tol      = 1e-6
%   opts.zero_tie_abs_tol = 3e-5
%
% Rules:
%   1. Near-equal reduced costs choose the smallest variable index.
%   2. Near-equal ratios choose the smallest current basic variable index.
%   3. Near equality in minrat vs h(qel) takes the leaving-variable branch.

if nargin < 7 || isempty(opts)
    opts = struct();
end
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

tol = 1e-10;
inD = setdiff(1:n, inB);

% Adjust signs if variables are initialized at upper bounds.
A(:, ~e) = -A(:, ~e);
ct(~e) = -ct(~e);
b = b + A(:, ~e) * h(~e);

y0 = A(:, inB) \ b;
done = false;
unbounded = false;

while (~done || ~unbounded) && (itlim > 0)
    itlim = itlim - 1;

    lamt = ct(inB) / A(:, inB);
    rdt = ct(inD) - lamt * A(:, inD);

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

    rat = compute_ratio(y0, yq, h(inB), tol);
    [minrat, p] = select_leaving(rat, inB, opts);

    if abs(minrat) <= max(tol, opts.zero_tie_abs_tol)
        [~, qind] = select_bland_negative(rdt, inD);
        qel = inD(qind);

        yq = A(:, inB) \ A(:, qel);
        if all(abs(yq) <= tol)
            unbounded = true;
            break;
        end

        rat = compute_ratio(y0, yq, h(inB), tol);
        [minrat, p] = select_leaving(rat, inB, opts);
    end

    bound_tol = compare_tol(minrat, h(qel), opts);
    if minrat > h(qel) + bound_tol
        % Entering variable goes to its opposite bound; basis is unchanged.
        e(qel) = ~e(qel);
        A(:, qel) = -A(:, qel);
        b = b + A(:, qel) * h(qel);
        ct(qel) = -ct(qel);
    elseif yq(p) > 0
        % Leaving variable returns to lower bound.
        pel = inB(p);
        inB(p) = qel;
        inD(qind) = pel;
    else
        % Leaving variable moves to upper bound.
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

function opts = default_tiebreak_opts(opts)
    if ~isfield(opts, 'tie_rel_tol')
        opts.tie_rel_tol = 1e-5;
    end
    if ~isfield(opts, 'tie_abs_tol')
        opts.tie_abs_tol = 1e-6;
    end
    if ~isfield(opts, 'zero_tie_abs_tol')
        opts.zero_tie_abs_tol = 3e-5;
    end
end

function rat = compute_ratio(y0, yq, hinB, tol)
    rat = y0 ./ yq;
    indm = yq < 0;
    rat(indm) = rat(indm) - hinB(indm) ./ yq(indm);
    rat(abs(yq) <= tol) = inf;
end

function [minr, qind] = select_entering(rdt, inD, opts)
    minr_raw = min(rdt);
    candidates = find(abs(rdt - minr_raw) <= vector_tol(rdt, opts));
    [~, local] = min(inD(candidates));
    qind = candidates(local);
    minr = rdt(qind);
end

function [minr, qind] = select_bland_negative(rdt, inD)
    candidates = find(rdt < 0);
    [~, local] = min(inD(candidates));
    qind = candidates(local);
    minr = rdt(qind);
end

function [minrat, p] = select_leaving(rat, inB, opts)
    minrat_raw = min(rat);
    finite_rat = rat(isfinite(rat));
    tol = vector_tol(finite_rat, opts);

    if abs(minrat_raw) <= opts.zero_tie_abs_tol
        tol = max(tol, opts.zero_tie_abs_tol);
        candidates = find(abs(rat) <= tol);
    else
        candidates = find(abs(rat - minrat_raw) <= tol);
    end

    [~, local] = min(inB(candidates));
    p = candidates(local);
    minrat = rat(p);
end

function tol = vector_tol(x, opts)
    if isempty(x)
        scale = 1;
    else
        scale = max(1, max(abs(x(:))));
    end
    tol = opts.tie_abs_tol + opts.tie_rel_tol * scale;
end

function tol = compare_tol(a, b, opts)
    scale = max(1, max(abs([a, b])));
    tol = opts.tie_abs_tol + opts.tie_rel_tol * scale;
end
