function [y0, inB, e, itlim, errout] = simplxuprevsol_tiebreak(A, ct, b, inB, h, e, varargin)
% Bounded revised simplex with deterministic, guarded tie-break rules.
%
% Daily use should not tune simplex internals. This solver keeps only one
% useful behavior: deterministic tie-breaks. Entering-variable ties use the
% historical control-scale tolerances, while leaving-ratio ties use an
% internal single-precision eps scale because ratio tests protect feasibility.
%
% Rules:
%   1. Near-equal reduced costs choose the smallest variable index.
%   2. Near-equal leaving ratios choose the smallest current basic variable
%      index only when the ratios are indistinguishable at float precision.
%   3. If that guarded tie-break produces an infeasible pivot, fall back to
%      the raw minimum-ratio pivot.
%
% Parameter guide, matching alloc_cpp/src/ControlAllocation/ControlAllocation.h:
%   opts.machine_eps      Main hardware knob. Use eps('single') for MCU float
%                         builds, eps('double') only when the complete solver
%                         arithmetic is double, or a measured larger value to
%                         simulate a noisier platform in MATLAB.
%   opts.tie_abs_tol,
%   opts.tie_rel_tol      Entering reduced-cost tie window. These are
%                         control-scale tolerances; increasing them makes more
%                         near-equal reduced costs choose the smaller variable
%                         index, but too large a value can select a different
%                         optimal vertex.
%   leaving/bound/
%   feasibility/pivot     These are derived from machine_eps by default. Keep
%                         them at defaults unless diagnosing target numerics:
%                         leaving-ratio and feasibility decisions should be set
%                         by floating-point precision, not command magnitude.

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

tol = opts.pivot_zero_tol;
inD = setdiff(1:n, inB);

% Variables initialized at their upper bound are represented by sign flips.
A(:, ~e) = -A(:, ~e);
ct(~e) = -ct(~e);
b = b + A(:, ~e) * h(~e);

y0 = A(:, inB) \ b;
done = false;
unbounded = false;

while ~done && ~unbounded && (itlim > 0)
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
    [minrat, p, raw_p] = select_leaving(rat, inB, opts);

    if abs(minrat) <= opts.zero_step_tol
        [~, qind] = select_bland_negative(rdt, inD);
        qel = inD(qind);

        yq = A(:, inB) \ A(:, qel);
        if all(abs(yq) <= tol)
            unbounded = true;
            break;
        end

        rat = compute_ratio(y0, yq, h(inB), tol);
        [minrat, p, raw_p] = select_leaving(rat, inB, opts);
    end

    A_before = A;
    b_before = b;
    ct_before = ct;
    e_before = e;
    inB_before = inB;
    inD_before = inD;

    bound_tol = bound_compare_tol(minrat, h(qel), opts);
    [A, b, ct, e, inB, inD] = apply_pivot(A, b, ct, e, inB, inD, qind, qel, yq, p, h, minrat, bound_tol);
    y0 = A(:, inB) \ b;

    if p ~= raw_p && violates_box(y0, h(inB), opts)
        A = A_before;
        b = b_before;
        ct = ct_before;
        e = e_before;
        inB = inB_before;
        inD = inD_before;

        p = raw_p;
        minrat = rat(p);
        bound_tol = bound_compare_tol(minrat, h(qel), opts);
        [A, b, ct, e, inB, inD] = apply_pivot(A, b, ct, e, inB, inD, qind, qel, yq, p, h, minrat, bound_tol);
        y0 = A(:, inB) \ b;
    end

    if violates_box(y0, h(inB), opts)
        unbounded = true;
        break;
    end
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
        args = args(1:end-1);
        if isempty(opts)
            opts = struct();
        end
        return;
    end

    % Backward-compatible form used by older local code:
    %   simplxuprevsol_tiebreak(A,c,b,inB,h,e,opts,m,n,itlim)
    if isstruct(args{1})
        opts = args{1};
        args = args(2:end);
        if isempty(opts)
            opts = struct();
        end
    end
end

function opts = default_tiebreak_opts(opts)
    % These historical knobs only affect entering-variable ties. Daily tests
    % use the defaults; leaving-ratio protection is fixed below.
    if ~isfield(opts, 'machine_eps')
        opts.machine_eps = eps('single');
    end
    if ~isfinite(opts.machine_eps) || opts.machine_eps <= 0
        opts.machine_eps = eps('single');
    end
    if ~isfield(opts, 'tie_rel_tol')
        opts.tie_rel_tol = 1e-5;
    end
    if ~isfield(opts, 'tie_abs_tol')
        opts.tie_abs_tol = 1e-6;
    end
    numeric_eps = double(opts.machine_eps);
    if ~isfield(opts, 'leaving_tie_abs_tol')
        opts.leaving_tie_abs_tol = 16 * numeric_eps;
    end
    if ~isfield(opts, 'leaving_tie_rel_tol')
        opts.leaving_tie_rel_tol = 16 * numeric_eps;
    end
    if ~isfield(opts, 'zero_step_tol')
        opts.zero_step_tol = 16 * numeric_eps;
    end
    if ~isfield(opts, 'bound_abs_tol')
        opts.bound_abs_tol = 16 * numeric_eps;
    end
    if ~isfield(opts, 'bound_rel_tol')
        opts.bound_rel_tol = 16 * numeric_eps;
    end
    if ~isfield(opts, 'feasibility_abs_tol')
        opts.feasibility_abs_tol = 128 * numeric_eps;
    end
    if ~isfield(opts, 'feasibility_rel_tol')
        opts.feasibility_rel_tol = 128 * numeric_eps;
    end
    if ~isfield(opts, 'pivot_zero_tol')
        opts.pivot_zero_tol = max(1e-10, 8 * numeric_eps);
    end
end

function [A, b, ct, e, inB, inD] = apply_pivot(A, b, ct, e, inB, inD, qind, qel, yq, p, h, minrat, bound_tol)
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

function [minrat, p, raw_p] = select_leaving(rat, inB, opts)
    [minrat_raw, raw_p] = min(rat);

    if isfinite(minrat_raw)
        tol = ratio_tol(minrat_raw, opts);
        candidates = find(abs(rat - minrat_raw) <= tol);
    else
        candidates = find(rat == minrat_raw);
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

function tol = ratio_tol(minrat, opts)
    scale = max(1, abs(minrat));
    tol = opts.leaving_tie_abs_tol + opts.leaving_tie_rel_tol * scale;
end

function tol = bound_compare_tol(a, b, opts)
    scale = max(1, max(abs([a, b])));
    tol = opts.bound_abs_tol + opts.bound_rel_tol * scale;
end

function bad = violates_box(y, h, opts)
    if isempty(y)
        bad = false;
        return;
    end

    violation = max([0; -double(y(:)); double(y(:) - h(:))]);
    scale = max(1, max(abs(double(h(:)))));
    tol = opts.feasibility_abs_tol + opts.feasibility_rel_tol * scale;
    bad = violation > tol;
end
