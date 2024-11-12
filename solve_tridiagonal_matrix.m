function result = solve_tridiagonal_matrix(a, b, c, d)
%% 函数说明：
%   输入参数:
%       a -->> 下三角
%       b -->> 主对角
%       c -->> 上三角
%       d -->> 右边项
%   输出参数:
%       result -->> 方程组的解向量
%% 方程求解：
    % 用于控制循环
    n = max(size(b));

    % 初始化向量
    u = zeros(1, n);
    y = zeros(1, n);
    l = zeros(1, n);
    r = zeros(1, n);

    % 初始化 u(1) 和 y(1)
    u(1) = b(1);
    y(1) = d(1);

    % 向前消元（追的过程）
    for i = 2 : n
        l(i) = a(i) / u(i - 1);
        u(i) = b(i) - l(i) * c(i - 1);
        y(i) = d(i) - l(i) * y(i - 1);
    end

    % 回代求解（赶的过程）
    r(n) = 0;
    for i = (n - 1) : -1 : 1
        r(i) = (y(i) - c(i) * r(i + 1)) / u(i);
    end

    result = r;
end