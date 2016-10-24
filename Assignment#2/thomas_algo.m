function x = thomas_algo(a, b, c, d)
    n = length(b); u = b;
    l = [1, a];
    for ii = 2 : n
        l(ii) = a(ii - 1) / u(ii - 1);
        u(ii) = b(ii) - l(ii) * c(ii - 1);
    end
    y = d;
    for ii = 2 : n
        y(ii) = d(ii) - l(ii) * y(ii - 1);
    end
    x = y;
    x(end) = y(end) / u(end);
    for ii = n - 1: -1 : 1
        x(ii) = (y(ii) - c(ii) * x(ii + 1)) / u(ii);
    end
end