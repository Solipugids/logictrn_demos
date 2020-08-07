
function str = fun_x_range_calc(n)

code = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'};

if n <=26
    str = code(n);
    return;
else
    n1 = floor(n / 26);
    n2 = mod(n, 26);
    
    str = strcat(code(n1), code(n2));
    
    return;
end

