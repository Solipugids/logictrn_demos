function [st] = randstring(length)

    symbols = ['a':'z' 'A':'Z' '0':'9'];       
    nums = randi(numel(symbols),[1 length]);
    st = symbols (nums);
end

