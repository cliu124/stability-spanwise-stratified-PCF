function text_value = latex_scientific(value)
% Format powers of ten without a leading 1; retain 1 and 0.1 as decimals.
if value == 0
    text_value = '0';
elseif value == 1 || abs(value-0.1) <= 10*eps(0.1)
    text_value = sprintf('%g',value);
else
    exponent = floor(log10(abs(value)));
    mantissa = value/10^exponent;
    if abs(mantissa-1) <= 10*eps(1)
        text_value = ['10^{', sprintf('%d',exponent), '}'];
    else
        text_value = [sprintf('%g',mantissa), '\times10^{', sprintf('%d',exponent), '}'];
    end
end
end
