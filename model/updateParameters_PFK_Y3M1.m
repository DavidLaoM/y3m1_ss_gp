% 153 to 160          factor on top of PFK
if ((isfield(setup,'missing_regulation_PFK_option'))&&(setup.missing_regulation_PFK_option == 1))
    p.PFK_factor1 = 1 .* 10 .^ x(153);
    p.PFK_factor2 = 1 .* 10 .^ x(154);
    p.PFK_factor3 = 1 .* 10 .^ x(155);
    p.PFK_factor4 = 1 .* 10 .^ x(156);

    p.PFK_factor5 = 1 .* 10 .^ x(157);
    p.PFK_factor6 = 1 .* 10 .^ x(158);
    p.PFK_factor7 = 1 .* 10 .^ x(159);
    p.PFK_factor8 = 1 .* 10 .^ x(160);

    setup.experiment2 = experiment; % link to then recall which factor is taken
end
