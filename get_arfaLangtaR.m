function y = get_arfaLangtaR(afa,langta,r)
y=linspace(0,0,length(afa))';
for i=1:1:length(afa)
    y = y+afa(i)*langta(i)*r(:,1);
end
end