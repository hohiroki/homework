function ouputGrid(X,Y,xnum,ynum)
fid=fopen('grid4060.in','wt');
fprintf(fid,'%i',xnum-1);
fprintf(fid,'    ');
fprintf(fid,'%i',ynum-1);
fprintf(fid,'\r\n');
for i=1:1:xnum
    for j=1:1:ynum
        fprintf(fid,'%f',X(i,j));
        fprintf(fid,'    ');
        fprintf(fid,'%f\r\n',Y(i,j));
    end
end
fclose(fid);
end